#include "sim.h"
#include "controller.h"
#include "scenario.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/* ------------------------------------------------------------------------ */
/* Counter-based random numbers: the value depends only on (seed, vehicle,  */
/* step), never on which thread runs the vehicle or in which order, so      */
/* sensor noise is identical in sequential and parallel runs.               */
/* ------------------------------------------------------------------------ */
static inline uint64_t splitmix64(uint64_t z)
{
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline double uniform_hash(uint64_t seed, uint64_t a, uint64_t b)
{
    uint64_t h = splitmix64(seed ^ splitmix64(a * 0x100000001B3ULL ^ splitmix64(b)));
    return ((double)(h >> 11) + 0.5) * 0x1.0p-53;          /* (0,1) */
}

static inline double gauss_hash(uint64_t seed, uint64_t a, uint64_t b)
{
    const double u1 = uniform_hash(seed, a, b);
    const double u2 = uniform_hash(seed ^ 0xA5A5A5A5ULL, a, b);
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

static void *aligned_calloc(size_t n_bytes)
{
    size_t sz = (n_bytes + CACHE_LINE - 1) / CACHE_LINE * CACHE_LINE;
    void *p = aligned_alloc(CACHE_LINE, sz);
    if (p) memset(p, 0, sz);
    return p;
}

static inline const PubData *pub_r(const Sim *S, int step, int i)
{
    return &S->ring[(size_t)(step % S->ring_len) * S->N + i].d;
}

static inline PubData *pub_w(Sim *S, int step, int i)
{
    return &S->ring[(size_t)(step % S->ring_len) * S->N + i].d;
}

/* ------------------------------------------------------------------------ */

int sim_init(Sim *S, const SimConfig *cfg)
{
    memset(S, 0, sizeof *S);
    S->cfg      = *cfg;
    S->N        = cfg->n_vehicles;
    S->nsteps   = (int)lround(cfg->duration / cfg->dt);
    S->K        = cfg->sync_window < 1 ? 1 : cfg->sync_window;
    S->D        = (int)lround(cfg->delay_s / cfg->dt);
    S->ring_len = S->K + S->D + 2;

    const size_t nlog = (size_t)S->N * (size_t)(S->nsteps + 1);
    S->veh  = aligned_calloc(sizeof(Vehicle) * S->N);
    S->ring = aligned_calloc(sizeof(PubSlot) * S->N * S->ring_len);
    S->lx   = aligned_calloc(sizeof(double) * nlog);
    S->lv   = aligned_calloc(sizeof(double) * nlog);
    S->la   = aligned_calloc(sizeof(double) * nlog);
    if (!S->veh || !S->ring || !S->lx || !S->lv || !S->la) { sim_free(S); return -1; }

    /* initial conditions: equilibrium cruise at v0, desired spacing */
    const double v0  = cfg->v0;
    const double gap0 = cfg->d0 + cfg->headway * v0;
    for (int i = 0; i < S->N; ++i) {
        Vehicle *V = &S->veh[i];
        double mass = 1500.0;
        if (cfg->hetero) mass = 1200.0 + 1200.0 * uniform_hash(cfg->seed, 0xBEEF, (uint64_t)i);
        vehicle_params_init(&V->p, mass);

        const double drag = (i == 0) ? 1.0 : drafting_factor(gap0);
        V->s.x  = -(double)i * (gap0 + cfg->veh_length);
        V->s.v  = v0;
        V->s.wf = v0 / V->p.R;
        V->s.wr = v0 / V->p.R;
        V->s.T  = V->p.R * vehicle_resistance(&V->p, V->s.x, v0, drag, 0.0, 1);
        V->a_meas = V->a_cmd = V->integ = 0.0;

        PubData *o = pub_w(S, 0, i);
        o->x = V->s.x; o->v = v0; o->a = 0.0;
        const size_t k = sim_log_index(S, i, 0);
        S->lx[k] = V->s.x; S->lv[k] = v0; S->la[k] = 0.0;
    }
    return 0;
}

void sim_free(Sim *S)
{
    free(S->veh); free(S->ring); free(S->lx); free(S->lv); free(S->la);
    S->veh = NULL; S->ring = NULL; S->lx = S->lv = S->la = NULL;
}

/*
 * Advance vehicle i from step s to s+1. Only reads published data from steps
 * <= s0 (the last synchronisation point) and writes only vehicle i's own
 * memory, so vehicles inside one window can run in any order or in parallel.
 */
static void vehicle_advance(Sim *S, int i, int s, int s0)
{
    const SimConfig *c = &S->cfg;
    Vehicle *V = &S->veh[i];
    const double dt = c->dt;
    const double t  = s * dt;
    const double x  = V->s.x, v = V->s.v;
    double a_des, drag = 1.0;

    if (i == 0) {
        a_des = leader_accel(c, t, v);
    } else {
        const PubData *P = pub_r(S, s0, i - 1);
        double xp = P->x, vp = P->v;
        if (s != s0) {                       /* relaxed: extrapolate snapshot */
            const double tau = (s - s0) * dt;
            xp = P->x + P->v * tau + 0.5 * P->a * tau * tau;
            vp = P->v + P->a * tau;
        }
        /* V2V: predecessor acceleration from D steps ago (if already synced) */
        int r = s - S->D;
        if (r < 0) r = 0;
        const double ap = (r <= s0) ? pub_r(S, r, i - 1)->a : P->a;

        double gap = xp - x - c->veh_length;
        if (c->gap_noise_std > 0.0)
            gap += c->gap_noise_std * gauss_hash(c->seed, (uint64_t)i, (uint64_t)s);

        drag  = drafting_factor(gap);
        a_des = cacc_accel(c, gap, v, vp, V->a_meas, ap);
    }

    V->a_cmd = shape_accel(c, a_des, V->a_cmd);
    const double T_cmd = lowlevel_torque(c, &V->p, &V->integ, x, v, V->a_cmd, V->a_meas, drag);
    vehicle_integrate(&V->p, &V->s, T_cmd, drag, V->a_meas, t, dt, c->substeps);
    V->a_meas = (V->s.v - v) / dt;

    PubData *o = pub_w(S, s + 1, i);
    o->x = V->s.x; o->v = V->s.v; o->a = V->a_meas;

    const size_t k = sim_log_index(S, i, s + 1);
    S->lx[k] = V->s.x; S->lv[k] = V->s.v; S->la[k] = V->a_meas;
}

/* ---- Sequential baseline: one core, vehicles updated one after another -- */
double sim_run_seq(Sim *S)
{
    const int N = S->N, K = S->K, nsteps = S->nsteps;
    const double t0 = omp_get_wtime();
    for (int s0 = 0; s0 < nsteps; s0 += K) {
        const int s1 = (s0 + K < nsteps) ? s0 + K : nsteps;
        for (int i = 0; i < N; ++i)
            for (int s = s0; s < s1; ++s)
                vehicle_advance(S, i, s, s0);
    }
    return omp_get_wtime() - t0;
}

/*
 * ---- Parallel version (OpenMP) -------------------------------------------
 * Master/worker mapping: the team is created once (fork); in every window the
 * vehicles are distributed over the threads (`omp for`), and the implicit
 * barrier at the end of `omp for` is the per-window synchronisation in which
 * all results are collected before the next window starts.
 */
double sim_run_omp(Sim *S, int threads)
{
    const int N = S->N, K = S->K, nsteps = S->nsteps;
    const double t0 = omp_get_wtime();
    #pragma omp parallel num_threads(threads) default(none) shared(S) firstprivate(N, K, nsteps)
    {
        for (int s0 = 0; s0 < nsteps; s0 += K) {
            const int s1 = (s0 + K < nsteps) ? s0 + K : nsteps;
            #pragma omp for schedule(runtime)
            for (int i = 0; i < N; ++i)
                for (int s = s0; s < s1; ++s)
                    vehicle_advance(S, i, s, s0);
            /* implicit barrier = synchronisation point */
        }
    }
    return omp_get_wtime() - t0;
}
