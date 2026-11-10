#include "metrics.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static inline double gap_at(const Sim *S, int i, int k)
{
    return S->lx[sim_log_index(S, i - 1, k)] - S->lx[sim_log_index(S, i, k)] - S->cfg.veh_length;
}

void metrics_behaviour(const Sim *S, Behaviour *b)
{
    const SimConfig *c = &S->cfg;
    const int N = S->N, n = S->nsteps;
    double min_gap = INFINITY, max_err = 0.0, sum_sq = 0.0, peak_first = 0.0, peak_last = 0.0;
    long   unsafe = 0, count = 0;
    int    collision = 0;

    for (int i = 1; i < N; ++i) {
        double peak = 0.0;
        for (int k = 0; k <= n; ++k) {
            const double gap = gap_at(S, i, k);
            const double e   = gap - (c->d0 + c->headway * S->lv[sim_log_index(S, i, k)]);
            if (gap < min_gap) min_gap = gap;
            if (gap < c->safe_gap) ++unsafe;
            if (gap <= 0.0) collision = 1;
            if (fabs(e) > peak) peak = fabs(e);
            sum_sq += e * e;
            ++count;
        }
        if (peak > max_err) max_err = peak;
        if (i == 1)     peak_first = peak;
        if (i == N - 1) peak_last  = peak;
    }

    b->min_gap     = (N > 1) ? min_gap : NAN;
    b->max_abs_err = max_err;
    b->rms_err     = count ? sqrt(sum_sq / count) : 0.0;
    b->string_amp  = (N > 2 && peak_first > 0.0) ? peak_last / peak_first : NAN;
    b->unsafe_time = unsafe * c->dt;
    b->collision   = collision;
}

void metrics_deviation(const Sim *ref, const Sim *run, Deviation *d)
{
    const size_t nlog = (size_t)ref->N * (size_t)(ref->nsteps + 1);
    d->bit_identical = memcmp(ref->lx, run->lx, nlog * sizeof(double)) == 0 &&
                       memcmp(ref->lv, run->lv, nlog * sizeof(double)) == 0 &&
                       memcmp(ref->la, run->la, nlog * sizeof(double)) == 0;

    double mx = 0, mv = 0, ma = 0, mg = 0, sx = 0;
    for (size_t k = 0; k < nlog; ++k) {
        const double dx = fabs(run->lx[k] - ref->lx[k]);
        const double dv = fabs(run->lv[k] - ref->lv[k]);
        const double da = fabs(run->la[k] - ref->la[k]);
        if (dx > mx) mx = dx;
        if (dv > mv) mv = dv;
        if (da > ma) ma = da;
        sx += dx * dx;
    }
    for (int i = 1; i < ref->N; ++i)
        for (int k = 0; k <= ref->nsteps; ++k) {
            const double dg = fabs(gap_at(run, i, k) - gap_at(ref, i, k));
            if (dg > mg) mg = dg;
        }
    d->max_dx = mx; d->max_dv = mv; d->max_da = ma; d->max_dgap = mg;
    d->rms_dx = sqrt(sx / (double)nlog);
}

int write_trajectory_csv(const Sim *S, const char *path, int every)
{
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    if (every < 1) every = 1;
    const SimConfig *c = &S->cfg;
    fprintf(f, "t,vehicle,x,v,a,gap,spacing_error\n");
    for (int k = 0; k <= S->nsteps; k += every) {
        const double t = k * c->dt;
        for (int i = 0; i < S->N; ++i) {
            const size_t q = sim_log_index(S, i, k);
            if (i == 0) {
                fprintf(f, "%.3f,%d,%.17g,%.17g,%.17g,,\n", t, i, S->lx[q], S->lv[q], S->la[q]);
            } else {
                const double gap = gap_at(S, i, k);
                const double e   = gap - (c->d0 + c->headway * S->lv[q]);
                fprintf(f, "%.3f,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                        t, i, S->lx[q], S->lv[q], S->la[q], gap, e);
            }
        }
    }
    return fclose(f);
}
