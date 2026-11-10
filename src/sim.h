/*
 * sim.h - simulation engine: sequential baseline and OpenMP parallel version.
 *
 * Shared simulation state
 * -----------------------
 * Each vehicle publishes (x, v, a) at the end of every control step into a
 * ring buffer `ring[step % ring_len][vehicle]`. Other vehicles only read
 * published entries (radar: predecessor position/speed, V2V: predecessor
 * acceleration, possibly delayed). A vehicle's full dynamic state is private
 * to the vehicle, so no two threads ever write the same memory.
 *
 * Synchronisation rule (from the research design)
 * -----------------------------------------------
 * All vehicles read a consistent snapshot published at the last
 * synchronisation point; the next window starts only after every vehicle has
 * finished (barrier).
 *   K = 1 (exact):   barrier every control step, the snapshot is always the
 *                    previous step -> identical to the sequential baseline.
 *   K > 1 (relaxed): barrier every K steps; between barriers each vehicle
 *                    extrapolates its predecessor from the last snapshot
 *                    (constant acceleration). Fewer barriers, approximate data.
 */
#ifndef SIM_H
#define SIM_H

#include <stddef.h>
#include "config.h"
#include "vehicle.h"

typedef struct { double x, v, a; } PubData;

typedef struct {
    PubData d;
    char    pad[CACHE_LINE - sizeof(PubData)];
} PubSlot;

typedef struct {
    _Alignas(CACHE_LINE) VehParams p;
    DynState s;
    double   a_meas;   /* acceleration measured over the last control step */
    double   a_cmd;    /* jerk-limited acceleration command                */
    double   integ;    /* low-level integral state                         */
} Vehicle;

typedef struct {
    SimConfig cfg;
    int       N;          /* vehicles                                   */
    int       nsteps;     /* control steps                              */
    int       K;          /* sync window                                */
    int       D;          /* V2V delay in steps                         */
    int       ring_len;   /* >= K + D + 2                               */
    Vehicle  *veh;        /* [N]                                        */
    PubSlot  *ring;       /* [ring_len * N]                             */
    double   *lx, *lv, *la; /* trajectory log, vehicle-major [N * (nsteps+1)] */
} Sim;

int    sim_init(Sim *S, const SimConfig *cfg);
void   sim_free(Sim *S);

/* Run the whole simulation; returns wall-clock seconds. */
double sim_run_seq(Sim *S);
double sim_run_omp(Sim *S, int threads);

static inline size_t sim_log_index(const Sim *S, int vehicle, int step)
{
    return (size_t)vehicle * (size_t)(S->nsteps + 1) + (size_t)step;
}

#endif
