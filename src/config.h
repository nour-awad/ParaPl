/*
 * config.h - experiment configuration shared by every module.
 *
 * Everything listed under "HELD CONSTANT" in the research design lives here,
 * so the sequential baseline and every parallel run are guaranteed to use the
 * same vehicle model, controller, scenario, time step and duration. Only the
 * execution method (sequential / OpenMP threads / sync window) changes.
 */
#ifndef CONFIG_H
#define CONFIG_H

#include <stdint.h>

#define G_ACC   9.81    /* gravity [m/s^2]      */
#define RHO_AIR 1.225   /* air density [kg/m^3] */

/* Apple M-series cores use 128-byte cache lines (x86 uses 64). Padding
 * per-vehicle data to 128 bytes prevents false sharing on both. */
#define CACHE_LINE 128

typedef enum { SCEN_CRUISE = 0, SCEN_BRAKE = 1 } Scenario;

typedef struct {
    /* ---- platoon & scenario ------------------------------------------- */
    int      n_vehicles;   /* leader + followers                          */
    Scenario scenario;     /* steady cruise or leader braking              */
    double   duration;     /* simulated time [s]                           */
    double   dt;           /* control / V2V step [s]                       */
    int      substeps;     /* RK4 dynamics sub-steps per control step      */

    /* ---- execution / synchronisation ---------------------------------- */
    int      sync_window;  /* K = 1: exact per-step sync (barrier every step)
                              K > 1: relaxed sync, barrier every K steps   */

    /* ---- spacing policy & CACC controller ----------------------------- */
    double v0;             /* initial / cruise speed [m/s]                 */
    double d0;             /* standstill distance [m]                      */
    double headway;        /* constant time gap h [s]                      */
    double kp, kd, ka;     /* CACC gains                                   */
    double veh_length;     /* vehicle length [m]                           */
    double a_max, a_min;   /* acceleration limits [m/s^2]                  */
    double jerk_max;       /* jerk limit [m/s^3]                           */

    /* ---- communication & sensing -------------------------------------- */
    double   delay_s;      /* V2V delay of predecessor acceleration [s]    */
    double   gap_noise_std;/* radar gap noise std-dev [m]                  */
    uint64_t seed;         /* seed for noise and heterogeneous masses      */
    int      hetero;       /* 1 = heterogeneous vehicle masses             */

    /* ---- safety threshold (fixed BEFORE running experiments) ---------- */
    double safe_gap;       /* gap below this counts as a safety event [m]  */
} SimConfig;

void        config_defaults(SimConfig *c);
const char *scenario_name(Scenario s);
int         scenario_parse(const char *str, Scenario *out);

#endif
