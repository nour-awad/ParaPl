/*
 * metrics.h - behavioural metrics and parallel-vs-sequential deviation.
 * Computed after the run from the full-resolution trajectory log, outside
 * the timed region.
 */
#ifndef METRICS_H
#define METRICS_H

#include "sim.h"

typedef struct {
    double min_gap;          /* smallest bumper-to-bumper gap [m]               */
    double max_abs_err;      /* largest |spacing error| over followers [m]      */
    double rms_err;          /* RMS spacing error over followers and time [m]   */
    double string_amp;       /* peak |e| of last follower / peak |e| of first    */
    double unsafe_time;      /* vehicle-seconds with gap < safe_gap             */
    int    collision;        /* 1 if any gap <= 0                               */
} Behaviour;

typedef struct {
    int    bit_identical;    /* every logged x, v, a bitwise equal              */
    double max_dx;           /* max |x_run - x_ref| [m]                          */
    double max_dv;           /* max |v_run - v_ref| [m/s]                        */
    double max_da;           /* max |a_run - a_ref| [m/s^2]                      */
    double max_dgap;         /* max |gap_run - gap_ref| [m]                      */
    double rms_dx;           /* RMS position difference [m]                      */
} Deviation;

void metrics_behaviour(const Sim *S, Behaviour *b);
void metrics_deviation(const Sim *ref, const Sim *run, Deviation *d);

/* Write t, vehicle, x, v, a, gap, spacing_error every `every` steps. */
int  write_trajectory_csv(const Sim *S, const char *path, int every);

#endif
