/*
 * scenario.h - leader input and environment (road grade, wind).
 *
 * All functions are pure and deterministic (functions of time or position
 * only), so every execution mode sees exactly the same disturbances.
 */
#ifndef SCENARIO_H
#define SCENARIO_H

#include "config.h"

/* Leader reference speed and acceleration at time t. */
void   leader_reference(const SimConfig *c, double t, double *v_ref, double *a_ref);

/* Road grade angle [rad] at longitudinal position x [m]. */
double road_grade(double x);

/* Head-wind speed [m/s] at time t (positive = against the platoon). */
double wind_speed(double t);

#endif
