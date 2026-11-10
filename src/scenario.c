#include "scenario.h"
#include <math.h>

#define TWO_PI 6.283185307179586

/* Scenario B (leader braking) timeline */
#define BRAKE_T0     15.0   /* braking starts [s]               */
#define BRAKE_DECEL  -5.0   /* hard braking [m/s^2]             */
#define BRAKE_VLOW   10.0   /* speed after braking [m/s]        */
#define BRAKE_HOLD   30.0   /* low-speed phase ends [s]         */
#define BRAKE_ACCEL   1.2   /* recovery acceleration [m/s^2]    */

void leader_reference(const SimConfig *c, double t, double *v_ref, double *a_ref)
{
    const double v0 = c->v0;

    if (c->scenario == SCEN_CRUISE) {          /* Scenario A: steady cruise */
        *v_ref = v0;
        *a_ref = 0.0;
        return;
    }

    /* Scenario B: cruise -> hard brake -> hold -> recover -> cruise */
    const double t1 = BRAKE_T0 + (BRAKE_VLOW - v0) / BRAKE_DECEL;  /* end of braking  */
    const double t3 = BRAKE_HOLD + (v0 - BRAKE_VLOW) / BRAKE_ACCEL; /* end of recovery */

    if (t < BRAKE_T0)        { *v_ref = v0; *a_ref = 0.0; }
    else if (t < t1)         { *v_ref = v0 + BRAKE_DECEL * (t - BRAKE_T0); *a_ref = BRAKE_DECEL; }
    else if (t < BRAKE_HOLD) { *v_ref = BRAKE_VLOW; *a_ref = 0.0; }
    else if (t < t3)         { *v_ref = BRAKE_VLOW + BRAKE_ACCEL * (t - BRAKE_HOLD); *a_ref = BRAKE_ACCEL; }
    else                     { *v_ref = v0; *a_ref = 0.0; }
}

double road_grade(double x)
{
    /* rolling terrain, max grade about 2.3 % */
    return 0.015 * sin(TWO_PI * x / 900.0) + 0.008 * sin(TWO_PI * x / 260.0 + 0.9);
}

double wind_speed(double t)
{
    /* gusty head-wind, deterministic so all runs see the same disturbance */
    return 3.0 + 1.5 * sin(TWO_PI * t / 40.0) + 0.7 * sin(TWO_PI * t / 7.3 + 1.3);
}
