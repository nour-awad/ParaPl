/*
 * controller.h - upper-level platoon control and lower-level torque control.
 *
 * Leader: cruise control tracking the scenario speed profile.
 * Followers: CACC with a constant time-gap spacing policy
 *     d_des = d0 + h * v_i
 *     e     = gap - d_des
 *     a_des = ka * a_{i-1}(V2V, possibly delayed) + kp * e + kd * (v_{i-1} - v_i - h * a_i)
 * Gap and relative speed come from the on-board radar; the predecessor's
 * acceleration comes over V2V.
 * Lower level: inverse-model feed-forward plus integral action converts the
 * desired acceleration into a wheel torque command.
 */
#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "config.h"
#include "vehicle.h"

double leader_accel(const SimConfig *c, double t, double v);

double cacc_accel(const SimConfig *c, double gap, double v, double v_pred,
                  double a_self, double a_pred);

/* acceleration and jerk limits */
double shape_accel(const SimConfig *c, double a_des, double a_prev);

double lowlevel_torque(const SimConfig *c, const VehParams *p, double *integ,
                       double x, double v, double a_cmd, double a_meas, double drag);

#endif
