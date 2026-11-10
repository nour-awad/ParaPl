#include "controller.h"
#include "scenario.h"

#define LEADER_KV  0.8   /* leader speed-tracking gain [1/s]            */
#define LL_KI      1.0   /* low-level integral gain [1/s]               */
#define LL_IMAX    2.0   /* integral clamp [m/s^2]                      */

double leader_accel(const SimConfig *c, double t, double v)
{
    double v_ref, a_ref;
    leader_reference(c, t, &v_ref, &a_ref);
    return a_ref + LEADER_KV * (v_ref - v);
}

double cacc_accel(const SimConfig *c, double gap, double v, double v_pred,
                  double a_self, double a_pred)
{
    const double e    = gap - (c->d0 + c->headway * v);
    const double edot = (v_pred - v) - c->headway * a_self;
    return c->ka * a_pred + c->kp * e + c->kd * edot;
}

double shape_accel(const SimConfig *c, double a_des, double a_prev)
{
    if (a_des > c->a_max) a_des = c->a_max;
    if (a_des < c->a_min) a_des = c->a_min;
    const double dmax = c->jerk_max * c->dt;
    if (a_des > a_prev + dmax) a_des = a_prev + dmax;
    if (a_des < a_prev - dmax) a_des = a_prev - dmax;
    return a_des;
}

double lowlevel_torque(const SimConfig *c, const VehParams *p, double *integ,
                       double x, double v, double a_cmd, double a_meas, double drag)
{
    *integ += LL_KI * (a_cmd - a_meas) * c->dt;
    if (*integ >  LL_IMAX) *integ =  LL_IMAX;
    if (*integ < -LL_IMAX) *integ = -LL_IMAX;

    const double F = p->m * (a_cmd + *integ) + vehicle_resistance(p, x, v, drag, 0.0, 0);
    double T = p->R * F;
    const double Tmax = vehicle_drive_torque_limit(p, v);
    if (T > Tmax)            T = Tmax;
    if (T < -p->T_brake_max) T = -p->T_brake_max;
    return T;
}
