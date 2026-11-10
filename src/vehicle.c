#include "vehicle.h"
#include "config.h"
#include "scenario.h"
#include <math.h>

/* ---- PAC2002 Magic Formula, pure longitudinal slip (Wu et al. Table 1) --- */
#define PCX1  1.685
#define PDX1  1.0
#define PDX2 -0.037
#define PEX1  0.344
#define PEX2  0.095
#define PEX3 -0.020
#define PEX4  0.0
#define PKX1 21.51
#define PKX2 -0.163
#define PKX3  0.245

/* ---- other model constants --------------------------------------------- */
#define V_EPS       1.0    /* slip regularisation speed [m/s]           */
#define FRONT_BRAKE 0.65   /* brake torque share on the front axle      */
#define ABS_SLIP    0.10   /* ABS slip threshold                        */
#define ABS_GAIN    8.0    /* ABS P gain on slip error                  */
#define DRAFT_MAX   0.30   /* max drag reduction when very close        */
#define DRAFT_LEN   15.0   /* decay length of the drafting effect [m]   */

void vehicle_params_init(VehParams *p, double mass)
{
    const double s = mass / 1500.0;          /* scale components with mass */
    p->m           = mass;
    p->J           = 1.5 * s;
    p->R           = 0.31;
    p->wb          = 2.7;
    p->lf          = 1.2;
    p->lr          = 1.5;
    p->hcg         = 0.55;
    p->Cd          = 0.30;
    p->A           = 2.2;
    p->tau         = 0.15;
    p->T_drive_max = 2500.0 * s;
    p->T_brake_max = 8000.0 * s;
    p->P_max       = 90000.0 * s;
    p->Fz0f        = mass * G_ACC * p->lr / p->wb;
    p->Fz0r        = mass * G_ACC * p->lf / p->wb;
}

static inline double tyre_fx(double kappa, double Fz, double Fz0)
{
    const double dfz = (Fz - Fz0) / Fz0;
    const double Cx  = PCX1;
    const double Dx  = (PDX1 + PDX2 * dfz) * Fz;
    double Ex = (PEX1 + PEX2 * dfz + PEX3 * dfz * dfz) * (1.0 - PEX4 * (kappa >= 0 ? 1.0 : -1.0));
    if (Ex > 1.0) Ex = 1.0;
    const double Kxk = Fz * (PKX1 + PKX2 * dfz) * exp(PKX3 * dfz);
    const double Bx  = Kxk / (Cx * Dx + 1e-9);
    const double Bk  = Bx * kappa;
    return Dx * sin(Cx * atan(Bk - Ex * (Bk - atan(Bk))));
}

static inline double rolling_coeff(double v)
{
    const double r = (v * 3.6) / 100.0;      /* (v in km/h) / 100 */
    return 0.009 + 0.0025 * r + 0.0004 * r * r * r * r;
}

static inline double abs_factor(double kappa)
{
    if (kappa >= -ABS_SLIP) return 1.0;
    double f = 1.0 - ABS_GAIN * (-kappa - ABS_SLIP);
    return f < 0.0 ? 0.0 : f;
}

double drafting_factor(double gap)
{
    if (gap < 0.0) gap = 0.0;
    return 1.0 - DRAFT_MAX * exp(-gap / DRAFT_LEN);
}

double vehicle_resistance(const VehParams *p, double x, double v, double drag_factor,
                          double t, int include_wind)
{
    const double th   = road_grade(x);
    const double W    = p->m * G_ACC;
    const double vrel = v + (include_wind ? wind_speed(t) : 0.0);
    const double Fair = 0.5 * RHO_AIR * p->Cd * drag_factor * p->A * vrel * fabs(vrel);
    const double Frol = rolling_coeff(v) * W * cos(th);
    return Fair + Frol + W * sin(th);
}

double vehicle_drive_torque_limit(const VehParams *p, double v)
{
    const double Tp = p->P_max * p->R / (v > 5.0 ? v : 5.0);
    return Tp < p->T_drive_max ? Tp : p->T_drive_max;
}

/* State derivative: y = [x, v, wf, wr, T] */
static inline void deriv(const VehParams *p, const double y[5], double T_cmd, double drag,
                         double a_lt, double t, double dy[5])
{
    const double x = y[0], v = y[1], wf = y[2], wr = y[3], T = y[4];
    const double th = road_grade(x);
    const double ct = cos(th), st = sin(th);
    const double W  = p->m * G_ACC;

    /* axle loads with grade and longitudinal load transfer */
    double Fzf = (W * ct * p->lr - (p->m * a_lt + W * st) * p->hcg) / p->wb;
    double Fzr = W * ct - Fzf;
    if (Fzf < 0.05 * W) Fzf = 0.05 * W;
    if (Fzr < 0.05 * W) Fzr = 0.05 * W;

    /* longitudinal slip */
    const double vden = fabs(v) > V_EPS ? fabs(v) : V_EPS;
    const double kf = (wf * p->R - v) / vden;
    const double kr = (wr * p->R - v) / vden;

    /* torque distribution: front-wheel drive, 65/35 braking, ABS on brakes */
    double Tf, Tr;
    if (T >= 0.0) { Tf = T; Tr = 0.0; }
    else {
        Tf = FRONT_BRAKE * T * abs_factor(kf);
        Tr = (1.0 - FRONT_BRAKE) * T * abs_factor(kr);
    }

    const double Fxf = tyre_fx(kf, Fzf, p->Fz0f);
    const double Fxr = tyre_fx(kr, Fzr, p->Fz0r);

    const double vrel = v + wind_speed(t);
    const double Fair = 0.5 * RHO_AIR * p->Cd * drag * p->A * vrel * fabs(vrel);
    double Frol = rolling_coeff(v) * W * ct;
    if (v < 0.1) Frol *= (v > 0.0 ? v / 0.1 : 0.0);

    dy[0] = v;
    dy[1] = (Fxf + Fxr - Fair - Frol - W * st) / p->m;
    dy[2] = (Tf - Fxf * p->R) / p->J;
    dy[3] = (Tr - Fxr * p->R) / p->J;
    dy[4] = (T_cmd - T) / p->tau;
}

void vehicle_integrate(const VehParams *p, DynState *s, double T_cmd, double drag_factor,
                       double a_lt, double t0, double dt, int substeps)
{
    const double h = dt / substeps;
    double y[5] = { s->x, s->v, s->wf, s->wr, s->T };
    double k1[5], k2[5], k3[5], k4[5], tmp[5];

    for (int j = 0; j < substeps; ++j) {
        const double t = t0 + j * h;
        deriv(p, y, T_cmd, drag_factor, a_lt, t, k1);
        for (int q = 0; q < 5; ++q) tmp[q] = y[q] + 0.5 * h * k1[q];
        deriv(p, tmp, T_cmd, drag_factor, a_lt, t + 0.5 * h, k2);
        for (int q = 0; q < 5; ++q) tmp[q] = y[q] + 0.5 * h * k2[q];
        deriv(p, tmp, T_cmd, drag_factor, a_lt, t + 0.5 * h, k3);
        for (int q = 0; q < 5; ++q) tmp[q] = y[q] + h * k3[q];
        deriv(p, tmp, T_cmd, drag_factor, a_lt, t + h, k4);
        for (int q = 0; q < 5; ++q)
            y[q] += (h / 6.0) * (k1[q] + 2.0 * k2[q] + 2.0 * k3[q] + k4[q]);
        if (y[2] < 0.0) y[2] = 0.0;   /* wheels do not spin backwards */
        if (y[3] < 0.0) y[3] = 0.0;
    }
    s->x = y[0]; s->v = y[1]; s->wf = y[2]; s->wr = y[3]; s->T = y[4];
}
