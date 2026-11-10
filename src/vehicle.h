/*
 * vehicle.h - longitudinal vehicle dynamics model (the per-vehicle workload).
 *
 * Reduced-order version of the vehicle model in Wu et al. (2024):
 *   - car body: longitudinal position and speed
 *   - front and rear axles: wheel spin dynamics (rotational inertia)
 *   - PAC2002 Magic Formula longitudinal tyre force (Wu et al. Table 1 coefficients)
 *   - static + dynamic axle load transfer
 *   - speed-dependent rolling resistance (Wu et al. eq. 9)
 *   - aerodynamic drag with head-wind (Wu et al. eq. 10) and drafting reduction
 *   - road grade
 *   - first-order traction/brake actuator lag
 *   - ABS as a unilateral P controller on wheel slip (as in Wu et al.)
 *
 * The model is integrated with RK4 using `substeps` sub-steps per control step,
 * so `substeps` is the model-fidelity / workload knob for the experiments.
 */
#ifndef VEHICLE_H
#define VEHICLE_H

typedef struct {
    double m;            /* total mass [kg]                       */
    double J;            /* axle rotational inertia [kg m^2]      */
    double R;            /* wheel radius [m]                      */
    double wb, lf, lr;   /* wheelbase, CG to front/rear axle [m]  */
    double hcg;          /* CG height [m]                         */
    double Cd, A;        /* drag coefficient, frontal area [m^2]  */
    double tau;          /* actuator time constant [s]            */
    double T_drive_max;  /* max traction torque at wheels [N m]   */
    double T_brake_max;  /* max brake torque at wheels [N m]      */
    double P_max;        /* max traction power [W]                */
    double Fz0f, Fz0r;   /* nominal (static) axle loads [N]       */
} VehParams;

typedef struct {
    double x;            /* position [m]                          */
    double v;            /* speed [m/s]                           */
    double wf, wr;       /* front/rear wheel speed [rad/s]        */
    double T;            /* actual wheel torque (actuator) [N m]  */
} DynState;

void   vehicle_params_init(VehParams *p, double mass);

/* Resistive force [N] (rolling + air + grade). include_wind = 0 gives the
 * controller's nominal estimate (the wind is an unknown disturbance). */
double vehicle_resistance(const VehParams *p, double x, double v, double drag_factor,
                          double t, int include_wind);

/* Traction torque limit at speed v (torque and power limited). */
double vehicle_drive_torque_limit(const VehParams *p, double v);

/* Advance the vehicle state by dt using `substeps` RK4 steps. T_cmd, the
 * drafting factor and the load-transfer acceleration are held over dt. */
void   vehicle_integrate(const VehParams *p, DynState *s, double T_cmd, double drag_factor,
                         double a_loadtransfer, double t0, double dt, int substeps);

/* Drag reduction from driving in the predecessor's wake (1 = no reduction). */
double drafting_factor(double gap);

#endif
