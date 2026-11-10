#include "config.h"
#include <string.h>

void config_defaults(SimConfig *c)
{
    c->n_vehicles    = 11;
    c->scenario      = SCEN_BRAKE;
    c->duration      = 60.0;
    c->dt            = 0.01;   /* 100 Hz control and V2V                  */
    c->substeps      = 20;     /* 0.5 ms dynamics step (stiff tyre slip)  */

    c->sync_window   = 1;      /* exact synchronisation                   */

    c->v0            = 25.0;   /* 90 km/h                                 */
    c->d0            = 5.0;
    c->headway       = 0.6;
    c->kp            = 0.45;
    c->kd            = 1.20;
    c->ka            = 1.00;
    c->veh_length    = 4.5;
    c->a_max         = 2.0;
    c->a_min         = -6.0;
    c->jerk_max      = 15.0;

    c->delay_s       = 0.0;
    c->gap_noise_std = 0.0;
    c->seed          = 42;
    c->hetero        = 0;

    c->safe_gap      = 2.0;
}

const char *scenario_name(Scenario s)
{
    return s == SCEN_CRUISE ? "cruise" : "brake";
}

int scenario_parse(const char *str, Scenario *out)
{
    if (strcmp(str, "cruise") == 0) { *out = SCEN_CRUISE; return 0; }
    if (strcmp(str, "brake")  == 0) { *out = SCEN_BRAKE;  return 0; }
    return -1;
}
