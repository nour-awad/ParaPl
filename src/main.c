/*
 * ParaPl - parallel vehicle platoon simulation (C + OpenMP).
 *
 *   platoon run      single run (+ sequential baseline for comparison)
 *   platoon bench    full experiment matrix -> CSV (one row per repeat)
 *   platoon selftest determinism check: exact sync must be bit-identical
 *
 * Run `platoon --help` for options.
 */
#ifdef __APPLE__
#define _DARWIN_C_SOURCE   /* pthread QoS API */
#endif
#include "config.h"
#include "metrics.h"
#include "sim.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#ifdef __APPLE__
#include <pthread.h>
#include <pthread/qos.h>
#endif

#define MAX_LIST 32

typedef struct {
    SimConfig base;
    int  vehicles[MAX_LIST], n_vehicles_list;
    int  threads[MAX_LIST],  n_threads_list;
    int  windows[MAX_LIST],  n_windows_list;
    Scenario scenarios[4];   int n_scenarios;
    int  repeats, warmup, csv_every;
    const char *out_path, *csv_path;
} Options;

static void usage(FILE *f)
{
    fprintf(f,
"Usage: platoon <run|bench|selftest> [options]\n"
"\n"
"Experiment factors (comma-separated lists in bench mode; first value used by run):\n"
"  -n, --vehicles LIST     platoon sizes                      (run: 11 | bench: 6,11,20,40)\n"
"  -p, --threads LIST      OpenMP threads, 0 = sequential     (run: 4  | bench: 1,2,4,8)\n"
"  -K, --sync-window LIST  1 = exact sync, >1 = relaxed       (run: 1  | bench: 1,5,10)\n"
"  -s, --scenario LIST     cruise | brake                     (run: brake | bench: cruise,brake)\n"
"\n"
"Held-constant model settings:\n"
"      --duration S        simulated time [s]                 (60)\n"
"      --dt S              control / V2V step [s]             (0.01)\n"
"      --substeps N        RK4 sub-steps per control step     (20)\n"
"      --headway S         CACC time gap h [s]                (0.6)\n"
"      --kp/--kd/--ka X    CACC gains                         (0.45/1.2/1.0)\n"
"      --delay-ms MS       V2V delay of predecessor accel     (0)\n"
"      --noise M           radar gap noise std-dev [m]        (0)\n"
"      --hetero            heterogeneous vehicle masses (1.2-2.4 t)\n"
"      --seed N            seed for noise / masses            (42)\n"
"      --safe-gap M        safety-event gap threshold [m]     (2.0)\n"
"\n"
"Measurement:\n"
"  -r, --repeats N         timed repeats per configuration    (run: 1 | bench: 5)\n"
"  -w, --warmup N          untimed warm-up runs               (1)\n"
"  -o, --out FILE          bench: results CSV                 (results/bench.csv)\n"
"      --csv FILE          run: write trajectory CSV of the parallel run\n"
"      --csv-every N       trajectory CSV decimation in steps (10)\n"
"\n"
"Scheduling: set OMP_SCHEDULE (default static). OMP_WAIT_POLICY=active lowers barrier latency.\n");
}

static int parse_int_list(const char *s, int *out, int max)
{
    int n = 0;
    char *buf = strdup(s), *save = NULL;
    for (char *tok = strtok_r(buf, ",", &save); tok && n < max; tok = strtok_r(NULL, ",", &save))
        out[n++] = atoi(tok);
    free(buf);
    return n;
}

static int parse_scenarios(const char *s, Scenario *out, int max)
{
    int n = 0;
    char *buf = strdup(s), *save = NULL;
    for (char *tok = strtok_r(buf, ",", &save); tok && n < max; tok = strtok_r(NULL, ",", &save)) {
        if (scenario_parse(tok, &out[n]) != 0) {
            fprintf(stderr, "unknown scenario '%s'\n", tok);
            free(buf);
            exit(2);
        }
        ++n;
    }
    free(buf);
    return n;
}

static void parse_options(int argc, char **argv, Options *o, int bench)
{
    config_defaults(&o->base);
    o->n_vehicles_list = bench ? parse_int_list("6,11,20,40", o->vehicles, MAX_LIST)
                               : parse_int_list("11", o->vehicles, MAX_LIST);
    o->n_threads_list  = bench ? parse_int_list("1,2,4,8", o->threads, MAX_LIST)
                               : parse_int_list("4", o->threads, MAX_LIST);
    o->n_windows_list  = bench ? parse_int_list("1,5,10", o->windows, MAX_LIST)
                               : parse_int_list("1", o->windows, MAX_LIST);
    o->n_scenarios     = bench ? parse_scenarios("cruise,brake", o->scenarios, 4)
                               : parse_scenarios("brake", o->scenarios, 4);
    o->repeats   = bench ? 5 : 1;
    o->warmup    = 1;
    o->csv_every = 10;
    o->out_path  = "results/bench.csv";
    o->csv_path  = NULL;

    enum { O_DURATION = 1000, O_DT, O_SUBSTEPS, O_HEADWAY, O_KP, O_KD, O_KA, O_DELAY,
           O_NOISE, O_HETERO, O_SEED, O_SAFE, O_CSV, O_CSV_EVERY };
    static const struct option longopts[] = {
        {"vehicles",    required_argument, 0, 'n'},
        {"threads",     required_argument, 0, 'p'},
        {"sync-window", required_argument, 0, 'K'},
        {"scenario",    required_argument, 0, 's'},
        {"repeats",     required_argument, 0, 'r'},
        {"warmup",      required_argument, 0, 'w'},
        {"out",         required_argument, 0, 'o'},
        {"duration",    required_argument, 0, O_DURATION},
        {"dt",          required_argument, 0, O_DT},
        {"substeps",    required_argument, 0, O_SUBSTEPS},
        {"headway",     required_argument, 0, O_HEADWAY},
        {"kp",          required_argument, 0, O_KP},
        {"kd",          required_argument, 0, O_KD},
        {"ka",          required_argument, 0, O_KA},
        {"delay-ms",    required_argument, 0, O_DELAY},
        {"noise",       required_argument, 0, O_NOISE},
        {"hetero",      no_argument,       0, O_HETERO},
        {"seed",        required_argument, 0, O_SEED},
        {"safe-gap",    required_argument, 0, O_SAFE},
        {"csv",         required_argument, 0, O_CSV},
        {"csv-every",   required_argument, 0, O_CSV_EVERY},
        {"help",        no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "n:p:K:s:r:w:o:h", longopts, NULL)) != -1) {
        switch (c) {
        case 'n': o->n_vehicles_list = parse_int_list(optarg, o->vehicles, MAX_LIST); break;
        case 'p': o->n_threads_list  = parse_int_list(optarg, o->threads,  MAX_LIST); break;
        case 'K': o->n_windows_list  = parse_int_list(optarg, o->windows,  MAX_LIST); break;
        case 's': o->n_scenarios     = parse_scenarios(optarg, o->scenarios, 4); break;
        case 'r': o->repeats   = atoi(optarg); break;
        case 'w': o->warmup    = atoi(optarg); break;
        case 'o': o->out_path  = optarg; break;
        case O_DURATION:  o->base.duration      = atof(optarg); break;
        case O_DT:        o->base.dt            = atof(optarg); break;
        case O_SUBSTEPS:  o->base.substeps      = atoi(optarg); break;
        case O_HEADWAY:   o->base.headway       = atof(optarg); break;
        case O_KP:        o->base.kp            = atof(optarg); break;
        case O_KD:        o->base.kd            = atof(optarg); break;
        case O_KA:        o->base.ka            = atof(optarg); break;
        case O_DELAY:     o->base.delay_s       = atof(optarg) / 1000.0; break;
        case O_NOISE:     o->base.gap_noise_std = atof(optarg); break;
        case O_HETERO:    o->base.hetero        = 1; break;
        case O_SEED:      o->base.seed          = strtoull(optarg, NULL, 10); break;
        case O_SAFE:      o->base.safe_gap      = atof(optarg); break;
        case O_CSV:       o->csv_path           = optarg; break;
        case O_CSV_EVERY: o->csv_every          = atoi(optarg); break;
        case 'h': usage(stdout); exit(0);
        default:  usage(stderr); exit(2);
        }
    }
    if (o->base.substeps < 1 || o->base.dt <= 0 || o->base.duration <= 0 || o->repeats < 1) {
        fprintf(stderr, "invalid settings\n");
        exit(2);
    }
}

/* Run one configuration on a fresh simulation. threads == 0 -> sequential. */
static double run_once(Sim *S, const SimConfig *cfg, int threads)
{
    if (sim_init(S, cfg) != 0) { fprintf(stderr, "out of memory\n"); exit(1); }
    return threads == 0 ? sim_run_seq(S) : sim_run_omp(S, threads);
}

static void warmup(const SimConfig *cfg, int threads, int n)
{
    for (int w = 0; w < n; ++w) { Sim t; run_once(&t, cfg, threads); sim_free(&t); }
}

/* ------------------------------------------------------------------------ */

static int cmd_run(Options *o)
{
    SimConfig cfg = o->base;
    cfg.n_vehicles  = o->vehicles[0];
    cfg.scenario    = o->scenarios[0];
    const int threads = o->threads[0];

    /* sequential exact baseline */
    SimConfig bcfg = cfg; bcfg.sync_window = 1;
    Sim base;
    warmup(&bcfg, 0, o->warmup);
    const double t_seq = run_once(&base, &bcfg, 0);

    cfg.sync_window = o->windows[0];
    warmup(&cfg, threads, o->warmup);
    double best = INFINITY, sum = 0.0;
    Sim S;
    for (int r = 0; r < o->repeats; ++r) {
        if (r) sim_free(&S);
        const double t = run_once(&S, &cfg, threads);
        sum += t; if (t < best) best = t;
    }
    const double t_par = sum / o->repeats;

    Behaviour bb, bp; Deviation d;
    metrics_behaviour(&base, &bb);
    metrics_behaviour(&S, &bp);
    metrics_deviation(&base, &S, &d);

    printf("scenario=%s vehicles=%d steps=%d substeps=%d procs=%d\n",
           scenario_name(cfg.scenario), cfg.n_vehicles, S.nsteps, cfg.substeps, omp_get_num_procs());
    printf("sequential     : %.4f s  (RTF %.4f)\n", t_seq, t_seq / cfg.duration);
    printf("%-15s: %.4f s  (RTF %.4f)  [%s, K=%d]\n",
           threads ? "parallel" : "sequential(2)", t_par, t_par / cfg.duration,
           cfg.sync_window == 1 ? "exact" : "relaxed", cfg.sync_window);
    if (threads)
        printf("speedup        : %.3f   efficiency %.1f %%  (threads=%d)\n",
               t_seq / t_par, 100.0 * t_seq / t_par / threads, threads);
    printf("behaviour      : min_gap %.3f m  max|e| %.3f m  rms(e) %.4f m  string_amp %.3f  unsafe %.2f veh*s  collision %d\n",
           bp.min_gap, bp.max_abs_err, bp.rms_err, bp.string_amp, bp.unsafe_time, bp.collision);
    printf("vs sequential  : bit_identical=%s  max|dx| %.3e m  max|dv| %.3e m/s  max|da| %.3e  max|dgap| %.3e m\n",
           d.bit_identical ? "yes" : "no", d.max_dx, d.max_dv, d.max_da, d.max_dgap);

    if (o->csv_path) {
        if (write_trajectory_csv(&S, o->csv_path, o->csv_every) != 0)
            fprintf(stderr, "could not write %s\n", o->csv_path);
        else
            printf("trajectory     : %s\n", o->csv_path);
    }
    sim_free(&S);
    sim_free(&base);
    return 0;
}

static const char *schedule_name(void)
{
    static char buf[32];
    omp_sched_t kind; int chunk;
    omp_get_schedule(&kind, &chunk);
    const char *k = "other";
    switch ((int)kind & ~(int)omp_sched_monotonic) {
    case omp_sched_static:  k = "static";  break;
    case omp_sched_dynamic: k = "dynamic"; break;
    case omp_sched_guided:  k = "guided";  break;
    case omp_sched_auto:    k = "auto";    break;
    }
    if (chunk > 0) snprintf(buf, sizeof buf, "%s:%d", k, chunk);
    else           snprintf(buf, sizeof buf, "%s", k);
    return buf;
}

static void bench_row(FILE *f, const SimConfig *cfg, const char *mode, int threads, int rep,
                      double wall, const Behaviour *b, const Deviation *d)
{
    fprintf(f, "%s,%d,%s,%s,%d,%d,%s,%d,%d,%.6f,%.3f,%.6f,%d,%.6e,%.6e,%.6e,%.6e,%.6f,%.6f,%.6f,%.6f,%.4f,%d\n",
            scenario_name(cfg->scenario), cfg->n_vehicles, mode,
            cfg->sync_window == 1 ? "exact" : "relaxed", cfg->sync_window, threads,
            threads ? schedule_name() : "none", cfg->substeps, rep,
            wall, cfg->duration, wall / cfg->duration,
            d->bit_identical, d->max_dx, d->max_dv, d->max_da, d->max_dgap,
            b->min_gap, b->max_abs_err, b->rms_err, b->string_amp, b->unsafe_time, b->collision);
    fflush(f);
}

static int cmd_bench(Options *o)
{
    FILE *f = fopen(o->out_path, "w");
    if (!f) { perror(o->out_path); return 1; }
    fprintf(f, "scenario,vehicles,mode,sync,K,threads,schedule,substeps,repeat,wall_s,sim_s,rtf,bit_identical,"
               "max_dx,max_dv,max_da,max_dgap,min_gap,max_abs_err,rms_err,string_amp,unsafe_time,collision\n");

    for (int si = 0; si < o->n_scenarios; ++si)
    for (int vi = 0; vi < o->n_vehicles_list; ++vi) {
        SimConfig cfg = o->base;
        cfg.scenario    = o->scenarios[si];
        cfg.n_vehicles  = o->vehicles[vi];
        cfg.sync_window = 1;

        /* sequential baseline: the reference for timing and behaviour */
        Sim base;
        warmup(&cfg, 0, o->warmup);
        for (int r = 0; r < o->repeats; ++r) {
            Sim S; Behaviour b; Deviation d;
            const double t = run_once(&S, &cfg, 0);
            metrics_behaviour(&S, &b);
            if (r == 0) { base = S; metrics_deviation(&base, &S, &d); }
            else        { metrics_deviation(&base, &S, &d); sim_free(&S); }
            bench_row(f, &cfg, "seq", 0, r, t, &b, &d);
        }
        fprintf(stderr, "[%s N=%d] sequential done\n", scenario_name(cfg.scenario), cfg.n_vehicles);

        for (int ki = 0; ki < o->n_windows_list; ++ki)
        for (int ti = 0; ti < o->n_threads_list; ++ti) {
            cfg.sync_window = o->windows[ki];
            const int p = o->threads[ti];
            warmup(&cfg, p, o->warmup);
            for (int r = 0; r < o->repeats; ++r) {
                Sim S; Behaviour b; Deviation d;
                const double t = run_once(&S, &cfg, p);
                metrics_behaviour(&S, &b);
                metrics_deviation(&base, &S, &d);
                bench_row(f, &cfg, p == 0 ? "seq" : "omp", p, r, t, &b, &d);
                sim_free(&S);
            }
            fprintf(stderr, "[%s N=%d] K=%d threads=%d done\n",
                    scenario_name(cfg.scenario), cfg.n_vehicles, cfg.sync_window, p);
        }
        sim_free(&base);
    }
    fclose(f);
    fprintf(stderr, "results written to %s\n", o->out_path);
    return 0;
}

/* H1 check: exact synchronisation must reproduce the sequential run bit for
 * bit, for any thread count, also with noise, V2V delay and mixed masses. */
static int cmd_selftest(Options *o)
{
    SimConfig cfg = o->base;
    cfg.n_vehicles = 11; cfg.scenario = SCEN_BRAKE; cfg.duration = 30.0;
    cfg.gap_noise_std = 0.05; cfg.delay_s = 0.05; cfg.hetero = 1;

    int fails = 0;
    Sim base;
    run_once(&base, &cfg, 0);

    const int ths[] = {1, 2, 3, 4, 8, 16};
    for (size_t q = 0; q < sizeof ths / sizeof *ths; ++q) {
        Sim S; Deviation d;
        run_once(&S, &cfg, ths[q]);
        metrics_deviation(&base, &S, &d);
        printf("exact   K=1  threads=%-2d bit_identical=%s\n", ths[q], d.bit_identical ? "yes" : "NO");
        fails += !d.bit_identical;
        sim_free(&S);
    }

    /* relaxed sync must be deterministic too (same result for any thread
     * count), but is expected to deviate from the exact baseline */
    SimConfig rc = cfg; rc.sync_window = 5;
    Sim r1, r4; Deviation d14, d0;
    run_once(&r1, &rc, 1);
    run_once(&r4, &rc, 4);
    metrics_deviation(&r1, &r4, &d14);
    metrics_deviation(&base, &r4, &d0);
    printf("relaxed K=5  threads 1 vs 4 bit_identical=%s\n", d14.bit_identical ? "yes" : "NO");
    printf("relaxed K=5  vs exact baseline: max|dx| %.3e m (expected > 0)\n", d0.max_dx);
    fails += !d14.bit_identical;
    fails += d0.bit_identical;   /* would mean relaxed mode is not doing anything */

    Behaviour b; metrics_behaviour(&base, &b);
    printf("baseline behaviour: min_gap %.3f m, collision %d\n", b.min_gap, b.collision);
    fails += b.collision;

    sim_free(&r1); sim_free(&r4); sim_free(&base);
    printf(fails ? "SELFTEST FAILED\n" : "SELFTEST PASSED\n");
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        usage(argc < 2 ? stderr : stdout);
        return argc < 2 ? 2 : 0;
    }
    const char *cmd = argv[1];
    const int bench = strcmp(cmd, "bench") == 0;

#ifdef __APPLE__
    /* Ask macOS to prefer performance cores for this thread; OpenMP worker
     * threads created from it inherit the QoS class. */
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif

    /* default schedule: static, unless the user sets OMP_SCHEDULE */
    if (!getenv("OMP_SCHEDULE")) omp_set_schedule(omp_sched_static, 0);

    Options o;
    parse_options(argc - 1, argv + 1, &o, bench);

    if (strcmp(cmd, "run") == 0)      return cmd_run(&o);
    if (bench)                        return cmd_bench(&o);
    if (strcmp(cmd, "selftest") == 0) return cmd_selftest(&o);
    usage(stderr);
    return 2;
}
