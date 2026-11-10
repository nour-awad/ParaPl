# ParaPl: Parallel Computing in Vehicle Platoon Simulation

CNM603 Net-Centric Computing, Group 4 (Parallel Systems × Platooning).

**Research question:** Does running a vehicle platoon simulation in parallel change how the simulated platoon behaves, and how much faster is it?

This repository holds the simulation for that study. A longitudinal platoon (leader plus CACC followers) is simulated twice with the same model, scenario and initial conditions:

* **Sequential baseline:** vehicles are updated one after another on one core.
* **Parallel (OpenMP):** vehicles are distributed over threads, with **exact** per-step synchronisation (barrier every control step) or **relaxed** synchronisation (barrier every *K* steps).

The program measures execution time, speedup, efficiency and real-time factor. It also measures how far the parallel run's trajectories, gaps and speeds deviate from the sequential run.

---

## Quick start

```bash
# macOS (Apple clang needs libomp)    | Linux
brew install libomp                    # sudo apt install build-essential
make            # builds build/platoon
make test       # determinism self-test (must print SELFTEST PASSED)
```

Single run (sequential baseline + 4-thread run, leader-braking scenario, 11 vehicles):

```bash
./build/platoon run -n 11 -p 4 -s brake --csv results/traj.csv
```

Full experiment plan (tables + figures in `results/`):

```bash
python3 -m pip install -r requirements.txt   # numpy + matplotlib, for plots only
./scripts/run_experiments.sh                 # ~15-20 min on an M1; QUICK=1 for a smoke test
```

For timing you will report, close other apps, keep the laptop on power, and leave it alone while the script runs. `results/meta.txt` records the hardware, compiler, OpenMP settings and system load.

---

## Why C + OpenMP (and not Python, MPI or CUDA)

The tool choice in Presentation 2 (Python + NumPy + multiprocessing) was marked *provisional*. We changed it for these reasons:

| Option | Fit for this study |
|---|---|
| **Python multiprocessing** | Every control step needs a synchronisation (exact sync = 6,000 barriers per 60 s run). Each one would cost inter-process pickling and messaging of roughly 10-100 µs, while one vehicle's update costs about 10 µs. The overhead would dominate, so the measurement would describe Python's IPC rather than parallel platoon simulation. The GIL also rules out plain threads. |
| **NumPy vectorisation** | This is data parallelism on one core (SIMD), not multi-core execution, so it cannot answer RQ1 about cores. |
| **C + OpenMP (chosen)** | Shared memory matches the design's "shared simulation state" directly. A barrier costs about 1 µs, and threads are created once. Running the same function sequentially and in parallel makes a bit-exact behaviour comparison possible (H1). It runs on every lab laptop. |
| **MPI** (used by Wu et al. 2024) | Distributed memory, needed across machines. On a single PC it adds message-passing cost with no benefit. It is the natural future-work step to a cluster. |
| **CUDA** | Requires an NVIDIA GPU, which Apple-silicon Macs do not have. A platoon of 6-40 vehicles is also too few parallel tasks for a GPU (thousands of threads), and a kernel launch plus sync per step costs about 5-10 µs. CUDA becomes relevant for thousands of vehicles or Monte-Carlo batches of whole platoons. |

Python is still used, but only for analysis and plots (`scripts/analyze.py`), outside every timed region.

---

## Mapping the research design to the code

| Research design (Presentation 2) | Implementation |
|---|---|
| Scenario A: steady cruise | `-s cruise`: leader at 25 m/s with road-grade and gusty head-wind disturbances |
| Scenario B: leader braking | `-s brake`: 25 m/s → brake at −5 m/s² at t = 15 s → 10 m/s → recover at 1.2 m/s² from t = 30 s |
| Shared simulation state at time *t* | `ring[step][vehicle]`: published (x, v, a) per vehicle (`src/sim.h`) |
| Coordinator / workers | OpenMP team; `omp for` distributes vehicles; its implicit barrier collects results before the next step (`sim_run_omp`) |
| Synchronisation rule: consistent snapshot of the previous step | Vehicles read only data published at the last sync point (Jacobi update); see `vehicle_advance` |
| Exact vs relaxed sync | `-K 1` (exact) vs `-K 5`, `-K 10` (barrier every K steps; predecessor extrapolated at constant acceleration) |
| Held constant | `SimConfig` in `src/config.h`: model, parameters, initial state, dt, duration, leader input, safety threshold |
| Varied | execution mode, threads (`-p`), platoon size (`-n`), scenario (`-s`), repeats (`-r`) |
| 1-worker parallel run to expose overhead | `-p 1` (OpenMP with one thread) vs `-p 0` (plain sequential loop) |
| Algorithm steps 1-6 (initialise → … → record) | `sim_init` → `vehicle_advance` (read snapshot, compute) → barrier → advance → log → `metrics.c` |

### Vehicle and control model (per-vehicle workload)

A reduced-order version of the vehicle model in Wu et al. (2024):

* body: longitudinal position and speed; front and rear axles with wheel-spin dynamics
* **PAC2002 Magic Formula** longitudinal tyre force with the coefficients of Wu et al. Table 1
* axle-load transfer, road grade, speed-dependent rolling resistance (Wu et al. eq. 9), aerodynamic drag with head-wind (eq. 10), and a drafting drag reduction for followers
* first-order traction/brake actuator (τ = 0.15 s), torque/power limits, ABS as a unilateral P controller on slip
* RK4 integration with `--substeps` sub-steps per 10 ms control step (default 20, i.e. 0.5 ms), which is the model-fidelity / workload knob

Control:

* **Leader:** cruise control tracking the scenario speed profile.
* **Followers:** CACC with a constant time-gap policy, `d_des = d0 + h·v` (d0 = 5 m, h = 0.6 s), `a_des = ka·a_pred + kp·e + kd·ė`. Gap and relative speed come from radar; the predecessor's acceleration comes over V2V (optional delay `--delay-ms`).
* **Lower level:** inverse-model feed-forward plus integral action, acceleration and jerk limits.

Optional realism switches: `--noise` (radar noise), `--delay-ms` (V2V delay), `--hetero` (1.2-2.4 t masses).

### Reproducibility guarantees

* Radar noise and heterogeneous masses come from a counter-based RNG keyed by (seed, vehicle, step). The random numbers are therefore identical regardless of which thread runs a vehicle.
* No `-ffast-math`, and `-ffp-contract=off`. The sequential and parallel paths call the same `vehicle_advance()`.
* Per-vehicle data is padded to 128-byte cache lines (Apple M-series), so there is no false sharing.
* `platoon selftest` checks that exact sync is **bit-identical** to the sequential run for 1, 2, 3, 4, 8 and 16 threads, with noise, V2V delay and mixed masses switched on.

---

## Experiments (`scripts/run_experiments.sh`)

| # | Experiment | Factors | Output |
|---|---|---|---|
| 0 | Correctness gate | selftest | console |
| 1 | Main matrix (RQ1, H1-H3) | N = 6, 11, 20, 40 × threads = 1, 2, 4, 8 × K = 1, 5, 10 × cruise/brake, 5 repeats + 1 warm-up | `bench_main.csv` |
| 2 | Workload granularity | substeps = 5, 10, 20, 40 (N = 11, 1-4 threads) | `bench_granularity.csv` |
| 3 | Scheduling on heterogeneous cores | static vs dynamic, 4/8 threads (N = 20) | `bench_schedule.csv` |
| 4 | Trajectories for plots | brake: seq / exact p = 4 / relaxed K = 10; cruise: seq | `traj_*.csv` |

`scripts/analyze.py` writes `results/summary.md` (all tables plus H1-H3 verdicts) and `results/figures/*.png`. The raw CSVs keep one row per repeat.

### Metrics

* **Computational:** wall time (mean ± sd), speedup `T_seq / T_p`, efficiency `speedup / p`, real-time factor `T_wall / T_sim`.
* **Behavioural:** spacing error `e = gap − (d0 + h·v)`, min gap, RMS/max |e|, string amplification (peak |e| of the last follower divided by that of the first), unsafe time (gap < 2 m), collisions.
* **Consistency vs sequential:** bit-identical flag, max |Δx|, |Δv|, |Δa|, |Δgap|.

Tolerances for "behaviour matches" are set at the top of `scripts/analyze.py` (`TOL`, currently 0.1 m / 0.1 m/s). The research design requires the group to **agree on these before the final runs** and not change them afterwards.

### Hypotheses

* **H1:** with exact sync, parallel and sequential behaviour match within tolerance. The selftest shows they are bit-identical.
* **H2:** speedup grows with cores, with efficiency below 100 %.
* **H3:** relaxed sync increases speedup at the cost of trajectory deviation.

---

## Notes and limitations (for the paper)

* **Apple M1 = 4 performance + 4 efficiency cores.** Up to 4 threads the platform is homogeneous. Beyond that, threads land on slower efficiency cores, and with a static schedule every barrier waits for the slowest thread. Expect efficiency to drop sharply at 8 threads. Experiment 3 shows how much dynamic scheduling recovers. Report this; do not hide it.
* `OMP_WAIT_POLICY=active` (set by the script) keeps threads spinning at barriers. Sleeping threads add wake-up latency comparable to a whole control step. When other programs compete for the cores, spinning threads can be preempted, so measure on an idle machine.
* **Load balance:** with N vehicles on p threads, the ideal speedup is `N / ceil(N/p)` (for example 6 vehicles on 4 threads is at most 3×).
* The vehicle model is far cheaper than Wu et al.'s multibody model, so absolute real-time factors are not comparable to theirs. Speedup and efficiency trends are comparable. `--substeps` scales the per-vehicle cost.
* This is simulated behaviour only. Consistency between execution modes says nothing about real-vehicle safety.

## Layout

```
src/        config.[ch]  scenario.[ch]  vehicle.[ch]  controller.[ch]  sim.[ch]  metrics.[ch]  main.c
scripts/    run_experiments.sh  analyze.py
results/    generated (git-ignored; commit the final dataset with git add -f)
```

## References

* Wu, Q., Ge, X., Han, Q.-L., Cole, C., & Spiryagin, M. (2024). A scalable parallel computing method for autonomous platoons. *Vehicle System Dynamics*, 62(9), 2283-2303.
* Wu, Q., Ge, X., Cole, C., & Spiryagin, M. (2024). Car and train platoon simulations. *Journal of Vibration and Control*, 30(1-2), 3-13.
* Dai, Y., Hu, L., Liu, Y., & Yang, A. (2025). Modeling and analyzing V2V communication limitations impacts on connected and automated vehicle platoons. *PLOS ONE*, 20(8).
* Munawar, A. et al. (2013). On-time data exchange in fully-parallelized co-simulation with conservative synchronization. *Winter Simulation Conference*.
* Ploeg, J. et al. (2011). Design and experimental evaluation of cooperative adaptive cruise control. *IEEE ITSC*.
