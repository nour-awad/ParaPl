#!/usr/bin/env bash
# Run the full experiment plan and produce tables + figures in results/.
#
#   ./scripts/run_experiments.sh            full matrix (about 10-15 min on an M1)
#   QUICK=1 ./scripts/run_experiments.sh    smoke test (short runs, 2 repeats)
#
# For trustworthy timings: close other apps, plug in the laptop, and do not
# use the machine while it runs. Every setting is recorded in results/meta.txt.
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=build/platoon
OUT=results
mkdir -p "$OUT"
make -s

# Busy-wait at barriers: per-step barriers are microseconds apart, so sleeping
# threads (the default "passive" policy) add large wake-up latency.
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-active}
export OMP_SCHEDULE=${OMP_SCHEDULE:-static}

if [[ "${QUICK:-0}" == "1" ]]; then
  COMMON=(--duration 25 -r 2 -w 1)
else
  COMMON=(--duration 60 -r 5 -w 1)
fi

# ---- metadata (reproducibility) -------------------------------------------
{
  echo "date:            $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git commit:      $(git rev-parse --short HEAD 2>/dev/null || echo n/a)"
  echo "os:              $(uname -srm)"
  if [[ "$(uname)" == "Darwin" ]]; then
    echo "cpu:             $(sysctl -n machdep.cpu.brand_string)"
    echo "cores (total):   $(sysctl -n hw.ncpu)"
    echo "perf cores:      $(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo n/a)"
    echo "eff. cores:      $(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null || echo n/a)"
    echo "memory (bytes):  $(sysctl -n hw.memsize)"
  else
    echo "cpu:             $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)"
    echo "cores (total):   $(nproc)"
  fi
  echo "compiler:        $($(make -s print-cc) --version 2>/dev/null | head -1)"
  echo "OMP_WAIT_POLICY: $OMP_WAIT_POLICY"
  echo "OMP_SCHEDULE:    $OMP_SCHEDULE"
  echo "settings:        ${COMMON[*]}"
  echo "load average:    $(uptime | sed 's/.*load average[s]*: //')"
} > "$OUT/meta.txt"
cat "$OUT/meta.txt"

# ---- 0. correctness gate (H1 precondition) ---------------------------------
"$BIN" selftest

# ---- 1. main matrix: size x threads x sync window x scenario (RQ1, H1-H3) --
"$BIN" bench -n 6,11,20,40 -p 1,2,4,8 -K 1,5,10 -s cruise,brake \
       "${COMMON[@]}" -o "$OUT/bench_main.csv"

# ---- 2. workload granularity: model fidelity vs parallel efficiency --------
: > "$OUT/bench_granularity.csv"
for SUB in 5 10 20 40; do
  "$BIN" bench -n 11 -p 1,2,4 -K 1 -s brake --substeps "$SUB" \
         "${COMMON[@]}" -o "$OUT/.gran_$SUB.csv"
  if [[ ! -s "$OUT/bench_granularity.csv" ]]; then
    cat "$OUT/.gran_$SUB.csv" > "$OUT/bench_granularity.csv"
  else
    tail -n +2 "$OUT/.gran_$SUB.csv" >> "$OUT/bench_granularity.csv"
  fi
  rm -f "$OUT/.gran_$SUB.csv"
done

# ---- 3. scheduling on heterogeneous cores (performance + efficiency) -------
: > "$OUT/bench_schedule.csv"
for SCHED in static dynamic,1; do
  OMP_SCHEDULE=$SCHED "$BIN" bench -n 20 -p 4,8 -K 1 -s brake \
         "${COMMON[@]}" -o "$OUT/.sched.csv"
  if [[ ! -s "$OUT/bench_schedule.csv" ]]; then
    cat "$OUT/.sched.csv" > "$OUT/bench_schedule.csv"
  else
    tail -n +2 "$OUT/.sched.csv" >> "$OUT/bench_schedule.csv"
  fi
  rm -f "$OUT/.sched.csv"
done

# ---- 4. trajectories for behaviour plots -----------------------------------
DUR=$([[ "${QUICK:-0}" == "1" ]] && echo 25 || echo 60)
"$BIN" run -s brake -n 11 -p 0 --duration "$DUR" --csv "$OUT/traj_brake_seq.csv"
"$BIN" run -s brake -n 11 -p 4 --duration "$DUR" --csv "$OUT/traj_brake_exact_p4.csv"
"$BIN" run -s brake -n 11 -p 4 -K 10 --duration "$DUR" --csv "$OUT/traj_brake_relaxed_K10_p4.csv"
"$BIN" run -s cruise -n 11 -p 0 --duration "$DUR" --csv "$OUT/traj_cruise_seq.csv"

# ---- 5. tables and figures --------------------------------------------------
PY=${PYTHON:-python3}
"$PY" scripts/analyze.py "$OUT"
echo "Done. See $OUT/summary.md and $OUT/figures/."
