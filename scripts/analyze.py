#!/usr/bin/env python3
"""Turn the raw benchmark CSVs into tables (summary.md) and figures (figures/*.png).

Usage: python3 scripts/analyze.py [results_dir]

Python is only used for analysis and plotting. All timed simulation runs are
done by the C/OpenMP binary, so Python overhead never enters the measurements.
"""
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# --------------------------------------------------------------------------
# Behavioural tolerances. Per the research design these must be agreed by the
# group BEFORE the final experiments are run, and not changed afterwards.
# --------------------------------------------------------------------------
TOL = {
    "max_dx": 0.10,    # m    - max trajectory (position) difference
    "max_dgap": 0.10,  # m    - max inter-vehicle gap difference
    "max_dv": 0.10,    # m/s  - max speed difference
}

# Categorical slots in fixed order (reference palette, light mode)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE_RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95",
             "#104281", "#0d366b"]

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2.0, "lines.markersize": 6, "font.size": 10,
    "legend.frameon": False, "figure.dpi": 150,
})

NUM = {"vehicles", "K", "threads", "substeps", "repeat", "wall_s", "sim_s", "rtf", "bit_identical",
       "max_dx", "max_dv", "max_da", "max_dgap", "min_gap", "max_abs_err", "rms_err", "string_amp",
       "unsafe_time", "collision"}


def load(path):
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            for k in NUM:
                if k in r:
                    r[k] = float(r[k]) if r[k] not in ("", "nan") else math.nan
            rows.append(r)
    return rows


def group(rows, keys):
    g = defaultdict(list)
    for r in rows:
        g[tuple(r[k] for k in keys)].append(r)
    return g


def seq_times(rows):
    """Mean/std sequential time per (scenario, vehicles, substeps)."""
    out = {}
    for key, rs in group([r for r in rows if r["mode"] == "seq" and r["K"] == 1],
                         ["scenario", "vehicles", "substeps"]).items():
        t = np.array([r["wall_s"] for r in rs])
        out[key] = (t.mean(), t.std(ddof=1) if len(t) > 1 else 0.0)
    return out


def perf_table(rows):
    """One record per configuration with mean/std time, speedup and efficiency."""
    base = seq_times(rows)
    recs = []
    for key, rs in sorted(group(rows, ["scenario", "vehicles", "substeps", "mode", "sync", "K",
                                       "threads", "schedule"]).items(),
                          key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], kv[0][5], kv[0][6])):
        scen, n, sub, mode, sync, K, p, sched = key
        t = np.array([r["wall_s"] for r in rs])
        tseq = base[(scen, n, sub)][0]
        sp = tseq / t                      # per-repeat speedup against mean baseline
        recs.append(dict(
            scenario=scen, vehicles=int(n), substeps=int(sub), mode=mode, sync=sync, K=int(K),
            threads=int(p), schedule=sched, repeats=len(t),
            t_mean=t.mean(), t_std=t.std(ddof=1) if len(t) > 1 else 0.0,
            rtf=t.mean() / rs[0]["sim_s"],
            speedup=tseq / t.mean(), speedup_std=sp.std(ddof=1) if len(sp) > 1 else 0.0,
            efficiency=(tseq / t.mean()) / p if p > 0 else 1.0,
            bit_identical=all(r["bit_identical"] == 1 for r in rs),
            max_dx=max(r["max_dx"] for r in rs), max_dv=max(r["max_dv"] for r in rs),
            max_da=max(r["max_da"] for r in rs), max_dgap=max(r["max_dgap"] for r in rs),
            min_gap=rs[0]["min_gap"], max_abs_err=rs[0]["max_abs_err"], rms_err=rs[0]["rms_err"],
            string_amp=rs[0]["string_amp"], unsafe_time=rs[0]["unsafe_time"],
            collision=int(rs[0]["collision"]),
        ))
    return recs


def within_tol(r):
    return all(r[k] <= v for k, v in TOL.items())


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(out)


def write_csv(path, recs):
    if not recs:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)


# ------------------------------------------------------------------ figures --

def label_ends(ax, items, min_frac=0.06):
    """Direct labels at line ends, nudged apart vertically so they never overlap.
    items: list of (x, y, text)."""
    if not items:
        return
    lo, hi = ax.get_ylim()
    gap = (hi - lo) * min_frac
    items = sorted(items, key=lambda it: it[1])
    ys = [it[1] for it in items]
    for k in range(1, len(ys)):
        ys[k] = max(ys[k], ys[k - 1] + gap)
    for (x, y, text), yl in zip(items, ys):
        ax.annotate(text, (x, y), xytext=(x, yl), textcoords="data", va="center", fontsize=9,
                    color=INK2, xycoords="data", annotation_clip=False,
                    bbox=None, ha="left").set_position((x + 0.15, yl))


def fig_speedup_eff(recs, figdir):
    exact = [r for r in recs if r["mode"] == "omp" and r["sync"] == "exact"]
    scens = sorted({r["scenario"] for r in exact})
    sizes = sorted({r["vehicles"] for r in exact})
    if not exact:
        return
    for metric, ylabel, fname in [("speedup", "Speedup  (T_seq / T_p)", "speedup.png"),
                                  ("efficiency", "Parallel efficiency  (speedup / p)", "efficiency.png")]:
        fig, axes = plt.subplots(1, len(scens), figsize=(5.2 * len(scens), 4), sharey=True, squeeze=False)
        for ax, scen in zip(axes[0], scens):
            pmax = max(r["threads"] for r in exact)
            if metric == "speedup":
                ax.plot([1, pmax], [1, pmax], ls="--", lw=1, color=MUTED, label="ideal")
            else:
                ax.axhline(1.0, ls="--", lw=1, color=MUTED, label="ideal")
            ends = []
            for ci, n in enumerate(sizes):
                rs = sorted([r for r in exact if r["scenario"] == scen and r["vehicles"] == n],
                            key=lambda r: r["threads"])
                x = [r["threads"] for r in rs]
                y = [r[metric] for r in rs]
                err = [r["speedup_std"] / (r["threads"] if metric == "efficiency" else 1) for r in rs]
                ax.errorbar(x, y, yerr=err, marker="o", color=SERIES[ci], label=f"{n} vehicles",
                            capsize=3, elinewidth=1)
                ends.append((x[-1], y[-1], f"N={n}"))
            ax.set_title(f"Scenario: {scen}")
            label_ends(ax, ends)
            ax.set_xlabel("Threads (cores)")
            ax.set_xticks(sorted({r["threads"] for r in exact}))
        axes[0][0].set_ylabel(ylabel)
        axes[0][-1].legend(loc="best")
        fig.suptitle(f"Exact synchronisation - {ylabel.split('(')[0].strip().lower()}", color=INK)
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, fname))
        plt.close(fig)


def fig_time_vs_n(recs, figdir):
    rs = [r for r in recs if r["scenario"] == "brake" and r["sync"] == "exact"]
    if not rs:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    configs = sorted({r["threads"] for r in rs})
    for ci, p in enumerate(configs):
        sel = sorted([r for r in rs if r["threads"] == p], key=lambda r: r["vehicles"])
        lab = "sequential" if p == 0 else f"{p} thread{'s' if p > 1 else ''}"
        x = [r["vehicles"] for r in sel]
        y = [r["t_mean"] for r in sel]
        ax.errorbar(x, y, yerr=[r["t_std"] for r in sel], marker="o", color=SERIES[ci],
                    label=lab, capsize=3, elinewidth=1)
    ax.set_xlabel("Platoon size (vehicles)")
    ax.set_ylabel("Wall-clock time per run [s]")
    ax.set_title("Execution time vs platoon size (leader braking, exact sync)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "time_vs_vehicles.png"))
    plt.close(fig)


def fig_relaxed(recs, figdir, threads=4):
    rs = [r for r in recs if r["mode"] == "omp" and r["threads"] == threads and r["scenario"] == "brake"]
    if not rs or len({r["K"] for r in rs}) < 2:
        return
    sizes = sorted({r["vehicles"] for r in rs})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4))
    for ci, n in enumerate(sizes):
        sel = sorted([r for r in rs if r["vehicles"] == n], key=lambda r: r["K"])
        K = [r["K"] for r in sel]
        a1.plot(K, [r["speedup"] for r in sel], marker="o", color=SERIES[ci], label=f"{n} vehicles")
        a2.plot(K, [r["max_dx"] for r in sel], marker="o", color=SERIES[ci], label=f"{n} vehicles")
    a2.axhline(TOL["max_dx"], ls="--", lw=1, color=MUTED)
    a2.annotate(f"tolerance {TOL['max_dx']} m", (min(r["K"] for r in rs), TOL["max_dx"]), xytext=(4, 4),
                textcoords="offset points", fontsize=9, color=INK2)
    a2.set_ylim(bottom=0)
    for a in (a1, a2):
        a.set_xlabel("Sync window K (control steps between barriers; 1 = exact)")
        a.set_xticks(sorted({r["K"] for r in rs}))
    a1.set_ylabel(f"Speedup on {threads} threads")
    a1.set_title("Performance gain")
    a2.set_ylabel("Max trajectory difference vs sequential [m]")
    a2.set_title("Behaviour cost (K = 1 is bit-identical: 0 m)")
    a1.legend()
    fig.suptitle("Relaxed synchronisation trade-off (leader braking)", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "relaxed_tradeoff.png"))
    plt.close(fig)


def load_traj(path):
    if not os.path.exists(path):
        return None
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d


def fig_spacing(results, figdir, name, title):
    d = load_traj(os.path.join(results, name))
    if d is None:
        return
    veh = d["vehicle"].astype(int)
    n = veh.max() + 1
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for i in range(n):
        m = veh == i
        col = BLUE_RAMP[min(i, len(BLUE_RAMP) - 1)] if i > 0 else INK
        a1.plot(d["t"][m], d["v"][m], lw=1.5 if i else 2.2, color=col)
        if i > 0:
            a2.plot(d["t"][m], d["spacing_error"][m], lw=1.5, color=col)
    a1.set_ylabel("Speed [m/s]")
    a1.set_title(title)
    a2.set_ylabel("Spacing error [m]")
    a2.set_xlabel("Time [s]")
    a1.annotate("leader (black); followers light -> dark = front -> back of platoon",
                (0.01, 0.04), xycoords="axes fraction", fontsize=8, color=INK2)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, name.replace(".csv", ".png")))
    plt.close(fig)


def fig_traj_diff(results, figdir):
    ref = load_traj(os.path.join(results, "traj_brake_seq.csv"))
    rel = load_traj(os.path.join(results, "traj_brake_relaxed_K10_p4.csv"))
    exa = load_traj(os.path.join(results, "traj_brake_exact_p4.csv"))
    if ref is None or rel is None:
        return
    veh = ref["vehicle"].astype(int)
    n = veh.max() + 1
    fig, ax = plt.subplots(figsize=(8, 3.8))
    for i in range(1, n):
        m = veh == i
        ax.plot(ref["t"][m], rel["gap"][m] - ref["gap"][m], lw=1.4, color=BLUE_RAMP[min(i, 9)])
    if exa is not None:
        dmax = np.nanmax(np.abs(exa["x"] - ref["x"]))
        ax.annotate(f"exact sync (4 threads): max |dx| = {dmax:g} m (bit-identical)",
                    (0.01, 0.92), xycoords="axes fraction", fontsize=9, color=INK2)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Gap difference vs sequential [m]")
    ax.set_title("Relaxed sync (K = 10, 4 threads): gap deviation per follower")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "relaxed_gap_difference.png"))
    plt.close(fig)


def fig_granularity(recs, figdir):
    rs = [r for r in recs if r["mode"] == "omp"]
    if len({r["substeps"] for r in rs}) < 2:
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    for ci, p in enumerate(sorted({r["threads"] for r in rs})):
        sel = sorted([r for r in rs if r["threads"] == p], key=lambda r: r["substeps"])
        ax.plot([r["substeps"] for r in sel], [r["efficiency"] for r in sel], marker="o",
                color=SERIES[ci], label=f"{p} thread{'s' if p > 1 else ''}")
    ax.set_xscale("log", base=2)
    ax.set_xticks(sorted({r["substeps"] for r in rs}))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("RK4 sub-steps per control step (work per vehicle per step)")
    ax.set_ylabel("Parallel efficiency")
    ax.set_title("Workload granularity vs efficiency (11 vehicles, braking)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "granularity.png"))
    plt.close(fig)


# ------------------------------------------------------------------- report --

def fmt(x, nd=3):
    if isinstance(x, float) and (math.isnan(x)):
        return "-"
    return f"{x:.{nd}f}" if isinstance(x, float) else str(x)


def main():
    results = sys.argv[1] if len(sys.argv) > 1 else "results"
    figdir = os.path.join(results, "figures")
    os.makedirs(figdir, exist_ok=True)

    main_rows = load(os.path.join(results, "bench_main.csv"))
    gran_rows = load(os.path.join(results, "bench_granularity.csv"))
    sched_rows = load(os.path.join(results, "bench_schedule.csv"))
    if not main_rows:
        sys.exit(f"no {results}/bench_main.csv - run scripts/run_experiments.sh first")

    recs = perf_table(main_rows)
    write_csv(os.path.join(results, "summary_main.csv"), recs)
    gran = perf_table(gran_rows) if gran_rows else []
    sched = perf_table(sched_rows) if sched_rows else []

    fig_speedup_eff(recs, figdir)
    fig_time_vs_n(recs, figdir)
    fig_relaxed(recs, figdir)
    fig_granularity(gran, figdir)
    fig_spacing(results, figdir, "traj_brake_seq.csv", "Leader braking - sequential baseline (11 vehicles)")
    fig_spacing(results, figdir, "traj_cruise_seq.csv", "Steady cruise - sequential baseline (11 vehicles)")
    fig_traj_diff(results, figdir)

    L = ["# Results summary", ""]
    meta = os.path.join(results, "meta.txt")
    if os.path.exists(meta):
        L += ["## Setup", "", "```", open(meta).read().strip(), "```", ""]
    L += ["Behavioural tolerances (fixed before the runs): " +
          ", ".join(f"{k} <= {v}" for k, v in TOL.items()), ""]

    # Table 1: performance, exact sync
    L += ["## 1. Computational performance (exact synchronisation)", ""]
    rows = []
    for r in recs:
        if r["sync"] != "exact":
            continue
        rows.append([r["scenario"], r["vehicles"], "seq" if r["threads"] == 0 else r["threads"],
                     f"{r['t_mean']:.4f} ± {r['t_std']:.4f}", fmt(r["speedup"], 2),
                     f"{100 * r['efficiency']:.1f} %" if r["threads"] else "-", f"{r['rtf']:.4f}"])
    L += [md_table(["scenario", "vehicles", "threads", "time [s] (mean ± sd)", "speedup",
                    "efficiency", "RTF"], rows), ""]

    # H1
    exact_omp = [r for r in recs if r["mode"] == "omp" and r["sync"] == "exact"]
    n_ident = sum(r["bit_identical"] for r in exact_omp)
    h1 = n_ident == len(exact_omp) and len(exact_omp) > 0
    L += ["## 2. Behaviour consistency (H1)", "",
          f"Exact-sync parallel configurations bit-identical to the sequential run: "
          f"**{n_ident} / {len(exact_omp)}** -> H1 **{'supported' if h1 else 'NOT supported'}**.", ""]

    # H2
    h2_ok = True
    for (scen, n), rs in group(exact_omp, ["scenario", "vehicles"]).items():
        rs = sorted(rs, key=lambda r: r["threads"])
        # judge only up to the number of performance cores where the platform is homogeneous
        if any(r["efficiency"] >= 1.0 for r in rs if r["threads"] > 1):
            h2_ok = False
    best = max(exact_omp, key=lambda r: r["speedup"]) if exact_omp else None
    L += ["## 3. Scaling (H2)", ""]
    if best:
        L += [f"Best exact-sync speedup: **{best['speedup']:.2f}x** "
              f"({best['vehicles']} vehicles, {best['threads']} threads, {best['scenario']}). "
              f"Efficiency below 100 % for all multi-thread runs: **{'yes' if h2_ok else 'no'}**. "
              "Check the speedup figure for whether speedup increases with cores "
              "(beyond the number of performance cores it may not).", ""]

    # H3: relaxed
    L += ["## 4. Relaxed synchronisation trade-off (H3)", ""]
    rows = []
    for r in recs:
        if r["sync"] != "relaxed":
            continue
        ex = next((e for e in recs if e["sync"] == "exact" and e["scenario"] == r["scenario"]
                   and e["vehicles"] == r["vehicles"] and e["threads"] == r["threads"]), None)
        rows.append([r["scenario"], r["vehicles"], r["threads"], r["K"], fmt(r["speedup"], 2),
                     fmt(ex["speedup"], 2) if ex else "-", f"{r['max_dx']:.3g}", f"{r['max_dgap']:.3g}",
                     f"{r['max_dv']:.3g}", "yes" if within_tol(r) else "**no**"])
    L += [md_table(["scenario", "vehicles", "threads", "K", "speedup (relaxed)", "speedup (exact)",
                    "max|dx| [m]", "max|dgap| [m]", "max|dv| [m/s]", "within tolerance"], rows), ""]

    # behaviour of the baseline
    L += ["## 5. Platoon behaviour (sequential baseline)", ""]
    rows = []
    for r in recs:
        if r["mode"] == "seq" and r["K"] == 1:
            rows.append([r["scenario"], r["vehicles"], fmt(r["min_gap"]), fmt(r["max_abs_err"]),
                         fmt(r["rms_err"], 4), fmt(r["string_amp"]), fmt(r["unsafe_time"], 2),
                         "yes" if r["collision"] else "no"])
    L += [md_table(["scenario", "vehicles", "min gap [m]", "max |e| [m]", "RMS e [m]",
                    "string amp.", "unsafe time [veh*s]", "collision"], rows), ""]

    if gran:
        L += ["## 6. Workload granularity (11 vehicles, braking, exact sync)", ""]
        rows = [[r["substeps"], "seq" if r["threads"] == 0 else r["threads"], f"{r['t_mean']:.4f}",
                 fmt(r["speedup"], 2), f"{100 * r['efficiency']:.1f} %" if r["threads"] else "-"]
                for r in gran]
        L += [md_table(["substeps", "threads", "time [s]", "speedup", "efficiency"], rows), ""]

    if sched:
        L += ["## 7. OpenMP schedule on heterogeneous cores (20 vehicles, braking, exact sync)", ""]
        rows = [[r["schedule"], r["threads"], f"{r['t_mean']:.4f} ± {r['t_std']:.4f}",
                 fmt(r["speedup"], 2), f"{100 * r['efficiency']:.1f} %", "yes" if r["bit_identical"] else "no"]
                for r in sched if r["mode"] == "omp"]
        L += [md_table(["schedule", "threads", "time [s]", "speedup", "efficiency", "bit-identical"], rows), ""]

    L += ["## Figures", ""] + [f"- `figures/{f}`" for f in sorted(os.listdir(figdir))]
    with open(os.path.join(results, "summary.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {results}/summary.md, {results}/summary_main.csv and {figdir}/")


if __name__ == "__main__":
    main()
