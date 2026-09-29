#!/usr/bin/env python3
"""Compare a benchmark run against a stored baseline.

Two modes:

``--expect same``
  * every sim metric must stay within ``max(k * SD_baseline, abs_floor)``
    (k = 3, per-metric absolute floors);
  * on real data the run-vs-baseline agreement must stay within the baseline
    replicate-vs-replicate agreement minus a margin (the noise floor).

``--expect improved``
  * mean sim matched accuracy over the baseline's sim datasets must increase,
    and no individual sim dataset may drop beyond its tolerance;
  * on real data ``total_admixture_rate`` must be <= baseline + tolerance for
    every dataset (skipped gracefully when the cellAdmix audit is unavailable).

Runtime and peak-RSS changes are reported; slowdowns or RSS growth beyond 20%
are warnings. Failures/timeouts and missing datasets fail the comparison.
Writes a Markdown and a JSON report; exit code 0 = pass, 1 = fail.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common    # noqa: E402
import metrics as m  # noqa: E402

# Per-metric absolute floors for --expect same (used when k*SD is smaller).
SIM_FLOORS = {
    "matched_accuracy": 0.01,
    "ari": 0.02,
    "ami": 0.02,
    "noise_precision": 0.02,
    "noise_recall": 0.02,
    "cell_count_ratio": 0.05,
    "over_segmentation_rate": 0.05,
    "under_segmentation_rate": 0.05,
    "recovery_rate": 0.05,
    "median_matched_jaccard": 0.05,
    "oracle_gap": 0.01,
}
DEFAULT_SIM_FLOOR = 0.02

# Margins subtracted from the baseline replicate agreement (higher is better).
REAL_MARGINS = {
    "molecule_ari": 0.05,
    "assigned_agreement": 0.02,
    "frac_cells_matched": 0.10,
    "median_jaccard": 0.10,
}
# Near-1 metrics: allowed absolute deviation, widened by k*SD of the
# baseline's replicate pairs when available.
REAL_RATIO_FLOORS = {
    "cell_count_ratio": 0.10,
    "median_mpc_rel_change": 0.10,
}
# Fallback absolute minima when a baseline has no replicate agreement at all.
REAL_ABS_MIN = {
    "molecule_ari": 0.90,
    "assigned_agreement": 0.95,
    "frac_cells_matched": 0.75,
    "median_jaccard": 0.75,
}

ADMIXTURE_TOLERANCE = 0.01
K_DEFAULT = 3.0
SLOWDOWN_WARN = 0.20
RSS_GROWTH_WARN = 0.20


def _finite(x) -> bool:
    try:
        return x is not None and not math.isnan(float(x))
    except (TypeError, ValueError):
        return False


def _mean(values: list[float]) -> float:
    arr = np.asarray([v for v in values if _finite(v)], dtype=float)
    return float(arr.mean()) if len(arr) else float("nan")


def _sd(values: list[float]) -> float:
    arr = np.asarray([v for v in values if _finite(v)], dtype=float)
    if len(arr) < 2:
        return 0.0
    return float(arr.std(ddof=1))


class Report:
    def __init__(self, run_id: str, baseline: str, expect: str):
        self.run_id = run_id
        self.baseline = baseline
        self.expect = expect
        self.checks: list[dict] = []
        self.warnings: list[str] = []
        self.run_failures: list[dict] = []
        self.runtime_rows: list[dict] = []
        self.meta: dict = {}

    def check(self, scope: str, dataset: str, metric: str, status: str, **kw):
        row = {"scope": scope, "dataset": dataset, "metric": metric,
               "status": status}
        row.update(kw)
        self.checks.append(row)

    def warn(self, msg: str):
        self.warnings.append(msg)

    @property
    def n_fail(self) -> int:
        return sum(1 for c in self.checks if c["status"] == "fail")

    @property
    def passed(self) -> bool:
        return self.n_fail == 0

    def to_json(self) -> dict:
        counts = {"pass": 0, "fail": 0, "skip": 0, "info": 0}
        for c in self.checks:
            counts[c["status"]] = counts.get(c["status"], 0) + 1
        return {
            "schema": 1,
            "run_id": self.run_id,
            "baseline": self.baseline,
            "expect": self.expect,
            "generated": common.utc_now(),
            "meta": self.meta,
            "checks": self.checks,
            "runtime": self.runtime_rows,
            "warnings": self.warnings,
            "run_failures": self.run_failures,
            "summary": {**counts, "pass_overall": self.passed},
        }


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------

def check_sim_dataset(rep: Report, ds_id: str, run_m: dict, base_m: dict,
                      k: float, mode: str):
    run_mean = (run_m.get("sim") or {}).get("mean") or {}
    base_mean = (base_m.get("sim") or {}).get("mean") or {}
    base_sd = (base_m.get("sim") or {}).get("sd") or {}
    if not base_mean:
        rep.check("sim", ds_id, "*", "fail", detail="baseline has no sim metrics")
        return
    if not run_mean:
        rep.check("sim", ds_id, "*", "fail", detail="run has no sim metrics (all reps failed?)")
        return

    if mode == "same":
        for metric, bval in base_mean.items():
            rval = run_mean.get(metric)
            floor = SIM_FLOORS.get(metric, DEFAULT_SIM_FLOOR)
            sd = base_sd.get(metric)
            sd = float(sd) if _finite(sd) else 0.0
            tol = max(k * sd, floor)
            if not _finite(bval) and not _finite(rval):
                rep.check("sim", ds_id, metric, "skip", detail="undefined for both")
                continue
            if not _finite(rval):
                rep.check("sim", ds_id, metric, "fail", baseline_mean=bval,
                          run_mean=rval, tolerance=tol, detail="run metric undefined")
                continue
            delta = float(rval) - float(bval)
            ok = abs(delta) <= tol
            rep.check("sim", ds_id, metric, "pass" if ok else "fail",
                      baseline_mean=float(bval), baseline_sd=sd, run_mean=float(rval),
                      delta=delta, tolerance=tol)
    else:  # improved: drop check on matched accuracy, info on everything else
        for metric, bval in base_mean.items():
            rval = run_mean.get(metric)
            if not _finite(bval) or not _finite(rval):
                rep.check("sim", ds_id, metric, "skip", detail="undefined")
                continue
            delta = float(rval) - float(bval)
            if metric == "matched_accuracy":
                sd = base_sd.get(metric)
                sd = float(sd) if _finite(sd) else 0.0
                tol = max(k * sd, SIM_FLOORS.get(metric, DEFAULT_SIM_FLOOR))
                ok = delta >= -tol
                rep.check("sim", ds_id, metric, "pass" if ok else "fail",
                          baseline_mean=float(bval), baseline_sd=sd,
                          run_mean=float(rval), delta=delta, tolerance=tol,
                          detail="drop tolerance")
            else:
                rep.check("sim", ds_id, metric, "info",
                          baseline_mean=float(bval), run_mean=float(rval), delta=delta)


def check_real_dataset(rep: Report, ds_id: str, run_m: dict, base_m: dict,
                       run_cells: list[np.ndarray], base_cells: list[np.ndarray],
                       k: float, mode: str, admixture_tolerance: float):
    if mode == "same":
        if not run_cells or not base_cells:
            rep.check("real", ds_id, "*", "fail",
                      detail="missing assignment tables for comparison")
            return
        pairs = [m.real_pair_metrics(rc, bc) for rc in run_cells for bc in base_cells]
        base_agr = ((base_m.get("real") or {}).get("rep_agreement") or {})
        b_mean = base_agr.get("mean") or {}
        b_sd = base_agr.get("sd") or {}
        has_floor = bool(b_mean)
        if not has_floor:
            rep.warn(f"{ds_id}: baseline has no replicate agreement "
                     f"(single replicate); using absolute fallbacks")
        for metric in ("molecule_ari", "assigned_agreement",
                       "frac_cells_matched", "median_jaccard"):
            r_val = _mean([p[metric] for p in pairs])
            if has_floor and _finite(b_mean.get(metric)):
                sd = b_sd.get(metric)
                sd = float(sd) if _finite(sd) else 0.0
                required = float(b_mean[metric]) - REAL_MARGINS[metric]
                rep.check("real", ds_id, metric,
                          "pass" if _finite(r_val) and r_val >= required else "fail",
                          baseline_rep_mean=float(b_mean[metric]),
                          baseline_rep_sd=sd, run_vs_baseline_mean=r_val,
                          required_min=required, margin=REAL_MARGINS[metric])
            else:
                floor = REAL_ABS_MIN[metric]
                rep.check("real", ds_id, metric,
                          "pass" if _finite(r_val) and r_val >= floor else "fail",
                          run_vs_baseline_mean=r_val, required_min=floor,
                          detail="absolute fallback (no baseline replicates)")
        for metric, abs_floor in REAL_RATIO_FLOORS.items():
            r_val = _mean([p[metric] for p in pairs])
            sd = _sd([v.get(metric) for v in (base_agr.get("per_pair") or [])])
            tol = max(abs_floor, k * sd)
            if not _finite(r_val):
                dev = float("inf")
            elif metric == "cell_count_ratio":
                dev = abs(r_val - 1.0)
            else:  # median_mpc_rel_change is already a relative deviation
                dev = abs(r_val)
            rep.check("real", ds_id, metric,
                      "pass" if dev <= tol else "fail",
                      baseline_rep_sd=sd, run_vs_baseline_mean=r_val,
                      deviation=dev, tolerance=tol,
                      detail="absolute deviation from 0 (no change)")
    else:  # improved: admixture audit only
        base_rate = ((base_m.get("real") or {}).get("celladmix") or {}).get("mean_total")
        run_rate = ((run_m.get("real") or {}).get("celladmix") or {}).get("mean_total")
        if not _finite(base_rate) or not _finite(run_rate):
            rep.check("real", ds_id, "total_admixture_rate", "skip",
                      detail="cellAdmix audit unavailable on "
                             + ("both sides" if not _finite(base_rate) and not _finite(run_rate)
                                else ("baseline" if not _finite(base_rate) else "run")))
            return
        allowed = float(base_rate) + admixture_tolerance
        ok = float(run_rate) <= allowed
        rep.check("real", ds_id, "total_admixture_rate", "pass" if ok else "fail",
                  baseline_mean=float(base_rate), run_mean=float(run_rate),
                  tolerance=admixture_tolerance, allowed_max=allowed)


def collect_run_failures(rep: Report, run_metrics: dict[str, dict]):
    for ds_id, mm in sorted(run_metrics.items()):
        for r in mm.get("reps", []):
            if r.get("status") != "ok":
                rep.run_failures.append({
                    "dataset": ds_id, "rep": r.get("rep"),
                    "status": r.get("status"), "exit_code": r.get("exit_code"),
                    "command": r.get("command_str"),
                })
                rep.check("run", ds_id, f"rep{r.get('rep')}", "fail",
                          detail=f"status={r.get('status')} "
                                 f"exit_code={r.get('exit_code')}")
                tail = (r.get("stderr_tail") or "").strip().splitlines()
                if tail:
                    rep.warn(f"{ds_id} rep{r.get('rep')} stderr tail: {tail[-1]}")


def compare_runtime(rep: Report, run_metrics: dict[str, dict],
                    base_metrics: dict[str, dict]):
    for ds_id in sorted(set(run_metrics) & set(base_metrics)):
        r = run_metrics[ds_id].get("runtime") or {}
        b = base_metrics[ds_id].get("runtime") or {}
        row = {"dataset": ds_id,
               "wall_base_s": b.get("wall_s_mean"), "wall_run_s": r.get("wall_s_mean"),
               "rss_base_kb": b.get("peak_rss_kb_mean"),
               "rss_run_kb": r.get("peak_rss_kb_mean")}
        row["wall_change"] = _ratio_change(r.get("wall_s_mean"), b.get("wall_s_mean"))
        row["rss_change"] = _ratio_change(r.get("peak_rss_kb_mean"),
                                          b.get("peak_rss_kb_mean"))
        rep.runtime_rows.append(row)
        if _finite(row["wall_change"]) and row["wall_change"] > SLOWDOWN_WARN:
            rep.warn(f"{ds_id}: runtime +{row['wall_change'] * 100:.1f}% "
                     f"({row['wall_base_s']:.2f}s -> {row['wall_run_s']:.2f}s)")
        if _finite(row["rss_change"]) and row["rss_change"] > RSS_GROWTH_WARN:
            rep.warn(f"{ds_id}: peak RSS +{row['rss_change'] * 100:.1f}% "
                     f"({row['rss_base_kb']:.0f}kB -> {row['rss_run_kb']:.0f}kB)")


def _ratio_change(new, base) -> float:
    if not _finite(new) or not _finite(base) or float(base) == 0:
        return float("nan")
    return float(new) / float(base) - 1.0


def check_improved_aggregate(rep: Report, run_metrics: dict[str, dict],
                             base_metrics: dict[str, dict]):
    """Mode 'improved': mean sim matched accuracy must increase."""
    sim_ids = [d for d, mm in base_metrics.items()
               if mm.get("kind") == "sim" or mm.get("dataset", {}).get("kind") == "sim"]
    if not sim_ids:
        rep.check("sim", "*", "aggregate_matched_accuracy", "fail",
                  detail="baseline contains no sim datasets")
        return

    def agg(metrics):
        vals = []
        for d in sim_ids:
            mm = metrics.get(d) or {}
            v = ((mm.get("sim") or {}).get("mean") or {}).get("matched_accuracy")
            if _finite(v):
                vals.append(float(v))
        return _mean(vals) if vals else float("nan")

    run_agg = agg(run_metrics)
    base_agg = agg(base_metrics)
    ok = _finite(run_agg) and _finite(base_agg) and run_agg > base_agg
    rep.check("sim", "*", "aggregate_matched_accuracy", "pass" if ok else "fail",
              baseline_mean=base_agg, run_mean=run_agg,
              delta=(run_agg - base_agg) if _finite(run_agg) and _finite(base_agg)
              else None,
              detail="mean matched accuracy over baseline sim datasets must increase")


# ---------------------------------------------------------------------------
# reports
# ---------------------------------------------------------------------------

def _fmt(v, nd=4) -> str:
    if v is None:
        return "-"
    if isinstance(v, float):
        if math.isnan(v):
            return "nan"
        if abs(v) >= 1000:
            return f"{v:.1f}"
        return f"{v:.{nd}f}"
    return str(v)


def render_markdown(rep: Report) -> str:
    lines = []
    lines.append(f"# Benchmark comparison: run `{rep.run_id}` vs baseline "
                 f"`{rep.baseline}` (expect: **{rep.expect}**)")
    lines.append("")
    lines.append(f"Generated: {common.utc_now()}")
    for k, v in rep.meta.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    verdict = "**PASS**" if rep.passed else "**FAIL**"
    counts = {"pass": 0, "fail": 0, "skip": 0, "info": 0}
    for c in rep.checks:
        counts[c["status"]] = counts.get(c["status"], 0) + 1
    lines.append(f"## Summary: {verdict} — "
                 f"{counts['pass']} passed, {counts['fail']} failed, "
                 f"{counts['skip']} skipped, {counts['info']} info; "
                 f"{len(rep.warnings)} warning(s)")
    lines.append("")

    def section(title, scope):
        rows = [c for c in rep.checks if c["scope"] == scope]
        if not rows:
            return
        lines.append(f"## {title}")
        lines.append("")
        lines.append("| dataset | metric | baseline | run | delta | tol/req | status |")
        lines.append("|---|---|---|---|---|---|---|")
        for c in rows:
            b = c.get("baseline_mean", c.get("baseline_rep_mean"))
            r = c.get("run_mean", c.get("run_vs_baseline_mean"))
            tol = c.get("tolerance", c.get("required_min"))
            delta = c.get("delta", c.get("deviation"))
            status = c["status"].upper()
            if c.get("detail") and c["status"] in ("fail", "skip"):
                status += f" ({c['detail']})"
            lines.append(f"| {c['dataset']} | {c['metric']} | {_fmt(b)} | {_fmt(r)} "
                         f"| {_fmt(delta)} | {_fmt(tol)} | {status} |")
        lines.append("")

    section("Sim datasets (vs truth)", "sim")
    section("Real datasets (vs baseline segmentation / admixture)", "real")
    section("Run health", "run")

    if rep.runtime_rows:
        lines.append("## Runtime & memory")
        lines.append("")
        lines.append("| dataset | wall base (s) | wall run (s) | Δwall | "
                     "RSS base (kB) | RSS run (kB) | ΔRSS |")
        lines.append("|---|---|---|---|---|---|---|")
        for row in rep.runtime_rows:
            wc, rc = row["wall_change"], row["rss_change"]
            wtxt = f"{wc * 100:+.1f}%" if _finite(wc) else "-"
            rtxt = f"{rc * 100:+.1f}%" if _finite(rc) else "-"
            lines.append(f"| {row['dataset']} | {_fmt(row['wall_base_s'], 2)} "
                         f"| {_fmt(row['wall_run_s'], 2)} | {wtxt} "
                         f"| {_fmt(row['rss_base_kb'], 0)} "
                         f"| {_fmt(row['rss_run_kb'], 0)} | {rtxt} |")
        lines.append("")

    if rep.run_failures:
        lines.append("## Failures / timeouts")
        lines.append("")
        for f in rep.run_failures:
            lines.append(f"- {f['dataset']} rep{f['rep']}: {f['status']} "
                         f"(exit {f['exit_code']}) `{f['command']}`")
        lines.append("")

    if rep.warnings:
        lines.append("## Warnings")
        lines.append("")
        for w in rep.warnings:
            lines.append(f"- {w}")
        lines.append("")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def load_assignments(paths: list[Path]) -> list[np.ndarray]:
    return [common.assignment_cells(p) for p in paths]


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--expect", required=True, choices=["same", "improved"])
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--baselines-dir", default=None,
                    help="default <repo>/benchmarks/baselines")
    ap.add_argument("--k", type=float, default=K_DEFAULT,
                    help="tolerance multiplier on baseline SD (default 3)")
    ap.add_argument("--admixture-tolerance", type=float,
                    default=ADMIXTURE_TOLERANCE,
                    help="allowed increase of total_admixture_rate (default 0.01)")
    ap.add_argument("--report-md", default=None)
    ap.add_argument("--report-json", default=None)
    args = ap.parse_args(argv)

    repo = common.repo_root()
    root = common.data_root(args.data_root)
    baselines_dir = Path(args.baselines_dir) if args.baselines_dir \
        else repo / "benchmarks" / "baselines"
    run_root = root / "runs" / args.run_id
    base_dir = baselines_dir / args.baseline

    if not run_root.is_dir():
        print(f"error: run '{args.run_id}' not found under {run_root}", file=sys.stderr)
        return 2
    if not base_dir.is_dir():
        print(f"error: baseline '{args.baseline}' not found at {base_dir}",
              file=sys.stderr)
        return 2

    run_metrics = {p.parent.name: common.read_json(p)
                   for p in sorted(run_root.glob("*/metrics.json"))}
    base_metrics = {p.stem: common.read_json(p)
                    for p in sorted(base_dir.glob("*.json"))}
    if not run_metrics:
        print(f"error: no metrics.json under {run_root}", file=sys.stderr)
        return 2
    if not base_metrics:
        print(f"error: baseline '{args.baseline}' has no metrics", file=sys.stderr)
        return 2

    rep = Report(args.run_id, args.baseline, args.expect)

    # provenance / metadata
    any_run = next(iter(run_metrics.values()))
    any_base = next(iter(base_metrics.values()))
    run_sha = (any_run.get("binary") or {}).get("sha256")
    base_sha = (any_base.get("binary") or {}).get("sha256")
    rep.meta = {
        "run label": any_run.get("label"),
        "run binary sha256": run_sha,
        "baseline binary sha256": base_sha,
        "k": args.k,
    }
    if run_sha != base_sha:
        rep.warn("binary sha256 differs from baseline "
                 "(expected for algorithm changes, suspicious for 'same')")
    if any_run.get("threads") != any_base.get("threads"):
        rep.warn(f"thread count differs: run={any_run.get('threads')} "
                 f"baseline={any_base.get('threads')} (noise floor may not transfer)")

    # run health
    collect_run_failures(rep, run_metrics)

    # dataset coverage
    for ds_id in sorted(set(base_metrics) - set(run_metrics)):
        rep.check("run", ds_id, "presence", "fail",
                  detail="dataset in baseline but missing from run")
    for ds_id in sorted(set(run_metrics) - set(base_metrics)):
        rep.warn(f"{ds_id}: present in run but not in baseline (ignored)")

    # per-dataset checks
    for ds_id in sorted(set(run_metrics) & set(base_metrics)):
        run_m, base_m = run_metrics[ds_id], base_metrics[ds_id]
        kind = (run_m.get("dataset") or {}).get("kind") or run_m.get("kind")
        if kind == "sim":
            check_sim_dataset(rep, ds_id, run_m, base_m, args.k, args.expect)
        elif kind == "real":
            run_cells, base_cells = [], []
            if args.expect == "same":
                run_paths = [run_root / ds_id / f"rep{r['rep']}" / "assignment.parquet"
                             for r in run_m.get("reps", [])
                             if r.get("status") == "ok" and r.get("assignment")]
                base_assign = root / "baselines" / args.baseline / ds_id
                base_paths = sorted(base_assign.glob("rep*/assignment.parquet")) \
                    if base_assign.is_dir() else []
                if not base_paths:
                    rep.check("real", ds_id, "baseline_assignments", "fail",
                              detail=f"missing {base_assign} (recreate the baseline)")
                run_cells = load_assignments(run_paths)
                base_cells = load_assignments(base_paths)
            check_real_dataset(rep, ds_id, run_m, base_m, run_cells, base_cells,
                               args.k, args.expect, args.admixture_tolerance)

    if args.expect == "improved":
        check_improved_aggregate(rep, run_metrics, base_metrics)

    compare_runtime(rep, run_metrics, base_metrics)

    # reports
    md_path = Path(args.report_md) if args.report_md else \
        run_root / f"compare_{args.baseline}_{args.expect}.md"
    json_path = Path(args.report_json) if args.report_json else \
        run_root / f"compare_{args.baseline}_{args.expect}.json"
    md = render_markdown(rep)
    common.write_json(json_path, rep.to_json())
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text(md)

    print(md)
    print(f"report: {md_path}\n        {json_path}")
    return 0 if rep.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
