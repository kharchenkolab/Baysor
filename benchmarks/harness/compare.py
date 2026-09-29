#!/usr/bin/env python3
"""Compare a benchmark run against a stored baseline.

Three expectation modes:

``--expect identical`` (default; the refactor gate)
  Run and baseline must both be 1-thread with exactly 1 replicate. For every
  dataset the ``assignment_sha256`` recorded per replicate in ``metrics.json``
  must match; on a mismatch the comparison fails and reports the metric
  deltas. Baysor is bitwise-deterministic at 1 thread, so a no-behaviour-change
  refactor must reproduce the baseline assignments exactly.

``--expect same``
  The algorithm must not have changed beyond the measurement noise:
  * provenance gates fail on thread-count or binary/threads mismatches and on
    differing dataset content hashes (``inputs.molecules_sha256`` /
    ``inputs.meta_sha256``);
  * a real dataset in the baseline with fewer than 2 successful replicates
    is an error (exit 2): replicate-agreement gates cannot be evaluated
    against a single-segmentation baseline and no thresholds are invented;
  * only a few **primary** metrics gate the verdict —
    sim: one-to-one accuracy, ARI over assigned molecules, recovery fraction,
    cell-count ratio; real: ARI over assigned molecules, matched-cell
    fraction, cell-count ratio. All other metrics are informational;
  * tolerance is ``max(k * SD_pooled, floor)`` with k = 3 and the SD pooled
    per metric across the baseline's datasets of the same kind (3-replicate
    SDs alone are too unreliable), floors calibrated in the README;
  * a false-alarm budget (normal approximation) is reported with the gated
    checks.

``--expect improved``
  Improvement must exceed the noise:
  * the mean gain in primary sim accuracy (one-to-one) across sim datasets
    must exceed both 2 standard errors (from the per-dataset replicate SDs)
    and a minimum effect (0.005); when the baseline or the run has no sim
    datasets this gate is ``skip``;
  * no sim dataset may regress beyond its tolerance in one-to-one accuracy,
    assigned-ARI, recovery or the over-segmentation rate;
  * on real data the cellAdmix ``total_admixture_rate`` must be
    <= baseline + max(k*SD(baseline audit replicates), floor), floored by
    ``--admixture-tolerance`` (default 0.0025, calibrated from
    ``celladmix/results/harness_baysor_sd.json``), gated only when the audit
    status is ``ok`` on both sides, the dataset is ``admixture_capable`` and
    the baseline crop has >= 2000 cells; with fewer than 2 baseline audit
    replicates the floor alone is used (with a warning). Otherwise it is
    reported as ``unavailable`` (never as 0);
  * when neither gate evaluates anything, the comparison fails with a
    ``nothing to evaluate`` message.

Runtime and peak-RSS changes are reported; slowdowns or RSS growth beyond
20% are warnings. Failures/timeouts and missing datasets fail the comparison.
Writes a Markdown and a JSON report; exit code 0 = pass, 1 = fail, 2 = setup
error.
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

MODES = ("identical", "same", "improved")

# --- primary (gated) metrics ------------------------------------------------
SIM_PRIMARY_SAME = ("accuracy_1to1", "ari_assigned", "recovery_rate",
                    "cell_count_ratio")
SIM_REGRESSION_IMPROVED = ("accuracy_1to1", "ari_assigned", "recovery_rate",
                           "over_segmentation_rate")
REAL_PRIMARY_SAME = ("ari_assigned", "frac_cells_matched",
                     "cell_count_ratio")
# informational metrics (never gated)
REAL_INFO = ("molecule_ari", "assigned_agreement", "assigned_fraction_candidate",
             "assigned_fraction_reference", "noise_precision", "noise_recall",
             "median_jaccard", "median_mpc_rel_change",
             "n_cells_reference", "n_cells_candidate")

# --- floors for tolerance = max(k * SD_pooled, floor) ------------------------
# Calibrated with recompute_metrics.py on harness-val{1..4} and
# rev-{same,scale09}-{sim,real}; see benchmarks/harness/README.md.
SIM_FLOORS = {
    "accuracy_1to1": 0.01,
    "ari_assigned": 0.02,
    "recovery_rate": 0.05,
    "cell_count_ratio": 0.03,
    "over_segmentation_rate": 0.05,
}
DEFAULT_SIM_FLOOR = 0.02
REAL_FLOORS = {
    "ari_assigned": 0.02,       # >= 0.02 by design (calibration: scale x0.9
    "frac_cells_matched": 0.05, # drops are caught without false alarms)
    "cell_count_ratio": 0.03,
}
DEFAULT_REAL_FLOOR = 0.05

K_DEFAULT = 3.0
MIN_EFFECT_ACCURACY = 0.005     # --expect improved: minimum mean gain
ADMIXTURE_MIN_CELLS = 2000      # admixture gate only on crops with >= cells
# calibrated floor for --admixture-tolerance: 3 x SD of the Baysor replicate
# audit on xenium_lung_cancer_quick (0.000844, celladmix/results/
# harness_baysor_sd.json); also covers baselines with < 2 audit replicates
DEFAULT_ADMIXTURE_TOLERANCE = 0.0025
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


def _pooled_sd(pairs: list[tuple[float, int]]) -> float:
    """Within-dataset pooled SD: sqrt(sum((n-1)*sd^2) / sum(n-1))."""
    num = sum((n - 1) * s * s for s, n in pairs
              if n > 1 and _finite(s) and s >= 0)
    den = sum(n - 1 for s, n in pairs if n > 1 and _finite(s) and s >= 0)
    return math.sqrt(num / den) if den > 0 else 0.0


def _two_sided_alpha(tol: float, sd: float) -> float:
    """Two-sided normal tail probability beyond ``tol`` given ``sd``."""
    if not _finite(sd) or sd <= 0 or not _finite(tol) or tol <= 0:
        return 0.0
    z = tol / sd
    return math.erfc(z / math.sqrt(2))     # = 2 * (1 - Phi(z))


def _one_sided_alpha(tol: float, sd: float) -> float:
    if not _finite(sd) or sd <= 0 or not _finite(tol) or tol <= 0:
        return 0.0
    z = tol / sd
    return 0.5 * math.erfc(z / math.sqrt(2))   # = 1 - Phi(z)


def pooled_sd_by_metric(base_metrics: dict[str, dict], kind: str) -> dict[str, float]:
    """Pool per-metric SDs across the baseline's datasets of one kind.

    Sim: replicate SDs of the dataset means (``sim.sd``, n = replicate
    count). Real: replicate-pair agreement SDs (``real.rep_agreement.sd``,
    n = pair count). Pooling guards against the unreliability of a single
    3-replicate SD.
    """
    per_metric: dict[str, list[tuple[float, int]]] = {}
    for mm in base_metrics.values():
        ds_kind = (mm.get("dataset") or {}).get("kind") or mm.get("kind")
        if ds_kind != kind:
            continue
        if kind == "sim":
            block = mm.get("sim") or {}
            sd = block.get("sd") or {}
            n = int(block.get("n_metric_reps") or 0)
        else:
            block = ((mm.get("real") or {}).get("rep_agreement")) or {}
            sd = block.get("sd") or {}
            n = len(block.get("per_pair") or [])
        if n < 2:
            continue
        for key, val in sd.items():
            if _finite(val):
                per_metric.setdefault(key, []).append((float(val), n))
    return {key: _pooled_sd(v) for key, v in per_metric.items()}


class Report:
    def __init__(self, run_id: str, baseline: str, expect: str):
        self.run_id = run_id
        self.baseline = baseline
        self.expect = expect
        self.checks: list[dict] = []
        self.warnings: list[str] = []
        self.run_failures: list[dict] = []
        self.runtime_rows: list[dict] = []
        self.gate_info: list[dict] = []
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

    @property
    def gated(self) -> list[dict]:
        return [c for c in self.checks if "false_alarm_p" in c]

    def finalize_budget(self, k: float):
        """Record the false-alarm budget of the gated (primary) checks."""
        gated = self.gated
        if not gated:
            return
        total = sum(float(c.get("false_alarm_p") or 0.0) for c in gated)
        self.meta["gated checks"] = len(gated)
        self.meta["false-alarm budget"] = (
            f"~{total:.3f} expected false failures across {len(gated)} "
            f"gated checks (normal approximation, k={k:g})")
        if total > 0.5:
            self.warn(f"false-alarm budget is high (~{total:.2f} expected "
                      f"false failures); consider a larger floor or k")

    def to_json(self) -> dict:
        counts = {"pass": 0, "fail": 0, "skip": 0, "info": 0}
        for c in self.checks:
            counts[c["status"]] = counts.get(c["status"], 0) + 1
        return {
            "schema": 2,
            "run_id": self.run_id,
            "baseline": self.baseline,
            "expect": self.expect,
            "generated": common.utc_now(),
            "meta": self.meta,
            "gates": self.gate_info,
            "checks": self.checks,
            "runtime": self.runtime_rows,
            "warnings": self.warnings,
            "run_failures": self.run_failures,
            "summary": {**counts, "pass_overall": self.passed},
        }


# ---------------------------------------------------------------------------
# provenance / content checks
# ---------------------------------------------------------------------------

def check_provenance(rep: Report, ds_id: str, run_m: dict, base_m: dict,
                     mode: str):
    """Threads / binary sha per dataset (and identical-mode preconditions)."""
    if mode == "identical":
        for side, mm in (("run", run_m), ("baseline", base_m)):
            ok_threads = mm.get("threads") == 1
            rep.check("provenance", ds_id, f"{side}_threads_1",
                      "pass" if ok_threads else "fail",
                      threads=mm.get("threads"),
                      detail="identical requires threads == 1"
                             if not ok_threads else None)
            ok_reps = mm.get("replicates") == 1
            rep.check("provenance", ds_id, f"{side}_replicates_1",
                      "pass" if ok_reps else "fail",
                      replicates=mm.get("replicates"),
                      detail="identical requires exactly 1 replicate"
                             if not ok_reps else None)
        return
    # same / improved: a thread mismatch invalidates the noise floor;
    # in 'same' a binary mismatch means the comparison is not an unchanged run
    r_threads, b_threads = run_m.get("threads"), base_m.get("threads")
    same_threads = r_threads == b_threads
    rep.check("provenance", ds_id, "threads", "pass" if same_threads else "fail",
              run_value=r_threads, baseline_value=b_threads,
              detail=None if same_threads else
              "thread count differs; the baseline noise floor does not transfer")
    if mode == "same":
        r_sha = (run_m.get("binary") or {}).get("sha256")
        b_sha = (base_m.get("binary") or {}).get("sha256")
        same_sha = r_sha is not None and r_sha == b_sha
        rep.check("provenance", ds_id, "binary_sha256",
                  "pass" if same_sha else "fail",
                  run_value=r_sha, baseline_value=b_sha,
                  detail=None if same_sha else
                  "binary sha256 differs from the baseline ('same' requires "
                  "the unchanged binary)")


def check_content_hashes(rep: Report, ds_id: str, run_m: dict, base_m: dict):
    """Dataset content hashes recorded by the runner (``inputs``).

    A mismatch fails in every mode: comparing against different input data
    is meaningless. Runs/baselines created before the hashes were recorded
    produce a ``skip`` (with a hint) instead.
    """
    ri = run_m.get("inputs") or {}
    bi = base_m.get("inputs") or {}
    for key in ("molecules_sha256", "meta_sha256"):
        rv, bv = ri.get(key), bi.get(key)
        if rv and bv:
            ok = rv == bv
            rep.check("provenance", ds_id, key, "pass" if ok else "fail",
                      run_value=rv, baseline_value=bv,
                      detail=None if ok else "dataset content differs")
        elif rv or bv:
            rep.check("provenance", ds_id, key, "skip",
                      detail="content hash recorded on one side only "
                             "(regenerate the baseline with recompute_metrics.py)")
            rep.warn(f"{ds_id}: {key} recorded on one side only")
        else:
            rep.check("provenance", ds_id, key, "skip",
                      detail="content hash not recorded (older metrics.json; "
                             "rerun with recompute_metrics.py)")


def _metric_delta_rows(rep: Report, scope: str, ds_id: str, kind: str,
                       run_m: dict, base_m: dict, detail: str):
    """Informational metric-delta rows (used by 'identical' on sha mismatch)."""
    if kind == "sim":
        run_mean = ((run_m.get("sim") or {}).get("mean")) or {}
        base_mean = ((base_m.get("sim") or {}).get("mean")) or {}
        for key in sorted(set(run_mean) | set(base_mean)):
            rv, bv = run_mean.get(key), base_mean.get(key)
            delta = (float(rv) - float(bv)
                     if _finite(rv) and _finite(bv) else None)
            rep.check(scope, ds_id, key, "info",
                      baseline_mean=bv, run_mean=rv, delta=delta, detail=detail)
    else:
        run_agr = ((run_m.get("real") or {}).get("rep_agreement") or {}).get("mean") or {}
        base_agr = ((base_m.get("real") or {}).get("rep_agreement") or {}).get("mean") or {}
        for key in sorted(set(run_agr) | set(base_agr)):
            rv, bv = run_agr.get(key), base_agr.get(key)
            delta = (float(rv) - float(bv)
                     if _finite(rv) and _finite(bv) else None)
            rep.check(scope, ds_id, key, "info",
                      baseline_rep_mean=bv, run_vs_baseline_mean=rv,
                      delta=delta, detail=detail)


# ---------------------------------------------------------------------------
# identical mode
# ---------------------------------------------------------------------------

def check_identical_dataset(rep: Report, ds_id: str, run_m: dict, base_m: dict,
                            run_cells: list[np.ndarray],
                            base_cells: list[np.ndarray]):
    """Compare per-replicate assignment sha256; emit metric deltas on mismatch."""
    run_reps = [r for r in run_m.get("reps", [])
                if r.get("status") == "ok" and r.get("assignment_sha256")]
    base_reps = [r for r in base_m.get("reps", [])
                 if r.get("status") == "ok" and r.get("assignment_sha256")]
    kind = (run_m.get("dataset") or {}).get("kind") or run_m.get("kind")
    detail = "assignment sha256 mismatch; metric deltas below"
    mismatch = False
    if not run_reps or not base_reps:
        rep.check("identical", ds_id, "assignment_sha256", "fail",
                  detail="no successful replicate with an assignment sha "
                         "on the " + ("run side" if not run_reps else "baseline side"))
        mismatch = True
    elif len(run_reps) != len(base_reps):
        rep.check("identical", ds_id, "assignment_sha256", "fail",
                  detail=f"replicate counts differ ({len(run_reps)} vs "
                         f"{len(base_reps)})")
        mismatch = True
    else:
        for rb, bb in zip(run_reps, base_reps):
            ok = rb["assignment_sha256"] == bb["assignment_sha256"]
            mismatch = mismatch or not ok
            rep.check("identical", ds_id,
                      f"rep{rb['rep']}_assignment_sha256",
                      "pass" if ok else "fail",
                      run_value=rb["assignment_sha256"],
                      baseline_value=bb["assignment_sha256"],
                      detail=None if ok else detail)
    if not mismatch:
        return
    _metric_delta_rows(rep, "identical", ds_id, kind, run_m, base_m, detail)
    if kind == "real":
        if run_cells and base_cells:
            pair = m.real_pair_metrics(run_cells[0], base_cells[0])
            for key, val in pair.items():
                rep.check("identical", ds_id, f"run_vs_baseline_{key}", "info",
                          run_vs_baseline=val, detail=detail)
        else:
            rep.warn(f"{ds_id}: assignment tables unavailable for pair "
                     f"metrics on sha mismatch")


# ---------------------------------------------------------------------------
# same mode
# ---------------------------------------------------------------------------

def check_sim_dataset(rep: Report, ds_id: str, run_m: dict, base_m: dict,
                      k: float, pooled: dict[str, float], mode: str = "same"):
    run_mean = (run_m.get("sim") or {}).get("mean") or {}
    base_mean = (base_m.get("sim") or {}).get("mean") or {}
    base_sd = (base_m.get("sim") or {}).get("sd") or {}
    if not base_mean:
        rep.check("sim", ds_id, "*", "fail", detail="baseline has no sim metrics")
        return
    if not run_mean:
        rep.check("sim", ds_id, "*", "fail",
                  detail="run has no sim metrics (all reps failed?)")
        return

    gated = SIM_REGRESSION_IMPROVED if mode == "improved" else SIM_PRIMARY_SAME
    floors = SIM_FLOORS
    lower_better = {"over_segmentation_rate"}
    for metric in gated:
        if metric not in base_mean:
            rep.check("sim", ds_id, metric, "skip",
                      detail="metric not in baseline (regenerate it with "
                             "recompute_metrics.py)")
            rep.warn(f"{ds_id}: baseline lacks primary metric {metric!r}")
            continue
        bval = base_mean[metric]
        rval = run_mean.get(metric)
        sd = pooled.get(metric, 0.0)
        floor = floors.get(metric, DEFAULT_SIM_FLOOR)
        tol = max(k * sd, floor)
        if not _finite(rval):
            rep.check("sim", ds_id, metric, "fail", baseline_mean=bval,
                      run_mean=rval, tolerance=tol,
                      detail="run metric undefined")
            continue
        if not _finite(bval):
            rep.check("sim", ds_id, metric, "skip", detail="undefined in baseline")
            continue
        delta = float(rval) - float(bval)
        sd_base = base_sd.get(metric)
        sd_base = float(sd_base) if _finite(sd_base) else 0.0
        if mode == "improved":
            if metric in lower_better:
                ok = delta <= tol
            else:
                ok = delta >= -tol
            alpha = _one_sided_alpha(tol, sd)
            detail = "regression beyond tolerance"
        else:
            ok = abs(delta) <= tol
            alpha = _two_sided_alpha(tol, sd)
            detail = None
        rep.check("sim", ds_id, metric, "pass" if ok else "fail",
                  baseline_mean=float(bval), baseline_sd=sd_base,
                  run_mean=float(rval), delta=delta, tolerance=tol,
                  false_alarm_p=alpha,
                  detail=None if ok else detail)

    # everything else is informational
    for metric in sorted(set(run_mean) | set(base_mean)):
        if metric in gated:
            continue
        bval, rval = base_mean.get(metric), run_mean.get(metric)
        delta = (float(rval) - float(bval)
                 if _finite(rval) and _finite(bval) else None)
        rep.check("sim", ds_id, metric, "info",
                  baseline_mean=bval, run_mean=rval, delta=delta)


def load_pair_cells(run_root: Path, root: Path, baseline: str, ds_id: str,
                    run_m: dict) -> tuple[list[np.ndarray], list[np.ndarray]]:
    run_paths = [run_root / ds_id / f"rep{r['rep']}" / "assignment.parquet"
                 for r in run_m.get("reps", [])
                 if r.get("status") == "ok" and r.get("assignment")]
    base_assign = root / "baselines" / baseline / ds_id
    base_paths = sorted(base_assign.glob("rep*/assignment.parquet")) \
        if base_assign.is_dir() else []
    return load_assignments(run_paths), load_assignments(base_paths)


def check_real_dataset(rep: Report, ds_id: str, run_m: dict, base_m: dict,
                       run_cells: list[np.ndarray], base_cells: list[np.ndarray],
                       k: float, mode: str, pooled: dict[str, float],
                       admixture_tolerance_floor: float):
    if mode == "same":
        if not run_cells or not base_cells:
            rep.check("real", ds_id, "*", "fail",
                      detail="missing assignment tables for comparison")
            return
        pairs = [m.real_pair_metrics(rc, bc) for rc in run_cells for bc in base_cells]
        base_agr = ((base_m.get("real") or {}).get("rep_agreement") or {})
        b_mean = base_agr.get("mean") or {}
        has_floor = bool(b_mean)

        def gated_row(metric, r_val, center, tol, ok, deviation, *,
                      two_sided=True, **extra):
            sd = pooled.get(metric, 0.0)
            alpha = (_two_sided_alpha(tol, sd) if two_sided
                     else _one_sided_alpha(tol, sd))
            rep.check("real", ds_id, metric, "pass" if ok else "fail",
                      baseline_rep_mean=center, run_vs_baseline_mean=r_val,
                      deviation=deviation, tolerance=tol,
                      false_alarm_p=alpha, **extra)

        for metric in REAL_PRIMARY_SAME:
            r_val = _mean([p.get(metric) for p in pairs])
            floor = REAL_FLOORS.get(metric, DEFAULT_REAL_FLOOR)
            sd = pooled.get(metric, 0.0)
            tol = max(k * sd, floor)
            if not _finite(r_val):
                rep.check("real", ds_id, metric, "fail",
                          run_vs_baseline_mean=r_val, tolerance=tol,
                          detail="run metric undefined")
                continue
            if metric == "cell_count_ratio":
                # run and baseline must segment about the same number of
                # cells: deviation of the mean ratio from 1 (two-sided —
                # both over- and under-segmentation are regressions)
                dev = abs(r_val - 1.0)
                gated_row(metric, r_val, 1.0, tol, dev <= tol, dev,
                          detail="|ratio - 1| vs tolerance")
            elif has_floor and _finite(b_mean.get(metric)):
                center = float(b_mean[metric])
                # one-sided: the run must not agree *less* with the baseline
                # than the baseline's own replicates agree with each other.
                # Agreement above that level is fine — a self-comparison
                # (run = the baseline's source) contains identity pairs and
                # systematically exceeds the centre, which a two-sided gate
                # would fail.
                dev = center - r_val
                gated_row(metric, r_val, center, tol, dev <= tol, dev,
                          two_sided=False,
                          detail="run-vs-baseline agreement vs baseline "
                                 "replicate agreement (one-sided)")
            else:
                # unreachable: main() rejects same-mode comparisons whose
                # real baseline has < 2 successful replicates (exit 2)
                rep.check("real", ds_id, metric, "fail",
                          run_vs_baseline_mean=r_val,
                          detail="baseline has no replicate agreement "
                                 "(needs >=2 successful replicates; should "
                                 "have been rejected with exit 2)")

        for metric in REAL_INFO:
            if metric not in pairs[0]:
                continue
            r_val = _mean([p.get(metric) for p in pairs])
            rep.check("real", ds_id, metric, "info", run_vs_baseline_mean=r_val)
    else:  # improved: admixture audit only
        check_real_admixture(rep, ds_id, run_m, base_m, k,
                             admixture_tolerance_floor)


def check_real_admixture(rep: Report, ds_id: str, run_m: dict, base_m: dict,
                         k: float, tolerance_floor: float):
    base_cam = ((base_m.get("real") or {}).get("celladmix")) or {}
    run_cam = ((run_m.get("real") or {}).get("celladmix")) or {}
    base_rate = base_cam.get("mean_total")
    run_rate = run_cam.get("mean_total")

    def unavailable(reason: str):
        rep.check("real", ds_id, "total_admixture_rate", "skip",
                  detail=f"unavailable: {reason}")

    if base_cam.get("status") != "ok" or run_cam.get("status") != "ok":
        which = []
        if run_cam.get("status") != "ok":
            which.append(f"run={run_cam.get('status', 'missing')}")
        if base_cam.get("status") != "ok":
            which.append(f"baseline={base_cam.get('status', 'missing')}")
        unavailable("cellAdmix audit not ok (" + ", ".join(which) + ")")
        return
    if not _finite(base_rate) or not _finite(run_rate):
        unavailable("audit rate missing on "
                    + ("both sides" if not _finite(base_rate) and not _finite(run_rate)
                       else ("baseline" if not _finite(base_rate) else "run")))
        return
    # the audit must be able to detect admixture at all (FIX-A2's flag;
    # falls back to the cell-count rule when the flag was not recorded)
    if base_cam.get("admixture_capable") is False \
            or run_cam.get("admixture_capable") is False:
        unavailable("dataset is not admixture_capable "
                    f"(crop < {ADMIXTURE_MIN_CELLS} cells)")
        return
    # gate only on crops large enough for the audit to have statistical power
    base_cells = [r.get("n_cells") for r in base_m.get("reps", [])
                  if r.get("status") == "ok" and _finite(r.get("n_cells"))]
    n_cells = int(round(_mean(base_cells))) if base_cells else 0
    if n_cells < ADMIXTURE_MIN_CELLS:
        unavailable(f"baseline crop has {n_cells} cells "
                    f"(< {ADMIXTURE_MIN_CELLS})")
        return
    rates = [t for t in (base_cam.get("per_rep_total") or [])
             if _finite(t)]
    if len(rates) >= 2:
        sd = _sd(rates)
        tol = max(k * sd, tolerance_floor)
    else:
        # <2 baseline audit replicates: no SD can be estimated — use the
        # calibrated floor only, and say so
        sd = None
        tol = tolerance_floor
        rep.warn(f"{ds_id}: baseline audit has {len(rates)} replicate(s) "
                 f"(<2); using --admixture-tolerance floor "
                 f"{tolerance_floor:g}")
    allowed = float(base_rate) + tol
    ok = float(run_rate) <= allowed
    rep.check("real", ds_id, "total_admixture_rate", "pass" if ok else "fail",
              baseline_mean=float(base_rate), run_mean=float(run_rate),
              baseline_sd=sd, tolerance=tol, allowed_max=allowed,
              n_cells=n_cells, false_alarm_p=_one_sided_alpha(tol, sd or 0.0),
              detail=(None if ok else
                      "admixture rate rose beyond tolerance "
                      "(max(k*SD of baseline audit replicates, floor))"))


# ---------------------------------------------------------------------------
# improved: aggregate gain
# ---------------------------------------------------------------------------

def check_improved_aggregate(rep: Report, run_metrics: dict[str, dict],
                             base_metrics: dict[str, dict], k: float):
    """Mean one-to-one accuracy gain must exceed 2 SE and a minimum effect.

    Skips (rather than fails) when either side has no sim datasets — e.g. a
    real-only baseline whose improvement is judged on the admixture gate.
    """
    def sim_ids(metrics):
        return [d for d, mm in metrics.items()
                if (mm.get("dataset") or {}).get("kind") == "sim"
                or mm.get("kind") == "sim"]

    base_sim = sim_ids(base_metrics)
    run_sim = sim_ids(run_metrics)
    if not base_sim:
        rep.check("sim", "*", "aggregate_accuracy_1to1", "skip",
                  detail="baseline contains no sim datasets")
        return
    if not run_sim:
        rep.check("sim", "*", "aggregate_accuracy_1to1", "skip",
                  detail="run contains no sim datasets")
        return
    sim_ids_eval = [d for d in base_sim if d in run_sim]

    def stats(metrics, ds):
        block = ((metrics.get(ds) or {}).get("sim")) or {}
        mean = (block.get("mean") or {}).get("accuracy_1to1")
        sd = (block.get("sd") or {}).get("accuracy_1to1")
        n = int(block.get("n_metric_reps") or 0)
        return (float(mean) if _finite(mean) else float("nan"),
                float(sd) if _finite(sd) else 0.0, max(n, 1))

    gains, variances = [], []
    used = []
    for ds in sim_ids_eval:
        rm, rs, rn = stats(run_metrics, ds)
        bm, bs, bn = stats(base_metrics, ds)
        if not (_finite(rm) and _finite(bm)):
            continue
        gains.append(rm - bm)
        variances.append(rs * rs / rn + bs * bs / bn)
        used.append(ds)
    if not gains:
        rep.check("sim", "*", "aggregate_accuracy_1to1", "fail",
                  detail="no comparable sim accuracy_1to1 values")
        return
    gain = float(np.mean(gains))
    se = math.sqrt(sum(variances)) / len(gains)
    threshold = max(2.0 * se, MIN_EFFECT_ACCURACY)
    ok = gain > threshold
    bound = "2*SE" if 2.0 * se >= MIN_EFFECT_ACCURACY else "min effect"
    rep.check("sim", "*", "aggregate_accuracy_1to1", "pass" if ok else "fail",
              delta=gain, se=se, threshold=threshold,
              n_datasets=len(used),
              detail=f"mean gain must exceed 2*SE={2 * se:.4f} and "
                     f"min effect={MIN_EFFECT_ACCURACY:.3f} "
                     f"(binding: {bound})")


# ---------------------------------------------------------------------------
# run health / runtime
# ---------------------------------------------------------------------------

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
            b = c.get("baseline_mean", c.get("baseline_rep_mean",
                                             c.get("baseline_value")))
            r = c.get("run_mean", c.get("run_vs_baseline_mean",
                                        c.get("run_value")))
            if b is None and c.get("baseline_value") is not None:
                b = c["baseline_value"]
            tol = c.get("tolerance", c.get("required_min"))
            delta = c.get("delta", c.get("deviation"))
            status = c["status"].upper()
            if c.get("detail") and c["status"] in ("fail", "skip"):
                status += f" ({c['detail']})"
            lines.append(f"| {c['dataset']} | {c['metric']} | {_fmt(b)} | {_fmt(r)} "
                         f"| {_fmt(delta)} | {_fmt(tol)} | {status} |")
        lines.append("")

    if rep.gate_info:
        lines.append("## Gated metrics (tolerance = max(k·SD_pooled, floor))")
        lines.append("")
        lines.append("| kind | metric | pooled SD | floor | tolerance |")
        lines.append("|---|---|---|---|---|")
        for g in rep.gate_info:
            lines.append(f"| {g['kind']} | {g['metric']} | {_fmt(g['pooled_sd'], 5)} "
                         f"| {_fmt(g['floor'])} | {_fmt(g['tolerance'])} |")
        lines.append("")

    if rep.expect == "identical":
        section("Identity (per-replicate assignment sha256)", "identical")
    section("Provenance / content hashes", "provenance")
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


def load_run_selection(run_root: Path) -> set[str]:
    """Dataset ids the run intended to execute (run.py ``_selection.json``).

    Returns an empty set for legacy runs; an empty set means "unknown
    selection" and keeps the strict presence behaviour (every baseline
    dataset must be present in the run).
    """
    p = run_root / "_selection.json"
    if not p.is_file():
        return set()
    try:
        sel = common.read_json(p)
    except (OSError, ValueError):
        return set()
    ids = sel.get("datasets") if isinstance(sel, dict) else None
    return set(ids) if isinstance(ids, list) else set()


def record_gate_info(pooled_sim: dict[str, float], pooled_real: dict[str, float],
                     k: float) -> list[dict]:
    """Tolerance table for the report (same / improved gates)."""
    out = []
    for kind, pooled, floors, default_floor, primary in (
            ("sim", pooled_sim, SIM_FLOORS, DEFAULT_SIM_FLOOR, SIM_PRIMARY_SAME),
            ("real", pooled_real, REAL_FLOORS, DEFAULT_REAL_FLOOR,
             REAL_PRIMARY_SAME)):
        for metric in primary:
            sd = pooled.get(metric, 0.0)
            floor = floors.get(metric, default_floor)
            out.append({"kind": kind, "metric": metric,
                        "pooled_sd": sd, "floor": floor,
                        "tolerance": max(k * sd, floor)})
    return out


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--expect", choices=list(MODES), default="identical",
                    help="identical = bitwise check at 1 thread/1 replicate "
                         "(default, refactor gate); same = unchanged "
                         "algorithm within the noise floor; improved = "
                         "measurably better")
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--baselines-dir", default=None,
                    help="default <repo>/benchmarks/baselines")
    ap.add_argument("--k", type=float, default=K_DEFAULT,
                    help="tolerance multiplier on the pooled baseline SD "
                         "(default 3)")
    ap.add_argument("--admixture-tolerance", type=float,
                    default=DEFAULT_ADMIXTURE_TOLERANCE,
                    help="minimum admixture tolerance floor; the actual "
                         "tolerance is max(k*SD of the baseline audit "
                         "replicates, floor), or the floor alone with "
                         "<2 baseline audit replicates (default "
                         f"{DEFAULT_ADMIXTURE_TOLERANCE}, = 3x the measured "
                         "Baysor replicate audit SD)")
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

    # 'same' needs a real baseline with replicate agreement: with a single
    # successful replicate there is no noise floor to compare against and no
    # threshold may be invented (the old absolute fallback wrongly failed
    # run-vs-baseline ARI ~0.83 against a required 0.90)
    if args.expect == "same":
        for ds_id in sorted(set(run_metrics) & set(base_metrics)):
            bm = base_metrics[ds_id]
            bkind = (bm.get("dataset") or {}).get("kind") or bm.get("kind")
            if bkind != "real":
                continue
            ok_reps = [r for r in bm.get("reps", [])
                       if r.get("status") == "ok"]
            if len(ok_reps) < 2:
                print(f"error: {ds_id}: baseline needs >=3 replicates for "
                      f"real same-mode checks "
                      f"(found {len(ok_reps)} successful replicate(s))",
                      file=sys.stderr)
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
        "threads (run)": sorted({mm.get("threads") for mm in run_metrics.values()
                                 if mm.get("threads") is not None}),
        "threads (baseline)": sorted({mm.get("threads") for mm in base_metrics.values()
                                      if mm.get("threads") is not None}),
    }
    if args.expect == "improved" and run_sha != base_sha:
        rep.meta["binary sha note"] = ("differs from baseline (expected for an "
                                       "algorithm change)")

    # run health
    collect_run_failures(rep, run_metrics)

    # dataset coverage: a baseline dataset missing from the run fails —
    # unless the run recorded an explicit selection that excludes it
    # (deliberate subset run; run.py's _selection.json), which is reported
    # as skipped instead. Legacy runs without a selection file keep the
    # strict behaviour.
    selection = load_run_selection(run_root)
    for ds_id in sorted(set(base_metrics) - set(run_metrics)):
        if selection and ds_id not in selection:
            rep.warn(f"{ds_id}: in the baseline but not in this run's "
                     "dataset selection (skipped)")
            continue
        rep.check("run", ds_id, "presence", "fail",
                  detail="dataset in baseline but missing from run")
    for ds_id in sorted(set(run_metrics) - set(base_metrics)):
        rep.warn(f"{ds_id}: present in run but not in baseline (ignored)")

    # pooled SDs and tolerance table (same / improved)
    pooled_sim: dict[str, float] = {}
    pooled_real: dict[str, float] = {}
    if args.expect in ("same", "improved"):
        pooled_sim = pooled_sd_by_metric(base_metrics, "sim")
        pooled_real = pooled_sd_by_metric(base_metrics, "real")
        rep.gate_info = record_gate_info(pooled_sim, pooled_real, args.k)

    # per-dataset checks
    for ds_id in sorted(set(run_metrics) & set(base_metrics)):
        run_m, base_m = run_metrics[ds_id], base_metrics[ds_id]
        kind = (run_m.get("dataset") or {}).get("kind") or run_m.get("kind")
        check_content_hashes(rep, ds_id, run_m, base_m)
        if args.expect == "identical":
            check_provenance(rep, ds_id, run_m, base_m, "identical")
            run_cells, base_cells = ([], [])
            if kind == "real":
                run_cells, base_cells = load_pair_cells(
                    run_root, root, args.baseline, ds_id, run_m)
            check_identical_dataset(rep, ds_id, run_m, base_m,
                                    run_cells, base_cells)
            continue

        check_provenance(rep, ds_id, run_m, base_m, args.expect)
        if kind == "sim":
            check_sim_dataset(rep, ds_id, run_m, base_m, args.k,
                              pooled_sim, mode=args.expect)
        elif kind == "real":
            run_cells, base_cells = [], []
            if args.expect == "same":
                run_cells, base_cells = load_pair_cells(
                    run_root, root, args.baseline, ds_id, run_m)
                if not base_cells:
                    rep.check("real", ds_id, "baseline_assignments", "fail",
                              detail="missing baseline assignment tables "
                                     "(recreate the baseline)")
            check_real_dataset(rep, ds_id, run_m, base_m, run_cells, base_cells,
                               args.k, args.expect, pooled_real,
                               args.admixture_tolerance)

    if args.expect == "improved":
        check_improved_aggregate(rep, run_metrics, base_metrics, args.k)
        # both gates skipped -> nothing was actually evaluated
        evaluated = any(c["scope"] in ("sim", "real")
                        and c["status"] in ("pass", "fail")
                        for c in rep.checks)
        if not evaluated:
            rep.check("run", "*", "evaluation", "fail",
                      detail="nothing to evaluate: no sim datasets for the "
                             "accuracy gate and no evaluable "
                             "admixture-capable real datasets for the audit "
                             "gate")

    compare_runtime(rep, run_metrics, base_metrics)
    if args.expect in ("same", "improved"):
        rep.finalize_budget(args.k)

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
