"""End-to-end tests for compare.py on synthetic runs (no Baysor binary).

Covers all three expectation modes: identical (sha256 gate), same
(primary-metric gates, provenance/content-hash gates, pooled tolerances)
and improved (gain above noise + regression gates + admixture conditions),
plus self-comparison passes, degraded-run failures and reports.
"""
import numpy as np
import pandas as pd
import pytest

import baseline
import common
import compare
import metrics as m
from fixtures import (make_sim_dataset, make_real_dataset, make_run, corrupt)


def _setup(tmp_path, *, n_reps=3, sim_degrade=False, run_id_base="rbase",
           run_id_new="rnew", with_real=False, base_rates=None,
           run_rates=None, new_degrade=False, new_wall=1.0):
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"

    sim_dir = make_sim_dataset(root / "sim" / "sim_a")
    truth = pd.read_parquet(sim_dir / "molecules.parquet")["cell"].to_numpy(np.int64)
    base_reps = [corrupt(truth, 0.2, seed=20 + k) if sim_degrade
                 else truth.copy() for k in range(n_reps)]
    new_reps = [corrupt(truth, 0.2, seed=40 + k) if new_degrade
                else truth.copy() for k in range(n_reps)]
    make_run(root, run_id_base, sim_dir, base_reps)
    make_run(root, run_id_new, sim_dir, new_reps, wall_s=new_wall)

    if with_real:
        real_dir = make_real_dataset(root / "real" / "real_a")
        cells = pd.read_parquet(real_dir / "molecules.parquet")[
            "cell"].to_numpy(np.int64)
        make_run(root, run_id_base, real_dir, [cells] * n_reps,
                 celladmix_rates=base_rates)
        new_cells = ([corrupt(cells, 0.2, seed=60 + k) if new_degrade
                      else cells.copy() for k in range(n_reps)])
        make_run(root, run_id_new, real_dir, new_cells,
                 celladmix_rates=run_rates, wall_s=new_wall)

    assert baseline.create(run_id_base, "btest", root, baselines,
                           force=True) == 0
    return root, baselines


def _compare(root, baselines, expect, run_id_new="rnew", **kw):
    argv = ["--run-id", run_id_new, "--baseline", "btest",
            "--expect", expect, "--data-root", str(root),
            "--baselines-dir", str(baselines)]
    for k, v in kw.items():
        argv += [f"--{k.replace('_', '-')}", str(v)]
    return compare.main(argv)


def test_same_identical_passes(tmp_path):
    root, baselines = _setup(tmp_path)
    assert _compare(root, baselines, "same") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    assert report["summary"]["pass_overall"] is True
    sim_checks = [c for c in report["checks"] if c["scope"] == "sim"]
    assert sim_checks and all(c["status"] in ("pass", "skip", "info")
                              for c in sim_checks)
    gated = [c for c in sim_checks if c["status"] == "pass"]
    assert {c["metric"] for c in gated} == set(compare.SIM_PRIMARY_SAME)
    # a false-alarm budget is reported for the gated checks
    assert report["meta"]["false-alarm budget"]
    assert report["summary"]["fail"] == 0


def test_same_degraded_fails(tmp_path):
    root, baselines = _setup(tmp_path, new_degrade=True)
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    failed = [c for c in report["checks"] if c["status"] == "fail"]
    assert any(c["metric"] == "accuracy_1to1" for c in failed)
    assert any(c["metric"] == "ari_assigned" for c in failed)


def test_improved_passes(tmp_path):
    # baseline degraded, candidate better
    root, baselines = _setup(tmp_path, sim_degrade=True)
    assert _compare(root, baselines, "improved") == 0


def test_improved_fails_when_worse(tmp_path):
    root, baselines = _setup(tmp_path, new_degrade=True)
    assert _compare(root, baselines, "improved") == 1


def test_improved_requires_increase(tmp_path):
    # identical runs: mean accuracy does not increase
    root, baselines = _setup(tmp_path)
    assert _compare(root, baselines, "improved") == 1


def test_real_same_identical_passes(tmp_path):
    root, baselines = _setup(tmp_path, with_real=True)
    assert _compare(root, baselines, "same") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    real = [c for c in report["checks"] if c["scope"] == "real"]
    assert real and all(c["status"] in ("pass", "info") for c in real)
    gated = [c for c in real if c["status"] == "pass"]
    assert {c["metric"] for c in gated} == set(compare.REAL_PRIMARY_SAME)


def test_real_same_degraded_fails(tmp_path):
    root, baselines = _setup(tmp_path, with_real=True, new_degrade=True)
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    failed = {c["metric"] for c in report["checks"] if c["status"] == "fail"}
    assert "ari_assigned" in failed or "frac_cells_matched" in failed


def test_admixture_within_tolerance(tmp_path):
    root, baselines = _setup(tmp_path, sim_degrade=True, with_real=True,
                             base_rates=[0.05, 0.05, 0.05],
                             run_rates=[0.06, 0.06, 0.06])
    assert _compare(root, baselines, "improved",
                    admixture_tolerance=0.01) == 0


def test_admixture_above_tolerance_fails(tmp_path):
    root, baselines = _setup(tmp_path, sim_degrade=True, with_real=True,
                             base_rates=[0.05, 0.05, 0.05],
                             run_rates=[0.08, 0.08, 0.08])
    assert _compare(root, baselines, "improved",
                    admixture_tolerance=0.01) == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    adm = [c for c in report["checks"]
           if c["metric"] == "total_admixture_rate"]
    assert adm and adm[0]["status"] == "fail"


def test_admixture_unavailable_skips(tmp_path):
    root, baselines = _setup(tmp_path, sim_degrade=True, with_real=True)
    # no celladmix rates anywhere -> graceful skip, sim improvement decides
    assert _compare(root, baselines, "improved") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    adm = [c for c in report["checks"]
           if c["metric"] == "total_admixture_rate"]
    assert adm and adm[0]["status"] == "skip"
    # never reported as 0
    assert adm[0].get("run_mean") is None
    assert "unavailable" in (adm[0].get("detail") or "")


def test_runtime_slowdown_warns_but_passes(tmp_path):
    root, baselines = _setup(tmp_path, new_wall=2.0)
    assert _compare(root, baselines, "same") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    assert any("runtime" in w for w in report["warnings"])
    # identical RSS -> no RSS-growth warning
    assert not any("peak RSS" in w for w in report["warnings"])


def test_failed_rep_fails_comparison(tmp_path):
    root, baselines = _setup(tmp_path)
    # corrupt one replicate's status in the new run
    mpath = root / "runs" / "rnew" / "sim_a" / "metrics.json"
    m = common.read_json(mpath)
    m["reps"][1]["status"] = "timeout"
    m["reps"][1]["exit_code"] = 124
    m["failures"] = [m["reps"][1]]
    common.write_json(mpath, m)
    assert _compare(root, baselines, "same") == 1


def test_missing_dataset_fails(tmp_path):
    root, baselines = _setup(tmp_path)
    # add a second dataset to the baseline only
    sim_b = make_sim_dataset(root / "sim" / "sim_b")
    truth_b = pd.read_parquet(sim_b / "molecules.parquet")["cell"].to_numpy(np.int64)
    make_run(root, "rbase", sim_b, [truth_b] * 3)
    assert baseline.create("rbase", "btest", root, baselines, force=True) == 0
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    assert any(c["metric"] == "presence" and c["status"] == "fail"
               for c in report["checks"])


def test_missing_run_or_baseline_exit_2(tmp_path):
    root, baselines = _setup(tmp_path)
    assert _compare(root, baselines, "same", run_id_new="nope") == 2
    argv = ["--run-id", "rnew", "--baseline", "nope", "--expect", "same",
            "--data-root", str(root), "--baselines-dir", str(baselines)]
    assert compare.main(argv) == 2


def test_markdown_report_written(tmp_path):
    root, baselines = _setup(tmp_path, with_real=True)
    assert _compare(root, baselines, "same") == 0
    md = (root / "runs" / "rnew" / "compare_btest_same.md").read_text()
    assert "# Benchmark comparison" in md
    assert "## Sim datasets" in md
    assert "## Real datasets" in md
    assert "## Runtime & memory" in md
    assert "PASS" in md


# ---------------------------------------------------------------------------
# --expect identical (default refactor gate)
# ---------------------------------------------------------------------------

def _setup_identical(tmp_path, *, new_corrupt=False, threads=1, reps=1,
                     new_threads=None, new_binary=None):
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"
    sim_dir = make_sim_dataset(root / "sim" / "sim_a")
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    base_reps = [truth.copy() for _ in range(reps)]
    new_reps = ([corrupt(truth, 0.05, seed=99) for _ in range(reps)]
                if new_corrupt else [truth.copy() for _ in range(reps)])
    make_run(root, "rbase", sim_dir, base_reps, threads=threads)
    make_run(root, "rnew", sim_dir, new_reps,
             threads=new_threads if new_threads is not None else threads,
             binary=new_binary)
    assert baseline.create("rbase", "btest", root, baselines,
                           force=True, allow_incomplete=True) == 0
    return root, baselines


def _compare_raw(root, baselines, run_id="rnew", expect=None):
    argv = ["--run-id", run_id, "--baseline", "btest",
            "--data-root", str(root), "--baselines-dir", str(baselines)]
    if expect is not None:
        argv += ["--expect", expect]
    return compare.main(argv)


def test_identical_is_the_default_and_passes(tmp_path):
    root, baselines = _setup_identical(tmp_path)
    # no --expect given -> identical
    assert _compare_raw(root, baselines) == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_identical.json")
    assert report["expect"] == "identical"
    assert report["summary"]["pass_overall"] is True
    sha_rows = [c for c in report["checks"]
                if c["scope"] == "identical"]
    assert sha_rows and all(c["status"] == "pass" for c in sha_rows)
    assert any("threads_1" in c["metric"] for c in report["checks"])


def test_identical_sha_mismatch_fails_with_metric_deltas(tmp_path):
    root, baselines = _setup_identical(tmp_path, new_corrupt=True)
    assert _compare_raw(root, baselines, expect="identical") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_identical.json")
    failed = [c for c in report["checks"] if c["status"] == "fail"]
    assert any("assignment_sha256" in c["metric"] for c in failed)
    # metric deltas are shown alongside the mismatch
    deltas = [c for c in report["checks"]
              if c["scope"] == "identical" and c["status"] == "info"
              and c.get("delta") is not None]
    assert any(c["metric"] == "accuracy_1to1" for c in deltas)
    assert any(abs(c["delta"]) > 0 for c in deltas)


def test_identical_requires_one_thread(tmp_path):
    root, baselines = _setup_identical(tmp_path, threads=6)
    assert _compare_raw(root, baselines, expect="identical") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_identical.json")
    fails = {(c["dataset"], c["metric"]) for c in report["checks"]
             if c["status"] == "fail"}
    assert ("sim_a", "run_threads_1") in fails
    assert ("sim_a", "baseline_threads_1") in fails


def test_identical_requires_exactly_one_replicate(tmp_path):
    root, baselines = _setup_identical(tmp_path, reps=3)
    assert _compare_raw(root, baselines, expect="identical") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_identical.json")
    fails = {c["metric"] for c in report["checks"] if c["status"] == "fail"}
    assert "run_replicates_1" in fails and "baseline_replicates_1" in fails


# ---------------------------------------------------------------------------
# --expect same: provenance and content gates
# ---------------------------------------------------------------------------

def test_same_thread_mismatch_fails(tmp_path):
    root, baselines = _setup(tmp_path)
    # rebuild the candidate with a different thread count
    sim_dir = root / "sim" / "sim_a"
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rnew", sim_dir, [truth.copy()] * 3, threads=3)
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    row = [c for c in report["checks"] if c["metric"] == "threads"]
    assert row and row[0]["status"] == "fail"


def test_same_binary_mismatch_fails(tmp_path):
    root, baselines = _setup(tmp_path)
    sim_dir = root / "sim" / "sim_a"
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rnew", sim_dir, [truth.copy()] * 3,
             binary={"sha256": "1" * 64})
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    row = [c for c in report["checks"] if c["metric"] == "binary_sha256"]
    assert row and row[0]["status"] == "fail"


def test_same_content_hash_mismatch_fails(tmp_path):
    root, baselines = _setup(tmp_path)
    sim_dir = root / "sim" / "sim_a"
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rnew", sim_dir, [truth.copy()] * 3,
             inputs={"molecules_sha256": "0" * 64,
                     "meta_sha256": "0" * 64})
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    row = [c for c in report["checks"]
           if c["metric"] == "molecules_sha256"]
    assert row and row[0]["status"] == "fail"


def test_same_content_hash_missing_skips(tmp_path):
    root, baselines = _setup(tmp_path)
    sim_dir = root / "sim" / "sim_a"
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rnew", sim_dir, [truth.copy()] * 3, inputs=None)
    assert _compare(root, baselines, "same") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    row = [c for c in report["checks"]
           if c["metric"] == "molecules_sha256"]
    assert row and row[0]["status"] == "skip"


# ---------------------------------------------------------------------------
# --expect improved: gain above noise, regression and admixture gates
# ---------------------------------------------------------------------------

def test_improved_min_effect(tmp_path):
    """A gain below the minimum effect (noise) must not pass improved."""
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"
    sim_dir = make_sim_dataset(root / "sim" / "sim_a")
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    # baseline: one molecule mislabelled in every replicate (~0.003 gain)
    base = truth.copy()
    base[0] = base[1] if base[1] != base[0] else base[2]
    make_run(root, "rbase", sim_dir, [base.copy()] * 3)
    make_run(root, "rnew", sim_dir, [truth.copy()] * 3)
    assert baseline.create("rbase", "btest", root, baselines,
                           force=True) == 0
    assert _compare(root, baselines, "improved") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    agg = [c for c in report["checks"]
           if c["metric"] == "aggregate_accuracy_1to1"]
    assert agg and agg[0]["status"] == "fail"
    assert agg[0]["threshold"] >= compare.MIN_EFFECT_ACCURACY
    assert 0 < agg[0]["delta"] < compare.MIN_EFFECT_ACCURACY


def test_improved_aggregate_two_standard_errors():
    """Gain above the min effect but below 2 SE still fails."""
    def block(mean, sd):
        return {"dataset": {"kind": "sim"},
                "sim": {"mean": {"accuracy_1to1": mean},
                        "sd": {"accuracy_1to1": sd},
                        "n_metric_reps": 3}}
    base = {"d": block(0.90, 0.06)}
    run = {"d": block(0.95, 0.01)}
    rep = compare.Report("r", "b", "improved")
    compare.check_improved_aggregate(rep, run, base, 3.0)
    row = rep.checks[0]
    se = np.sqrt(0.06 ** 2 / 3 + 0.01 ** 2 / 3)
    assert row["se"] == pytest.approx(se)
    assert row["delta"] == pytest.approx(0.05)
    assert row["status"] == "fail"          # 0.05 <= 2*SE ~ 0.070
    # a larger gain passes
    run["d"]["sim"]["mean"]["accuracy_1to1"] = 0.99
    rep2 = compare.Report("r", "b", "improved")
    compare.check_improved_aggregate(rep2, run, base, 3.0)
    assert rep2.checks[0]["status"] == "pass"


def _adm_metrics(status, rate, n_cells):
    return {"dataset": {"kind": "real"},
            "reps": [{"status": "ok", "n_cells": n_cells}],
            "real": {"celladmix": {
                "status": status,
                "mean_total": rate,
                "per_rep_total": [rate] * 3 if rate is not None else []}}}


def test_admixture_gate_conditions():
    # audit not ok on the baseline -> unavailable
    rep = compare.Report("r", "b", "improved")
    compare.check_real_admixture(rep, "d",
                                 _adm_metrics("ok", 0.04, 2500),
                                 _adm_metrics("failed", None, 2500), 3.0, 0.0)
    assert rep.checks[0]["status"] == "skip"
    assert "unavailable" in rep.checks[0]["detail"]
    # audit ok but the crop is too small -> unavailable (never 0)
    rep = compare.Report("r", "b", "improved")
    compare.check_real_admixture(rep, "d",
                                 _adm_metrics("ok", 0.04, 615),
                                 _adm_metrics("ok", 0.04, 615), 3.0, 0.0)
    assert rep.checks[0]["status"] == "skip"
    assert "615" in rep.checks[0]["detail"]
    assert "run_mean" not in rep.checks[0]
    # ok on both sides, large crop: tolerance = k * SD(baseline reps)
    base = dict(_adm_metrics("ok", 0.05, 2500))
    base["real"]["celladmix"]["per_rep_total"] = [0.05, 0.06, 0.05]
    rep = compare.Report("r", "b", "improved")
    compare.check_real_admixture(rep, "d", _adm_metrics("ok", 0.065, 2500),
                                 base, 3.0, 0.0)
    row = rep.checks[0]
    sd = float(np.std([0.05, 0.06, 0.05], ddof=1))
    assert row["status"] == "pass"
    assert row["tolerance"] == pytest.approx(3 * sd)
    assert row["allowed_max"] == pytest.approx(0.05 + 3 * sd)
    # zero baseline spread -> tolerance 0 -> any rise fails
    rep = compare.Report("r", "b", "improved")
    compare.check_real_admixture(rep, "d", _adm_metrics("ok", 0.06, 2500),
                                 _adm_metrics("ok", 0.05, 2500), 3.0, 0.0)
    assert rep.checks[0]["status"] == "fail"


def _setup_big_real(tmp_path, *, base_rates, run_rates):
    """Baseline with degraded sim accuracy (so improved can pass) plus a
    >= 2000-cell real dataset carrying cellAdmix rates."""
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"
    sim_dir = make_sim_dataset(root / "sim" / "sim_a")
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rbase", sim_dir, [corrupt(truth, 0.2, seed=20 + k)
                                      for k in range(3)])
    make_run(root, "rnew", sim_dir, [truth.copy()] * 3)
    real_dir = make_real_dataset(root / "real" / "real_big",
                                 n_cells=2000, per_cell=1, noise=50)
    cells = pd.read_parquet(real_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rbase", real_dir, [cells] * 3, celladmix_rates=base_rates)
    make_run(root, "rnew", real_dir, [cells] * 3, celladmix_rates=run_rates)
    assert baseline.create("rbase", "btest", root, baselines,
                           force=True) == 0
    return root, baselines


def test_admixture_within_tolerance_passes(tmp_path):
    root, baselines = _setup_big_real(
        tmp_path, base_rates=[0.05, 0.06, 0.05], run_rates=[0.06] * 3)
    assert _compare(root, baselines, "improved") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    adm = [c for c in report["checks"]
           if c["metric"] == "total_admixture_rate"]
    assert adm and adm[0]["status"] == "pass"
    assert adm[0]["tolerance"] > 0


def test_admixture_above_tolerance_fails(tmp_path):
    root, baselines = _setup_big_real(
        tmp_path, base_rates=[0.05] * 3, run_rates=[0.06] * 3)
    assert _compare(root, baselines, "improved") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    adm = [c for c in report["checks"]
           if c["metric"] == "total_admixture_rate"]
    assert adm and adm[0]["status"] == "fail"


# ---------------------------------------------------------------------------
# tolerance helpers
# ---------------------------------------------------------------------------

def test_pooled_sd_by_metric():
    base = {
        "a": {"dataset": {"kind": "sim"},
              "sim": {"sd": {"accuracy_1to1": 0.03, "recovery_rate": 0.1},
                      "n_metric_reps": 3}},
        "b": {"dataset": {"kind": "sim"},
              "sim": {"sd": {"accuracy_1to1": 0.06},
                      "n_metric_reps": 3}},
        "c": {"dataset": {"kind": "real"},
              "real": {"rep_agreement": {"sd": {"ari_assigned": 0.01},
                                         "per_pair": [{}, {}, {}]}}},
        # fewer than 2 replicates -> excluded
        "d": {"dataset": {"kind": "sim"},
              "sim": {"sd": {"accuracy_1to1": 9.9}, "n_metric_reps": 1}},
    }
    pooled = compare.pooled_sd_by_metric(base, "sim")
    # sqrt(((3-1)*0.03^2 + (3-1)*0.06^2) / 4)
    assert pooled["accuracy_1to1"] == pytest.approx(
        np.sqrt((2 * 0.03 ** 2 + 2 * 0.06 ** 2) / 4))
    # metric present in only one dataset -> its SD passes through
    assert pooled["recovery_rate"] == pytest.approx(0.1)
    real = compare.pooled_sd_by_metric(base, "real")
    assert real["ari_assigned"] == pytest.approx(0.01)


def test_false_alarm_alpha():
    # 3-sigma two-sided normal tail ~0.0027
    assert compare._two_sided_alpha(3.0, 1.0) == pytest.approx(0.0027, rel=1e-3)
    assert compare._one_sided_alpha(3.0, 1.0) == pytest.approx(0.00135, rel=1e-3)
    assert compare._two_sided_alpha(0.02, 0.0) == 0.0


# ---------------------------------------------------------------------------
# round-2 integration fixes (single-replicate baselines, skip semantics)
# ---------------------------------------------------------------------------

def _setup_real_single_rep(tmp_path, *, base_rates=None, run_rates=None,
                           n_cells=2000):
    """Real-only baseline + candidate, each with a single replicate."""
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"
    real_dir = make_real_dataset(root / "real" / "real_a",
                                 n_cells=n_cells, per_cell=1, noise=50)
    cells = pd.read_parquet(real_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rbase", real_dir, [cells], celladmix_rates=base_rates)
    make_run(root, "rnew", real_dir, [cells], celladmix_rates=run_rates)
    assert baseline.create("rbase", "btest", root, baselines,
                           force=True, allow_incomplete=True) == 0
    return root, baselines


def test_same_real_baseline_without_replicates_exits_2(tmp_path, capsys):
    """A 1-replicate real baseline has no noise floor: same-mode must error
    instead of inventing an absolute threshold."""
    root, baselines = _setup_real_single_rep(tmp_path)
    rc = _compare(root, baselines, "same")
    assert rc == 2
    err = capsys.readouterr().err
    assert "baseline needs >=3 replicates for real same-mode checks" in err
    assert "real_a" in err


def test_same_two_replicate_real_baseline_runs(tmp_path):
    """2 successful replicates give one agreement pair -> no exit 2."""
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"
    real_dir = make_real_dataset(root / "real" / "real_a")
    cells = pd.read_parquet(real_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rbase", real_dir, [cells, cells.copy()])
    make_run(root, "rnew", real_dir, [cells, cells.copy()])
    assert baseline.create("rbase", "btest", root, baselines,
                           force=True, allow_incomplete=True) == 0
    assert _compare(root, baselines, "same") == 0


def test_same_single_rep_sim_baseline_still_runs(tmp_path):
    """Sim same-mode checks work from pooled SDs/floors without replicates."""
    root, baselines = _setup_identical(tmp_path, reps=1)
    # _setup_identical builds a 1-replicate baseline with --allow-incomplete
    assert _compare(root, baselines, "same") == 0


def test_admixture_single_baseline_replicate_uses_floor(tmp_path):
    """<2 baseline audit replicates: tolerance = --admixture-tolerance floor
    (default 0.0025) with a warning, not SD=0."""
    root, baselines = _setup_real_single_rep(
        tmp_path, base_rates=[0.05], run_rates=[0.052])
    assert _compare(root, baselines, "improved",
                    admixture_tolerance=compare.DEFAULT_ADMIXTURE_TOLERANCE) == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    adm = [c for c in report["checks"]
           if c["metric"] == "total_admixture_rate"]
    assert adm and adm[0]["status"] == "pass"
    assert adm[0]["tolerance"] == pytest.approx(0.0025)
    assert adm[0]["allowed_max"] == pytest.approx(0.05 + 0.0025)
    assert adm[0]["baseline_sd"] is None
    assert any("floor" in w for w in report["warnings"])
    # a rise beyond the floor still fails
    root, baselines = _setup_real_single_rep(
        tmp_path / "b", base_rates=[0.05], run_rates=[0.054])
    assert _compare(root, baselines, "improved",
                    admixture_tolerance=compare.DEFAULT_ADMIXTURE_TOLERANCE) == 1


def test_admixture_default_floor_is_calibrated():
    assert compare.DEFAULT_ADMIXTURE_TOLERANCE == 0.0025


def test_admixture_not_capable_unavailable():
    rep = compare.Report("r", "b", "improved")
    base = _adm_metrics("ok", 0.05, 2500)
    base["real"]["celladmix"]["admixture_capable"] = False
    compare.check_real_admixture(rep, "d", _adm_metrics("ok", 0.05, 2500),
                                 base, 3.0, 0.0025)
    assert rep.checks[0]["status"] == "skip"
    assert "not admixture_capable" in rep.checks[0]["detail"]


def test_improved_real_only_baseline_skips_sim_aggregate(tmp_path):
    """Real-only baseline: the sim aggregate gate is skip, the admixture
    gate decides (issue: it used to fail with 'no sim datasets')."""
    root, baselines = _setup_real_single_rep(
        tmp_path, base_rates=[0.05], run_rates=[0.05])
    assert _compare(root, baselines, "improved") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    agg = [c for c in report["checks"]
           if c["metric"] == "aggregate_accuracy_1to1"]
    assert agg and agg[0]["status"] == "skip"
    assert "no sim datasets" in agg[0]["detail"]
    adm = [c for c in report["checks"]
           if c["metric"] == "total_admixture_rate"]
    assert adm and adm[0]["status"] == "pass"


def test_improved_nothing_to_evaluate_fails(tmp_path):
    """Real-only baseline without any audit: nothing is evaluated ->
    clear failure instead of a green run."""
    root, baselines = _setup_real_single_rep(tmp_path)  # no audit rates
    assert _compare(root, baselines, "improved") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    agg = [c for c in report["checks"]
           if c["metric"] == "aggregate_accuracy_1to1"]
    assert agg and agg[0]["status"] == "skip"
    adm = [c for c in report["checks"]
           if c["metric"] == "total_admixture_rate"]
    assert adm and adm[0]["status"] == "skip"
    nothing = [c for c in report["checks"] if c["metric"] == "evaluation"]
    assert nothing and nothing[0]["status"] == "fail"
    assert "nothing to evaluate" in nothing[0]["detail"]


def test_improved_run_without_sim_skips_aggregate(tmp_path):
    """Baseline has sim datasets, run does not -> aggregate is skip."""
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"
    sim_dir = make_sim_dataset(root / "sim" / "sim_a")
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rbase", sim_dir, [corrupt(truth, 0.2, seed=5)] * 3)
    real_dir = make_real_dataset(root / "real" / "real_a")
    cells = pd.read_parquet(real_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    make_run(root, "rbase", real_dir, [cells] * 3)
    make_run(root, "rnew", real_dir, [cells] * 3)
    assert baseline.create("rbase", "btest", root, baselines,
                           force=True) == 0
    assert _compare(root, baselines, "improved") == 1  # presence fail: sim_a
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_improved.json")
    agg = [c for c in report["checks"]
           if c["metric"] == "aggregate_accuracy_1to1"]
    assert agg and agg[0]["status"] == "skip"
    assert "run contains no sim datasets" in agg[0]["detail"]


def test_same_self_comparison_with_noisy_replicates_passes(tmp_path):
    """Self-comparison (run = baseline's source) must pass even when the
    baseline's own replicates disagree.

    The run x baseline pair set contains identity pairs, so its mean sits
    systematically above the replicate-agreement centre by ~(1-centre)/3;
    the real agreement gates are therefore one-sided (a run must not agree
    *less* than the baseline's replicates agree with each other). The three
    replicates below disagree pairwise by a similar amount (low centre, low
    spread), which is exactly the situation the old two-sided gate failed
    and the root cause of the flaky pipeline self-comparison at >1 thread.
    """
    root = tmp_path / "data"
    baselines = tmp_path / "baselines"
    real_dir = make_real_dataset(root / "real" / "real_a")
    cells = pd.read_parquet(real_dir / "molecules.parquet")[
        "cell"].to_numpy(np.int64)
    # each replicate perturbs a disjoint ~10% of molecules: all three pairs
    # disagree by nearly the same amount -> low centre, tiny spread
    rng = np.random.default_rng(7)
    order = rng.permutation(len(cells))
    chunk = len(cells) // 10
    chunks = [order[k * chunk:(k + 1) * chunk] for k in range(3)]
    others = np.unique(cells[cells > 0])

    def variant(k):
        out = cells.copy()
        for i in chunks[k]:
            cand = others[others != out[i]]
            out[i] = int(rng.choice(cand if len(cand) else others))
        return out

    noisy = [variant(0), variant(1), variant(2)]
    make_run(root, "rbase", real_dir, noisy)
    assert baseline.create("rbase", "btest", root, baselines,
                           force=True) == 0
    m = common.read_json(root / "runs" / "rbase" / "real_a" / "metrics.json")
    center = m["real"]["rep_agreement"]["mean"]["ari_assigned"]
    tol = max(3.0 * m["real"]["rep_agreement"]["sd"]["ari_assigned"],
              compare.REAL_FLOORS["ari_assigned"])
    # the old two-sided gate would have failed: mean - centre = (1-c)/3
    assert (1.0 - center) / 3.0 > tol
    assert center < 0.95
    # the self-comparison still passes
    assert _compare(root, baselines, "same", run_id_new="rbase") == 0
    report = common.read_json(
        root / "runs" / "rbase" / "compare_btest_same.json")
    row = [c for c in report["checks"] if c["metric"] == "ari_assigned"]
    assert row and row[0]["status"] == "pass"
    assert row[0]["run_vs_baseline_mean"] > row[0]["baseline_rep_mean"]


# ---------------------------------------------------------------------------
# run selection (_selection.json): deliberate subset runs
# ---------------------------------------------------------------------------

def _write_selection(run_root, ids, spec="custom"):
    import json
    with open(run_root / "_selection.json", "w") as fh:
        json.dump({"datasets": sorted(ids),
                   "invocations": [{"spec": spec, "kind": None}]}, fh)


def test_subset_run_skips_unselected_baseline_datasets(tmp_path):
    root, baselines = _setup(tmp_path)
    sim_b = make_sim_dataset(root / "sim" / "sim_b")
    truth_b = pd.read_parquet(sim_b / "molecules.parquet")["cell"].to_numpy(np.int64)
    make_run(root, "rbase", sim_b, [truth_b] * 3)
    assert baseline.create("rbase", "btest", root, baselines, force=True) == 0
    # rnew covers only sim_a and says so: sim_b must be skipped, not failed
    _write_selection(root / "runs" / "rnew", ["sim_a"])
    assert _compare(root, baselines, "same") == 0
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    assert not any(c["metric"] == "presence" and c["status"] == "fail"
                   for c in report["checks"])
    assert any("not in this run's dataset selection" in w
               for w in report.get("warnings", []))


def test_selected_dataset_missing_from_run_still_fails(tmp_path):
    root, baselines = _setup(tmp_path)
    sim_b = make_sim_dataset(root / "sim" / "sim_b")
    truth_b = pd.read_parquet(sim_b / "molecules.parquet")["cell"].to_numpy(np.int64)
    make_run(root, "rbase", sim_b, [truth_b] * 3)
    assert baseline.create("rbase", "btest", root, baselines, force=True) == 0
    # the selection claims sim_b, but its metrics.json is absent -> fail
    _write_selection(root / "runs" / "rnew", ["sim_a", "sim_b"])
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    assert any(c["metric"] == "presence" and c["status"] == "fail"
               for c in report["checks"])


def test_record_selection_merges_invocations(tmp_path):
    import argparse
    import run as runner
    run_root = tmp_path / "runs" / "rx"
    run_root.mkdir(parents=True)

    class DS:
        def __init__(self, i):
            self.id = i

    args = argparse.Namespace(datasets="quick", kind=None, replicates=3,
                              threads=6, scale_factor=1.0)
    runner.record_selection(run_root, [DS("a"), DS("b")], args)
    args2 = argparse.Namespace(datasets="full", kind="real", replicates=1,
                               threads=1, scale_factor=0.9)
    runner.record_selection(run_root, [DS("b"), DS("c")], args2)
    sel = common.read_json(run_root / "_selection.json")
    assert sel["datasets"] == ["a", "b", "c"]
    assert [i["spec"] for i in sel["invocations"]] == ["quick", "full"]
    assert sel["invocations"][1]["scale_factor"] == 0.9
    assert compare.load_run_selection(run_root) == {"a", "b", "c"}
