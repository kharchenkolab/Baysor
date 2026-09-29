"""End-to-end tests for compare.py on synthetic runs (no Baysor binary).

Covers both expectation modes: self-comparison must pass, a deliberately
degraded run must fail, and the cellAdmix admixture check with tolerances.
"""
import numpy as np
import pandas as pd
import pytest

import baseline
import common
import compare
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
    assert sim_checks and all(c["status"] in ("pass", "skip") for c in sim_checks)


def test_same_degraded_fails(tmp_path):
    root, baselines = _setup(tmp_path, new_degrade=True)
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    failed = [c for c in report["checks"] if c["status"] == "fail"]
    assert any(c["metric"] == "matched_accuracy" for c in failed)


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
    assert real and all(c["status"] == "pass" for c in real)


def test_real_same_degraded_fails(tmp_path):
    root, baselines = _setup(tmp_path, with_real=True, new_degrade=True)
    assert _compare(root, baselines, "same") == 1
    report = common.read_json(
        root / "runs" / "rnew" / "compare_btest_same.json")
    failed = {c["metric"] for c in report["checks"] if c["status"] == "fail"}
    assert "molecule_ari" in failed or "frac_cells_matched" in failed


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
