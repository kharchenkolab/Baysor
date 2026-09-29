"""End-to-end pipeline test with the real Baysor binary.

Skipped when the binary is not on this machine (set ``BAYSOR_BIN`` to use a
different build). Flow: run 3 replicates -> create a temporary baseline ->
compare the run with itself (must pass) -> compare against a deliberately
degraded run (must fail).
"""
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import baseline
import common
import compare
import run as runner
from fixtures import make_sim_dataset, make_real_dataset, make_run, corrupt

DEFAULT_BIN = ("/home/vpetukhov/.bb/thread-storage/thr_cpwic2f6q3/"
               "baysor-bugfixes/build-rel/baysor")
BAYSOR = Path(os.environ.get("BAYSOR_BIN", DEFAULT_BIN))

pytestmark = pytest.mark.skipif(
    not BAYSOR.is_file(), reason=f"Baysor binary not available at {BAYSOR}")


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    root = tmp_path_factory.mktemp("bench_pipeline")
    sim_dir = make_sim_dataset(root / "sim" / "sim_pipe")
    real_dir = make_real_dataset(root / "real" / "real_pipe")
    rc = runner.main([
        "--baysor", str(BAYSOR),
        "--datasets", "quick",
        "--run-id", "p1",
        "--replicates", "3",
        "--threads", "2",
        "--data-root", str(root),
    ])
    assert rc == 0
    baselines = root / "repo_baselines"
    assert baseline.create("p1", "pipe", root, baselines) == 0
    return {"root": root, "baselines": baselines,
            "sim_dir": sim_dir, "real_dir": real_dir}


def test_run_outputs(pipeline):
    root = pipeline["root"]
    for ds in ("sim_pipe", "real_pipe"):
        m = common.read_json(root / "runs" / "p1" / ds / "metrics.json")
        assert m["replicates"] == 3
        assert not m["failures"]
        assert m["binary"]["sha256"]
        for k in range(3):
            rep = root / "runs" / "p1" / ds / f"rep{k}"
            assert (rep / "assignment.parquet").is_file()
            assert (rep / "run.json").is_file()
            assert (rep / "seg" / "molecules.parquet").is_file()
            rec = common.read_json(rep / "run.json")
            assert rec["status"] == "ok"
            assert rec["wall_s"] > 0
            assert rec["peak_rss_kb"] and rec["peak_rss_kb"] > 0
            assert "-s" in rec["command"]
    m = common.read_json(root / "runs" / "p1" / "sim_pipe" / "metrics.json")
    assert m["sim"]["mean"]["matched_accuracy"] > 0.85
    assert m["sim"]["n_metric_reps"] == 3
    m = common.read_json(root / "runs" / "p1" / "real_pipe" / "metrics.json")
    assert m["real"]["rep_agreement"]["mean"]["molecule_ari"] > 0.5


def test_assignment_alignment(pipeline):
    root = pipeline["root"]
    src = pd.read_parquet(pipeline["sim_dir"] / "molecules.parquet")
    a = pd.read_parquet(root / "runs" / "p1" / "sim_pipe" /
                        "rep0" / "assignment.parquet")
    assert len(a) == len(src)
    assert a["mol_index"].tolist() == list(range(len(src)))
    assert (a["cell"].to_numpy() >= 0).all()


def test_self_comparison_passes(pipeline):
    root, baselines = pipeline["root"], pipeline["baselines"]
    rc = compare.main(["--run-id", "p1", "--baseline", "pipe",
                       "--expect", "same",
                       "--data-root", str(root),
                       "--baselines-dir", str(baselines)])
    assert rc == 0
    report = common.read_json(
        root / "runs" / "p1" / "compare_pipe_same.json")
    assert report["summary"]["pass_overall"] is True


def test_degraded_run_fails(pipeline, tmp_path):
    """A run with 20% of molecules randomly reassigned must not compare 'same'."""
    root, baselines = pipeline["root"], pipeline["baselines"]
    # build a degraded candidate from the baseline's own assignments
    for ds_id, kind in (("sim_pipe", "sim"), ("real_pipe", "real")):
        cells = common.assignment_cells(
            root / "runs" / "p1" / ds_id / "rep0" / "assignment.parquet")
        degraded = [corrupt(cells, 0.2, seed=90 + k) for k in range(3)]
        make_run(root, "p_degraded", root / kind / ds_id, degraded,
                 binary={"sha256": common.sha256_file(BAYSOR),
                         "path": str(BAYSOR)})
    rc = compare.main(["--run-id", "p_degraded", "--baseline", "pipe",
                       "--expect", "same",
                       "--data-root", str(root),
                       "--baselines-dir", str(baselines)])
    assert rc == 1
    report = common.read_json(
        root / "runs" / "p_degraded" / "compare_pipe_same.json")
    failed = {c["metric"] for c in report["checks"] if c["status"] == "fail"}
    assert "matched_accuracy" in failed or "molecule_ari" in failed
