"""Tests for recompute_metrics.py (recompute metrics.json from assignments)."""
import pandas as pd
import pytest

import common
import recompute_metrics as recompute
from fixtures import make_sim_dataset, make_real_dataset, make_run


def test_recompute_sim_adds_inputs(tmp_path):
    root = tmp_path / "data"
    sim_dir = make_sim_dataset(root / "sim" / "sim_a")
    truth = pd.read_parquet(sim_dir / "molecules.parquet")[
        "cell"].to_numpy("int64")
    make_run(root, "r1", sim_dir, [truth.copy()] * 3, inputs=None)
    mpath = root / "runs" / "r1" / "sim_a" / "metrics.json"
    assert "inputs" not in common.read_json(mpath)

    assert recompute.main(["--run", "r1", "--data-root", str(root)]) == 0
    m = common.read_json(mpath)
    assert m["inputs"]["molecules_sha256"] == common.sha256_file(
        sim_dir / "molecules.parquet")
    assert m["inputs"]["meta_sha256"] == common.sha256_file(
        sim_dir / "meta.json")
    assert m["metrics_recomputed"]["script"] == "recompute_metrics.py"
    assert m["sim"]["mean"]["accuracy_1to1"] == 1.0
    # provenance preserved
    assert m["replicates"] == 3


def test_recompute_replaces_metric_blocks_keeps_audit(tmp_path):
    root = tmp_path / "data"
    real_dir = make_real_dataset(root / "real" / "real_a")
    cells = pd.read_parquet(real_dir / "molecules.parquet")[
        "cell"].to_numpy("int64")
    make_run(root, "r2", real_dir, [cells] * 3,
             celladmix_rates=[0.1, 0.1, 0.1])
    mpath = root / "runs" / "r2" / "real_a" / "metrics.json"
    m = common.read_json(mpath)
    # simulate stale metrics: wipe the agreement block, keep the audit
    m["real"]["rep_agreement"] = {}
    common.write_json(mpath, m)

    assert recompute.main(["--run", "r2", "--data-root", str(root)]) == 0
    m = common.read_json(mpath)
    ra = m["real"]["rep_agreement"]["mean"]
    assert ra["ari_assigned"] == 1.0
    assert ra["frac_cells_matched"] == 1.0
    assert m["real"]["celladmix"]["mean_total"] == pytest.approx(0.1)  # audit untouched
    assert m["inputs"]["molecules_sha256"]


def test_recompute_unknown_run_exits_2(tmp_path):
    assert recompute.main(["--run", "nope",
                           "--data-root", str(tmp_path)]) == 2
