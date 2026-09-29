"""Tests for baseline creation (copies metrics, keeps assignments outside git)."""
import numpy as np
import pandas as pd
import pytest

import baseline
import common
from fixtures import make_sim_dataset, make_real_dataset, make_run, corrupt


def _sim_run(root, ds_id="sim_a", n_reps=3, degrade=False):
    ds_dir = make_sim_dataset(root / "sim" / ds_id)
    molecules = pd.read_parquet(ds_dir / "molecules.parquet")
    base = molecules["cell"].to_numpy(np.int64)
    reps = [corrupt(base, 0.2, seed=10 + k) if degrade else base.copy()
            for k in range(n_reps)]
    make_run(root, "run1", ds_dir, reps)
    return ds_dir


def test_create_copies_metrics_and_assignments(tmp_path):
    root = tmp_path / "data"
    _sim_run(root)
    baselines = tmp_path / "baselines"
    rc = baseline.create("run1", "smoke", root, baselines)
    assert rc == 0
    out = baselines / "smoke" / "sim_a.json"
    assert out.is_file()
    m = common.read_json(out)
    assert m["run_id"] == "run1"
    assert m["baseline"]["name"] == "smoke"
    assert m["baseline"]["noise_floor_replicates"] == 3
    assert m["baseline"]["noise_floor_valid"] is True
    # assignment tables live under the data root, not the repo
    for k in range(3):
        p = root / "baselines" / "smoke" / "sim_a" / f"rep{k}" / "assignment.parquet"
        assert p.is_file()
        sha = common.sha256_file(p)
        assert m["baseline"]["assignment_sha256"][str(k)] == sha
    # sim metrics were carried over
    assert m["sim"]["mean"]["matched_accuracy"] == pytest.approx(1.0)
    assert m["sim"]["sd"]["matched_accuracy"] == 0.0


def test_create_refuses_fewer_than_three_replicates(tmp_path):
    root = tmp_path / "data"
    _sim_run(root, n_reps=1)
    rc = baseline.create("run1", "thin", root, tmp_path / "baselines")
    assert rc == 2
    assert not (tmp_path / "baselines" / "thin").exists()


def test_create_allow_incomplete(tmp_path):
    root = tmp_path / "data"
    _sim_run(root, n_reps=1)
    rc = baseline.create("run1", "thin", root, tmp_path / "baselines",
                         allow_incomplete=True)
    assert rc == 0
    m = common.read_json(tmp_path / "baselines" / "thin" / "sim_a.json")
    assert m["baseline"]["noise_floor_valid"] is False


def test_create_overwrite_requires_force(tmp_path):
    root = tmp_path / "data"
    _sim_run(root)
    baselines = tmp_path / "baselines"
    assert baseline.create("run1", "smoke", root, baselines) == 0
    assert baseline.create("run1", "smoke", root, baselines) == 2
    assert baseline.create("run1", "smoke", root, baselines, force=True) == 0


def test_create_missing_run_and_bad_name(tmp_path):
    root = tmp_path / "data"
    _sim_run(root)
    assert baseline.create("nope", "x", root, tmp_path / "b") == 2
    assert baseline.create("run1", "../evil", root, tmp_path / "b") == 2


def test_create_real_dataset(tmp_path):
    root = tmp_path / "data"
    ds_dir = make_real_dataset(root / "real" / "real_a")
    cells = pd.read_parquet(ds_dir / "molecules.parquet")["cell"].to_numpy(np.int64)
    make_run(root, "rreal", ds_dir, [cells, cells, cells])
    rc = baseline.create("rreal", "realsmoke", root, tmp_path / "baselines")
    assert rc == 0
    m = common.read_json(tmp_path / "baselines" / "realsmoke" / "real_a.json")
    # replicate agreement on identical reps is perfect
    assert m["real"]["rep_agreement"]["mean"]["molecule_ari"] == pytest.approx(1.0)


def test_list_baselines(tmp_path, capsys):
    root = tmp_path / "data"
    _sim_run(root)
    baselines = tmp_path / "baselines"
    baseline.create("run1", "smoke", root, baselines)
    assert baseline.list_baselines(baselines) == 0
    out = capsys.readouterr().out
    assert "smoke" in out and "sim_a" in out
