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


# ---------------------------------------------------------------------------
# successful-replicate counting, identical flavour, atomic swap, typing
# ---------------------------------------------------------------------------

def test_failed_replicates_do_not_count(tmp_path):
    """The noise floor must be built from *successful* replicates only."""
    root = tmp_path / "data"
    ds_dir = make_sim_dataset(root / "sim" / "sim_ok")
    cells = pd.read_parquet(ds_dir / "molecules.parquet")["cell"].to_numpy(np.int64)
    # three replicate records, two failed -> only 1 successful
    make_run(root, "runmix", ds_dir, [cells, cells, cells],
             statuses=["ok", "failed", "failed"])
    rc = baseline.create("runmix", "mix", root, tmp_path / "baselines")
    assert rc == 2
    assert not (tmp_path / "baselines" / "mix").exists()


def test_identical_flavour_requires_one_thread(tmp_path):
    root = tmp_path / "data"
    ds_dir = make_sim_dataset(root / "sim" / "sim_id")
    cells = pd.read_parquet(ds_dir / "molecules.parquet")["cell"].to_numpy(np.int64)

    make_run(root, "run6t", ds_dir, [cells], threads=6)
    rc = baseline.create("run6t", "ident", root, tmp_path / "baselines",
                         identical=True)
    assert rc == 2                       # 6 threads are not bitwise-stable
    assert not (tmp_path / "baselines" / "ident").exists()

    make_run(root, "run1t", ds_dir, [cells], threads=1)
    rc = baseline.create("run1t", "ident", root, tmp_path / "baselines",
                         identical=True)
    assert rc == 0                       # 1 thread, 1 replicate is enough
    m = common.read_json(tmp_path / "baselines" / "ident" / "sim_id.json")
    assert m["baseline"]["flavour"] == "identical"
    assert m["baseline"]["identical"] is True
    assert m["baseline"]["noise_floor_replicates"] == 1
    assert m["baseline"]["noise_floor_valid"] is False


def test_force_error_keeps_old_baseline(tmp_path):
    """A failing --force must neither delete the old baseline nor leave
    stale temp files behind."""
    root = tmp_path / "data"
    _sim_run(root)
    baselines = tmp_path / "baselines"
    assert baseline.create("run1", "keepme", root, baselines) == 0
    old_json = baselines / "keepme" / "sim_a.json"
    old_data = root / "baselines" / "keepme" / "sim_a" / "rep0" / "assignment.parquet"
    assert old_json.is_file() and old_data.is_file()
    old_sha = common.sha256_file(old_json)

    # break the run: remove one replicate's assignment table
    (root / "runs" / "run1" / "sim_a" / "rep2" / "assignment.parquet").unlink()
    rc = baseline.create("run1", "keepme", root, baselines, force=True)
    assert rc == 2
    # old baseline intact, no temp/old dirs left anywhere
    assert old_json.is_file()
    assert common.sha256_file(old_json) == old_sha
    assert old_data.is_file()
    leftovers = [d.name for d in baselines.iterdir() if d.name.startswith(".")]
    leftovers += [d.name for d in (root / "baselines").iterdir()
                  if d.name.startswith(".")]
    assert leftovers == []


def test_baseline_stores_celltypes_and_fixed_pairs(tmp_path):
    import json
    root = tmp_path / "data"
    ds_dir = make_real_dataset(root / "real" / "real_typ")
    cells = pd.read_parquet(ds_dir / "molecules.parquet")["cell"].to_numpy(np.int64)
    make_run(root, "runtyp", ds_dir, [cells, cells, cells])
    # simulate what run.py's audit wrote for replicate 0
    rep0 = root / "runs" / "runtyp" / "real_typ" / "rep0"
    pd.DataFrame({"cell": ["1", "2", "3"],
                  "celltype": ["a", "b", "a"]}).to_parquet(
        rep0 / "celltypes.parquet", index=False)
    with open(rep0 / "celladmix.json", "w") as fh:
        json.dump({"metrics": {"status": "ok"},
                   "pairs_top": [{"source": "a", "target": "b", "rate": 0.1},
                                 {"source": "b", "target": "a", "rate": 0.2}]},
                  fh)

    baselines = tmp_path / "baselines"
    assert baseline.create("runtyp", "typed", root, baselines) == 0
    m = common.read_json(baselines / "typed" / "real_typ.json")
    data = root / "baselines" / "typed" / "real_typ"
    assert m["baseline"]["celltypes_sha256"] == \
        common.sha256_file(data / "celltypes.parquet")
    assert m["baseline"]["fixed_pairs_sha256"] == \
        common.sha256_file(data / "fixed_pairs.json")
    assert m["baseline"]["n_fixed_pairs"] == 2
    pairs = json.loads((data / "fixed_pairs.json").read_text())
    assert pairs == [{"source": "a", "target": "b"},
                     {"source": "b", "target": "a"}]


def test_baseline_without_audit_records_no_typing(tmp_path):
    root = tmp_path / "data"
    _sim_run(root)
    baselines = tmp_path / "baselines"
    assert baseline.create("run1", "plain", root, baselines) == 0
    m = common.read_json(baselines / "plain" / "sim_a.json")
    assert m["baseline"]["celltypes_sha256"] is None
    assert m["baseline"]["fixed_pairs_sha256"] is None
