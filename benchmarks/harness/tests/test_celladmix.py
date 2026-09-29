"""Tests for the cellAdmix integration adapter.

The success path is a *contract* test: it runs the real
``benchmarks/celladmix/audit.py`` on a small synthetic fixture and checks that
the adapter reads the audit's actual JSON layout (``metrics.
total_admixture_rate``, ``pairs_top``, counts) — the mocked test that used to
live here hid a layout mismatch. The audit needs the cellAdmix bindings
(``pip install``-ed into the bench env, see ``celladmix/INSTALL.md``); the
contract test skips when they are missing.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import celladmix as camix
import common
from fixtures import make_audit_fixture

REPO = common.repo_root()
AUDIT_SCRIPT = REPO / "benchmarks" / "celladmix" / "audit.py"
TRANSFER_SCRIPT = REPO / "benchmarks" / "celladmix" / "transfer.py"

# audit parameters that make the default-sized fixture (48 cells, 48 genes)
# scorable: the defaults (>= 200 target cells per type) target real crops.
FIXTURE_AUDIT_ARGS = [
    "--min-target-cells", "8", "--min-reference-cells", "5",
    "--n-pool", "12", "--neighbor-k", "10", "--min-excess", "50",
]

# The audit subprocess needs the cellAdmix bindings in this interpreter's
# site-packages (see celladmix/INSTALL.md). Note: find_spec("celladmix") is
# useless here — the harness's own celladmix.py shadows the bindings on
# sys.path — so ask the package metadata instead.
def _bindings_installed() -> bool:
    from importlib import metadata
    try:
        metadata.version("celladmix")
    except metadata.PackageNotFoundError:
        return False
    return True


HAS_CELLADMIX = _bindings_installed()


def _make_repo(tmp_path, script_body):
    script = tmp_path / "benchmarks" / "celladmix" / "audit.py"
    script.parent.mkdir(parents=True)
    script.write_text(script_body)
    return tmp_path


# ---------------------------------------------------------------------------
# absence / command construction
# ---------------------------------------------------------------------------

def test_absent_when_module_missing(tmp_path):
    out = camix.run_audit(tmp_path / "mols.parquet", tmp_path / "a.parquet",
                          tmp_path / "out.json", repo=tmp_path)
    assert out["status"] == camix.STATUS_ABSENT


def test_format_command(tmp_path):
    cmd = camix.format_command("python", tmp_path / "audit.py",
                               tmp_path / "m.parquet", tmp_path / "a.parquet",
                               tmp_path / "o.json")
    assert "--image" not in cmd and "--threads" not in cmd
    cmd = camix.format_command("python", tmp_path / "audit.py",
                               tmp_path / "m.parquet", tmp_path / "a.parquet",
                               tmp_path / "o.json", image=tmp_path / "dapi.tif",
                               threads=4, celltypes=tmp_path / "ct.parquet",
                               fixed_pairs=tmp_path / "pairs.json",
                               save_celltypes=tmp_path / "saved.parquet")
    assert cmd[cmd.index("--image") + 1] == str(tmp_path / "dapi.tif")
    assert cmd[cmd.index("--threads") + 1] == "4"
    assert cmd[cmd.index("--celltypes") + 1] == str(tmp_path / "ct.parquet")
    assert cmd[cmd.index("--fixed-pairs") + 1] == str(tmp_path / "pairs.json")
    assert cmd[cmd.index("--save-celltypes") + 1] == str(tmp_path / "saved.parquet")


# ---------------------------------------------------------------------------
# normalize_audit: the harness-side contract
# ---------------------------------------------------------------------------

def _raw(metrics, **counts):
    return {
        "metrics": {"status": "ok", "total_admixture_rate": 0.123,
                    "total_admixture_molecules": 40.0,
                    "n_pairs_evaluated": 5, "n_pairs_detected": 2,
                    "admixture_capable": False, **metrics},
        "counts": {"n_cells": 48, **counts},
        "parameters": {"typing": "celltypes", "seed": 1, "threads": 2},
        "pairs_top": [{"source": "A", "target": "B", "rate": 0.4}],
    }


def test_normalize_reads_metrics_block_and_pairs_top():
    out = camix.normalize_audit(_raw({}))
    assert out["status"] == camix.STATUS_OK
    assert out["total_admixture_rate"] == pytest.approx(0.123)
    assert out["pairs"][0]["source"] == "A"
    assert out["n_pairs_evaluated"] == 5
    assert out["n_pairs_detected"] == 2
    assert out["n_cells"] == 48
    assert out["admixture_capable"] is False
    assert out["typing"] == "celltypes"


@pytest.mark.parametrize("metrics", [
    {"status": "no_detected_pairs", "total_admixture_rate": None},  # old bug
    {"status": "unavailable", "total_admixture_rate": None},
    {"status": "ok", "total_admixture_rate": None},
    {"status": "ok", "total_admixture_rate": 0.5, "n_pairs_evaluated": 0},
])
def test_normalize_never_yields_zero_for_unscorable_audits(metrics):
    out = camix.normalize_audit(_raw(metrics))
    assert out["status"] == camix.STATUS_UNAVAILABLE
    assert out["total_admixture_rate"] is None
    assert out["reason"]


def test_normalize_rejects_foreign_layout():
    with pytest.raises(ValueError):
        camix.normalize_audit({"something": 1})


# ---------------------------------------------------------------------------
# failure paths of the subprocess wrapper (fake scripts)
# ---------------------------------------------------------------------------

def test_audit_nonzero_exit(tmp_path):
    repo = _make_repo(tmp_path, "import sys; print('boom', file=sys.stderr); sys.exit(3)")
    out = camix.run_audit(tmp_path / "m.parquet", tmp_path / "a.parquet",
                          tmp_path / "o.json", repo=repo)
    assert out["status"] == camix.STATUS_FAILED
    assert out["returncode"] == 3
    assert "boom" in out["stderr_tail"]


def test_audit_bad_output(tmp_path):
    repo = _make_repo(tmp_path, """
import argparse, json
p = argparse.ArgumentParser()
p.add_argument("--molecules"); p.add_argument("--assignment")
p.add_argument("--out"); p.add_argument("--image", default=None)
a = p.parse_args()
json.dump({"something": 1}, open(a.out, "w"))
""")
    out = camix.run_audit(tmp_path / "m.parquet", tmp_path / "a.parquet",
                          tmp_path / "o.json", repo=repo)
    assert out["status"] == camix.STATUS_FAILED
    assert "metrics" in out["reason"]


def test_audit_unparsable_json(tmp_path):
    repo = _make_repo(tmp_path, """
import argparse
p = argparse.ArgumentParser()
p.add_argument("--molecules"); p.add_argument("--assignment")
p.add_argument("--out"); p.add_argument("--image", default=None)
a = p.parse_args()
open(a.out, "w").write("not json {{")
""")
    out = camix.run_audit(tmp_path / "m.parquet", tmp_path / "a.parquet",
                          tmp_path / "o.json", repo=repo)
    assert out["status"] == camix.STATUS_FAILED


# ---------------------------------------------------------------------------
# contract test: the real audit.py on a small fixture
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not HAS_CELLADMIX,
                    reason="cellAdmix bindings not installed (celladmix/INSTALL.md)")
def test_contract_with_real_audit(tmp_path):
    """run_audit() on the real audit.py: the adapter must read
    metrics.total_admixture_rate and pairs_top from the actual JSON."""
    fx = make_audit_fixture(tmp_path / "ds")
    out_json = tmp_path / "audit.json"
    out = camix.run_audit(
        fx["molecules"], fx["assignment"], out_json, repo=REPO,
        threads=2, celltypes=fx["celltypes"], extra_args=FIXTURE_AUDIT_ARGS)
    assert out["status"] == camix.STATUS_OK, out
    assert out_json.is_file()
    raw = json.loads(out_json.read_text())
    # the harness reads the metrics block, not a top-level rate
    assert raw["metrics"]["status"] == "ok"
    assert out["total_admixture_rate"] == pytest.approx(
        raw["metrics"]["total_admixture_rate"])
    assert out["total_admixture_rate"] > 0
    # pair data comes from pairs_top
    assert out["pairs"] == raw["pairs_top"]
    assert raw["pairs_top"], "fixture must yield detected pairs"
    assert out["n_pairs_evaluated"] == raw["metrics"]["n_pairs_evaluated"] > 0
    assert out["n_pairs_detected"] == raw["metrics"]["n_pairs_detected"] > 0
    assert out["n_cells"] == raw["counts"]["n_cells"] == 48
    # 48 cells: the audit runs, but the crop is not admixture-capable
    assert out["admixture_capable"] is False
    # --threads and --celltypes reached the audit
    assert out["audit_threads"] == 2
    assert out["typing"] == "celltypes"
    assert Path(out["audit_json"]) == out_json


@pytest.mark.skipif(not HAS_CELLADMIX,
                    reason="cellAdmix bindings not installed (celladmix/INSTALL.md)")
def test_contract_save_celltypes_and_unavailable_fallback(tmp_path):
    """--save-celltypes writes the typing used, and an audit run on a
    too-small crop degrades to status=unavailable with rate None."""
    fx = make_audit_fixture(tmp_path / "ds")
    saved = tmp_path / "saved_types.parquet"
    # default thresholds (>= 200 target cells): no pair survives -> unavailable
    out = camix.run_audit(fx["molecules"], fx["assignment"],
                          tmp_path / "audit.json", repo=REPO, threads=2,
                          celltypes=fx["celltypes"], save_celltypes=saved)
    assert out["status"] == camix.STATUS_UNAVAILABLE
    assert out["total_admixture_rate"] is None
    assert saved.is_file()
    types = pd.read_parquet(saved)
    assert set(types.columns) == {"cell", "celltype"}
    assert len(types) == 48


# ---------------------------------------------------------------------------
# transfer wrapper (real transfer.py)
# ---------------------------------------------------------------------------

def test_run_transfer_from_baseline_assignment(tmp_path):
    if not TRANSFER_SCRIPT.is_file():
        pytest.skip("transfer.py missing")
    fx = make_audit_fixture(tmp_path / "ds")
    # perturb the target assignment: merge every pair of cells 1..8 -> new ids
    assign = pd.read_parquet(fx["assignment"])
    target = assign.copy()
    remap = {1: 101, 2: 101, 3: 102, 4: 102}
    target["cell"] = target["cell"].map(lambda c: remap.get(c, c)).astype(np.int64)
    target_path = tmp_path / "target.parquet"
    target.to_parquet(target_path, index=False)

    out_ct = tmp_path / "transferred.parquet"
    report = tmp_path / "transfer.json"
    stats = camix.run_transfer(fx["molecules"], fx["assignment"],
                               fx["celltypes"], target_path, out_ct,
                               report=report, repo=REPO)
    assert stats["status"] == camix.STATUS_OK, stats
    assert report.is_file() and out_ct.is_file()
    assert stats["n_typed_target_cells"] > 0
    got = pd.read_parquet(out_ct)
    # merged cell 101 got the majority baseline type of cells 1 and 2
    by_id = dict(zip(got["cell"].astype(str), got["celltype"]))
    assert by_id["101"] in ("type0",)
    assert by_id.get("48") in ("type3",)


@pytest.mark.skipif(not HAS_CELLADMIX,
                    reason="cellAdmix bindings not installed (celladmix/INSTALL.md)")
def test_contract_fixed_pair_scoring(tmp_path):
    """Baseline audit -> fixed pair set -> every later audit scores the same
    pairs (detection-independent), through the real audit.py."""
    fx = make_audit_fixture(tmp_path / "ds")
    base_out = tmp_path / "base.json"
    base = camix.run_audit(
        fx["molecules"], fx["assignment"], base_out, repo=REPO, threads=2,
        celltypes=fx["celltypes"], extra_args=FIXTURE_AUDIT_ARGS,
        save_celltypes=tmp_path / "saved.parquet")
    assert base["status"] == camix.STATUS_OK, base
    assert base["pairs"], "baseline audit must detect pairs"
    fixed = tmp_path / "fixed_pairs.json"
    fixed.write_text(json.dumps([{"source": p["source"], "target": p["target"]}
                                 for p in base["pairs"]]))

    out_json = tmp_path / "run.json"
    out = camix.run_audit(
        fx["molecules"], fx["assignment"], out_json, repo=REPO, threads=2,
        celltypes=fx["celltypes"], fixed_pairs=fixed,
        extra_args=FIXTURE_AUDIT_ARGS)
    assert out["status"] == camix.STATUS_OK, out
    assert out["n_pairs_evaluated"] == len(base["pairs"])
    assert out["total_admixture_rate"] is not None
    assert out["total_admixture_rate"] > 0
    raw = json.loads(out_json.read_text())
    # provenance of the fixed set + the relaxed marker-pool gate (patch 0004)
    assert raw["inputs"]["fixed_pairs_sha256"] == common.sha256_file(fixed)
    assert raw["parameters"]["audit"]["min_pool_markers"] == 1
