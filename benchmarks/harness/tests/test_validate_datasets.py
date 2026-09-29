"""Tests for validate_datasets.py (contract validator) and inventory.py."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import common
import fixtures
import inventory as inv
import validate_datasets as vd


# ---------------------------------------------------------------------------
# fixture completion: bring the harness fixtures up to full contract shape
# ---------------------------------------------------------------------------

def _fix_meta(path: Path, kind: str) -> Path:
    """Rewrite a fixture dataset's meta.json so the validator passes cleanly.

    * molecules get contract dtypes (int32 cell ids, float64 coordinates);
    * bbox is fitted around the molecules (area = bbox area);
    * all stats are recomputed from the data;
    * sim sources/truth gain generator_version/seed; real sources gain the
      documented fields;
    * difficulty classes are recomputed from the stats.
    """
    meta = json.loads((path / "meta.json").read_text())
    df = pd.read_parquet(path / "molecules.parquet")
    astype = {"x": np.float64, "y": np.float64}
    if "cell" in df.columns:
        astype["cell"] = np.int32
    df = df.astype(astype)
    df.to_parquet(path / "molecules.parquet", index=False)
    x0, x1 = float(df.x.min() - 1.0), float(df.x.max() + 1.0)
    y0, y1 = float(df.y.min() - 1.0), float(df.y.max() + 1.0)
    area = (x1 - x0) * (y1 - y0)
    meta["crop"]["bbox_um"] = [x0, y0, x1, y1]
    stats = meta.setdefault("stats", {})
    stats["n_molecules"] = int(len(df))
    stats["n_genes"] = int(df["gene"].nunique())
    stats["area_um2"] = area
    stats["molecules_per_um2"] = len(df) / area
    if kind == "sim":
        n_true = int(df.loc[df["cell"] > 0, "cell"].nunique())
        stats.update({"n_true_cells": n_true,
                      "true_cells_per_mm2": n_true / area * 1e6,
                      "n_vendor_cells": 0, "vendor_cells_per_mm2": 0.0})
        meta["source"].update({"generator": "fixtures.make_sim_dataset",
                               "generator_version": 1, "seed": 0})
        meta["platform"] = "simulated"
        truth = meta.setdefault("truth", {})
        truth.setdefault("generator_version", 1)
        truth.setdefault("params", {"seed": 0})
    else:
        raw = df["cell_vendor"].fillna("").astype(str).str.strip()
        n_vendor = int(raw[raw != ""].nunique())
        stats.update({"n_vendor_cells": n_vendor,
                      "vendor_cells_per_mm2": n_vendor / area * 1e6})
        stats.pop("n_true_cells", None)
        stats.pop("true_cells_per_mm2", None)
        meta["source"].update({"url": "fixture", "license": "fixture",
                               "original_dataset": "fixture", "doi": None,
                               "retrieved": "2026-01-01"})
    diff = meta.setdefault("difficulty", {})
    diff["gene_panel"] = vd.gene_panel_class(stats["n_genes"])
    key = "n_true_cells" if kind == "sim" else "n_vendor_cells"
    diff["cell_density"] = vd.cell_density_class(stats[key] / area * 1e6)
    (path / "meta.json").write_text(json.dumps(meta, indent=2))
    return path


def _manifest_entry(path: Path, group: str = "sim") -> dict:
    meta = json.loads((path / "meta.json").read_text())
    return {
        "group": group,
        "path": Path(f"{group}.yaml"),
        "entry": {
            "id": meta["id"],
            "tier": meta["tier"],
            "platform": meta.get("platform"),
            "sha256": {
                "molecules.parquet": common.sha256_file(path / "molecules.parquet"),
                "meta.json": common.sha256_file(path / "meta.json"),
            },
        },
        "baysor": {},
    }


@pytest.fixture()
def sim_dataset(tmp_path: Path) -> Path:
    ds = fixtures.make_sim_dataset(tmp_path / "sim" / "sim_fixture_a")
    return _fix_meta(ds, "sim")


@pytest.fixture()
def real_dataset(tmp_path: Path) -> Path:
    ds = fixtures.make_real_dataset(tmp_path / "real" / "real_fixture_a")
    return _fix_meta(ds, "real")


def _validate(ds: Path, manifests: dict) -> vd.Report:
    report = vd.Report()
    vd.validate_dataset(ds, common.repo_root(), manifests, report)
    return report


# ---------------------------------------------------------------------------
# valid datasets
# ---------------------------------------------------------------------------

def test_valid_sim_dataset_passes(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset, "sim")}
    report = _validate(sim_dataset, manifests)
    assert report.issues == [], [i.message for i in report.issues]
    assert report.errors == 0


def test_valid_real_dataset_passes(real_dataset: Path):
    manifests = {real_dataset.name: _manifest_entry(real_dataset, "real_other")}
    report = _validate(real_dataset, manifests)
    assert report.issues == [], [i.message for i in report.issues]


# ---------------------------------------------------------------------------
# molecule-table (data) findings
# ---------------------------------------------------------------------------

def test_unsorted_rows_is_data_finding(sim_dataset: Path):
    df = pd.read_parquet(sim_dataset / "molecules.parquet")
    df = df.sample(frac=1.0, random_state=0).reset_index(drop=True)
    df.to_parquet(sim_dataset / "molecules.parquet", index=False)
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    report = _validate(sim_dataset, manifests)
    assert any(i.severity == "data" and i.check == "sorted"
               for i in report.issues)
    assert report.errors == 0


def test_wrong_dtype_is_data_finding(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    df = pd.read_parquet(sim_dataset / "molecules.parquet")
    df = df.astype({"y": np.float32})
    df.to_parquet(sim_dataset / "molecules.parquet", index=False)
    report = _validate(sim_dataset, manifests)
    assert any(i.check == "dtypes" and i.severity == "data"
               for i in report.issues)


def test_missing_required_column_is_data_finding(real_dataset: Path):
    df = pd.read_parquet(real_dataset / "molecules.parquet")
    df = df.drop(columns=["cell_vendor"])
    df.to_parquet(real_dataset / "molecules.parquet", index=False)
    manifests = {real_dataset.name: _manifest_entry(real_dataset, "real_other")}
    report = _validate(real_dataset, manifests)
    assert any(i.check == "columns" and "cell_vendor" in i.message
               for i in report.issues)
    # stats now reference a missing column: reported as data, not metadata
    assert all(i.severity == "data" for i in report.issues)


def test_vendor_sentinel_is_data_finding(real_dataset: Path):
    manifests = {real_dataset.name: _manifest_entry(real_dataset, "real_other")}
    df = pd.read_parquet(real_dataset / "molecules.parquet")
    df.loc[df.index[:5], "cell_vendor"] = "0"
    df.to_parquet(real_dataset / "molecules.parquet", index=False)
    report = _validate(real_dataset, manifests)
    assert any(i.check == "content" and "sentinel" in i.message
               and i.severity == "data" for i in report.issues)


def test_control_gene_is_data_finding(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    df = pd.read_parquet(sim_dataset / "molecules.parquet")
    df.loc[df.index[:3], "gene"] = "NegControlProbe_0001"
    df.to_parquet(sim_dataset / "molecules.parquet", index=False)
    report = _validate(sim_dataset, manifests)
    assert any("control/blank" in i.message and i.severity == "data"
               for i in report.issues)


def test_molecules_outside_bbox_is_data_finding(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    df = pd.read_parquet(sim_dataset / "molecules.parquet")
    df.loc[df.index[0], ["x", "y"]] = [1000.0, 1000.0]
    df.to_parquet(sim_dataset / "molecules.parquet", index=False)
    report = _validate(sim_dataset, manifests)
    assert any(i.check == "bbox" and i.severity == "data"
               for i in report.issues)


# ---------------------------------------------------------------------------
# metadata (error) findings
# ---------------------------------------------------------------------------

def test_missing_meta_field_is_error(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    meta = json.loads((sim_dataset / "meta.json").read_text())
    del meta["difficulty"]
    (sim_dataset / "meta.json").write_text(json.dumps(meta))
    report = _validate(sim_dataset, manifests)
    assert any(i.severity == "error" and "difficulty" in i.message
               for i in report.issues)


def test_bad_tier_enum_is_error(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    meta = json.loads((sim_dataset / "meta.json").read_text())
    meta["tier"] = "xlong"
    (sim_dataset / "meta.json").write_text(json.dumps(meta))
    report = _validate(sim_dataset, manifests)
    assert any("tier" in i.message and i.severity == "error"
               for i in report.issues)


def test_stats_mismatch_is_error(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    meta = json.loads((sim_dataset / "meta.json").read_text())
    meta["stats"]["n_molecules"] += 1
    (sim_dataset / "meta.json").write_text(json.dumps(meta))
    report = _validate(sim_dataset, manifests)
    assert any(i.check == "stats" and i.severity == "error"
               for i in report.issues)


def test_panel_class_mismatch_is_error(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    meta = json.loads((sim_dataset / "meta.json").read_text())
    meta["difficulty"]["gene_panel"] = "huge"
    (sim_dataset / "meta.json").write_text(json.dumps(meta))
    report = _validate(sim_dataset, manifests)
    assert any(i.check == "difficulty" and "gene_panel" in i.message
               for i in report.issues)


def test_density_class_mismatch_is_error(real_dataset: Path):
    manifests = {real_dataset.name: _manifest_entry(real_dataset, "real_other")}
    meta = json.loads((real_dataset / "meta.json").read_text())
    meta["difficulty"]["cell_density"] = "sparse"
    (real_dataset / "meta.json").write_text(json.dumps(meta))
    report = _validate(real_dataset, manifests)
    assert any(i.check == "difficulty" and "cell_density" in i.message
               for i in report.issues)


def test_missing_prior_image_is_error(real_dataset: Path):
    manifests = {real_dataset.name: _manifest_entry(real_dataset, "real_other")}
    meta = json.loads((real_dataset / "meta.json").read_text())
    meta["baysor"]["prior"] = "image:images/nope.tif"
    meta["baysor"]["prior_confidence"] = 0.5
    (real_dataset / "meta.json").write_text(json.dumps(meta))
    report = _validate(real_dataset, manifests)
    assert any(i.check == "files" and "prior image" in i.message
               for i in report.issues)


def test_prior_column_missing_is_data_finding(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    meta = json.loads((sim_dataset / "meta.json").read_text())
    meta["baysor"]["prior"] = "column"
    (sim_dataset / "meta.json").write_text(json.dumps(meta))
    report = _validate(sim_dataset, manifests)
    assert any(i.check == "columns" and "prior" in i.message
               and i.severity == "data" for i in report.issues)


def test_missing_image_file_is_error(real_dataset: Path, tmp_path: Path):
    meta = json.loads((real_dataset / "meta.json").read_text())
    meta["images"] = [{"name": "dapi", "file": "images/dapi.tif",
                       "pixel_size_um": 5.0, "origin_um": [0.0, 0.0]}]
    (real_dataset / "meta.json").write_text(json.dumps(meta))
    manifests = {real_dataset.name: _manifest_entry(real_dataset, "real_other")}
    report = _validate(real_dataset, manifests)
    assert any(i.check == "files" and "images[0]" in i.message
               for i in report.issues)


def test_image_covering_bbox_passes(real_dataset: Path):
    import tifffile
    meta = json.loads((real_dataset / "meta.json").read_text())
    bb = meta["crop"]["bbox_um"]
    (real_dataset / "images").mkdir(exist_ok=True)
    tifffile.imwrite(real_dataset / "images" / "dapi.tif",
                     np.zeros((8, 8), dtype=np.uint16))
    meta["images"] = [{"name": "dapi", "file": "images/dapi.tif",
                       "pixel_size_um": 5.0, "origin_um": [bb[0], bb[1]]}]
    (real_dataset / "meta.json").write_text(json.dumps(meta))
    manifests = {real_dataset.name: _manifest_entry(real_dataset, "real_other")}
    report = _validate(real_dataset, manifests)
    assert report.errors == 0, [i.message for i in report.issues]


# ---------------------------------------------------------------------------
# manifest checks
# ---------------------------------------------------------------------------

def test_manifest_sha_mismatch_is_error(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    manifests[sim_dataset.name]["entry"]["sha256"]["meta.json"] = "0" * 64
    report = _validate(sim_dataset, manifests)
    assert any(i.check == "manifest" and "sha256" in i.message
               for i in report.issues)


def test_dataset_not_in_manifest_is_error(sim_dataset: Path):
    report = _validate(sim_dataset, {})
    assert any(i.check == "manifest" and "not listed" in i.message
               for i in report.issues)


def test_manifest_tier_mismatch_is_error(sim_dataset: Path):
    manifests = {sim_dataset.name: _manifest_entry(sim_dataset)}
    manifests[sim_dataset.name]["entry"]["tier"] = "full"
    report = _validate(sim_dataset, manifests)
    assert any("tier" in i.message and i.severity == "error"
               for i in report.issues)


def test_manifest_missing_dataset_is_error(tmp_path: Path):
    root = tmp_path / "data"
    (root / "sim").mkdir(parents=True)
    manifests = {"ghost": {"group": "sim", "path": Path("sim.yaml"),
                           "entry": {"id": "ghost", "tier": "quick"},
                           "baysor": {}}}
    report = vd.validate_all(root, common.repo_root(), manifests=manifests)
    assert any(i.dataset == "ghost" and "missing" in i.message
               for i in report.issues)


# ---------------------------------------------------------------------------
# class helpers (contract boundaries)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n,cls", [
    (0, "tiny"), (49, "tiny"), (50, "small"), (250, "small"),
    (251, "medium"), (700, "medium"), (701, "large"), (2000, "large"),
    (2001, "huge"), (5101, "huge"),
])
def test_gene_panel_class_boundaries(n, cls):
    assert vd.gene_panel_class(n) == cls


@pytest.mark.parametrize("d,cls", [
    (0.0, "sparse"), (2499.9, "sparse"), (2500.0, "medium"),
    (7000.0, "medium"), (7000.1, "dense"), (40000.0, "dense"),
])
def test_cell_density_class_boundaries(d, cls):
    assert vd.cell_density_class(d) == cls


# ---------------------------------------------------------------------------
# CLI exit codes
# ---------------------------------------------------------------------------

def _write_repo_manifests(repo: Path, entries: list[dict]) -> None:
    d = repo / "benchmarks" / "datasets"
    d.mkdir(parents=True, exist_ok=True)
    (d / "sim.yaml").write_text(yaml.safe_dump({"datasets": entries}))


def test_cli_exit_codes(tmp_path: Path, capsys):
    root = tmp_path / "data"
    repo = tmp_path / "repo"
    ds = fixtures.make_sim_dataset(root / "sim" / "sim_fixture_cli")
    ds = _fix_meta(ds, "sim")
    meta = json.loads((ds / "meta.json").read_text())
    _write_repo_manifests(repo, [{
        "id": meta["id"], "tier": meta["tier"],
        "sha256": {
            "molecules.parquet": common.sha256_file(ds / "molecules.parquet"),
            "meta.json": common.sha256_file(ds / "meta.json"),
        },
    }])
    assert vd.main(["--data-root", str(root), "--repo", str(repo)]) == 0

    # data finding: passes without --strict, fails with it
    df = pd.read_parquet(ds / "molecules.parquet")
    df.loc[df.index[:2], "gene"] = "BLANK_0001"
    df.to_parquet(ds / "molecules.parquet", index=False)
    # keep the metadata consistent with the mutated table: only the
    # molecule content (control gene) is the finding under test
    m = json.loads((ds / "meta.json").read_text())
    m["stats"]["n_genes"] = int(df["gene"].nunique())
    (ds / "meta.json").write_text(json.dumps(m))
    assert vd.main(["--data-root", str(root), "--repo", str(repo),
                    "--no-hashes"]) == 0
    assert vd.main(["--data-root", str(root), "--repo", str(repo),
                    "--no-hashes", "--strict"]) == 1

    # metadata error: always fails
    m = json.loads((ds / "meta.json").read_text())
    m["tier"] = "bogus"
    (ds / "meta.json").write_text(json.dumps(m))
    assert vd.main(["--data-root", str(root), "--repo", str(repo)]) == 1
    out = capsys.readouterr()
    assert "metadata error" in out.out


def test_cli_bad_dataset_spec(tmp_path: Path):
    assert vd.main(["--data-root", str(tmp_path), "--datasets", "nope_*"]) == 2


# ---------------------------------------------------------------------------
# inventory
# ---------------------------------------------------------------------------

def _make_repo_with_manifests(repo: Path, entries: list[dict]) -> None:
    d = repo / "benchmarks" / "datasets"
    d.mkdir(parents=True, exist_ok=True)
    (d / "real_other.yaml").write_text(yaml.safe_dump({"datasets": entries}))


def test_inventory_collect_and_render(tmp_path: Path):
    root = tmp_path / "data"
    repo = tmp_path / "repo"
    sim = _fix_meta(fixtures.make_sim_dataset(root / "sim" / "sim_inv_a"), "sim")
    real = _fix_meta(fixtures.make_real_dataset(root / "real" / "real_inv_a"),
                     "real")
    entries = [
        {"id": "sim_inv_a", "tier": "quick", "scenario": "circles_gaps",
         "generator": "trivial"},
        {"id": "real_inv_a", "tier": "quick", "platform": "fixture",
         "tissue": "test tissue"},
    ]
    d = repo / "benchmarks" / "datasets"
    d.mkdir(parents=True, exist_ok=True)
    (d / "sim.yaml").write_text(yaml.safe_dump({"datasets": entries[:1]}))
    (d / "real_other.yaml").write_text(yaml.safe_dump({"datasets": entries[1:]}))

    rows = inv.collect_rows(root, repo)
    assert {r["id"] for r in rows} == {"sim_inv_a", "real_inv_a"}
    text = inv.render(rows, root)
    assert "`sim_inv_a`" in text and "`real_inv_a`" in text
    assert "test tissue" in text
    assert "circles_gaps" in text
    # coverage matrices sum to the row counts
    assert "| **total** |" in text
    # sim row: generator column, no admixture flag
    sim_row = next(r for r in rows if r["id"] == "sim_inv_a")
    assert sim_row["admixture_capable"] == "—"
    assert sim_row["dim"] == "2D"
    real_row = next(r for r in rows if r["id"] == "real_inv_a")
    assert real_row["admixture_capable"] == "no"    # 6 vendor cells < 2000


def test_inventory_matrix_counts():
    rows = [
        {"density": "dense", "panel": "tiny", "kind": "sim", "id": "a",
         "tier": "quick"},
        {"density": "dense", "panel": "tiny", "kind": "sim", "id": "b",
         "tier": "quick"},
        {"density": "sparse", "panel": "huge", "kind": "sim", "id": "c",
         "tier": "full"},
    ]
    out = "\n".join(inv._matrix(rows))
    assert "| **dense** | 2 | 0 | 0 | 0 | 0 | 2 |" in out
    assert "| **sparse** | 0 | 0 | 0 | 0 | 1 | 1 |" in out
    assert "| **total** | 2 | 0 | 0 | 0 | 1 | 3 |" in out


def test_inventory_cli_check(tmp_path: Path):
    root = tmp_path / "data"
    repo = tmp_path / "repo"
    _fix_meta(fixtures.make_sim_dataset(root / "sim" / "sim_inv_b"), "sim")
    d = repo / "benchmarks" / "datasets"
    d.mkdir(parents=True, exist_ok=True)
    (d / "sim.yaml").write_text(yaml.safe_dump({"datasets": [
        {"id": "sim_inv_b", "tier": "quick", "generator": "trivial",
         "scenario": "circles_gaps"}]}))
    out = tmp_path / "DATASETS.md"
    assert inv.main(["--data-root", str(root), "--repo", str(repo),
                     "--out", str(out)]) == 0
    assert out.is_file()
    assert inv.main(["--data-root", str(root), "--repo", str(repo),
                     "--out", str(out), "--check"]) == 0
    out.write_text(out.read_text() + "stale\n")
    assert inv.main(["--data-root", str(root), "--repo", str(repo),
                     "--out", str(out), "--check"]) == 1
