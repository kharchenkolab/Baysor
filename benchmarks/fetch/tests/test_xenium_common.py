"""Tests for the BENCH-REALX Xenium fetch/crop helpers.

Run from the repository root:

    .deps/bench/bin/python -m pytest benchmarks/fetch/tests
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import xenium_common as xc  # noqa: E402


# ---------------------------------------------------------------------------
# gene filtering
# ---------------------------------------------------------------------------


def _chunk(values, **cols) -> pa.Table:
    data = {"feature_name": pa.array(values, pa.string())}
    data.update({k: pa.array(v) for k, v in cols.items()})
    return pa.table(data)


def test_gene_mask_name_fallback_drops_controls():
    names = [
        "ACTB", "NegControlCodeword_0502", "NegControlProbe_0001",
        "UnassignedCodeword_0013", "DeprecatedCodeword_0007",
        "Intergenic_Region_0001", "BLANK_0006", "blank_0007", "GAPDH",
    ]
    mask = xc.gene_mask(_chunk(names))
    assert mask.tolist() == [True, False, False, False, False, False, False, False, True]


def test_gene_mask_prefers_is_gene():
    names = ["ACTB", "BLANK_0006"]  # would be dropped by the name filter
    tbl = _chunk(names, is_gene=[True, True])
    assert xc.gene_mask(tbl).all()


def test_gene_mask_codeword_category():
    tbl = _chunk(
        ["A", "B", "C"],
        codeword_category=["predesigned_gene", "custom_gene", "blank_codeword"],
    )
    assert xc.gene_mask(tbl).tolist() == [True, True, False]


def test_gene_mask_binary_feature_names():
    tbl = pa.table(
        {
            "feature_name": pa.array([b"ACTB", b"BLANK_0006"], pa.binary()),
        }
    )
    assert xc.gene_mask(tbl).tolist() == [True, False]


def test_filter_mask_combines_qv_and_gene():
    tbl = _chunk(["ACTB", "ACTB", "BLANK_0001"], qv=[30.0, 10.0, 30.0])
    assert xc.filter_mask(tbl, min_qv=20).tolist() == [True, False, False]


# ---------------------------------------------------------------------------
# prior / cell_vendor
# ---------------------------------------------------------------------------


def test_assign_vendor_prior_numeric_ids():
    ids = np.array([0, 5, 7, 3], dtype=np.int32)
    ov = np.array([0, 1, 0, 1], dtype=np.uint8)
    vendor, prior = xc.assign_vendor_prior(ids, ov)
    assert vendor.tolist() == ["", "5", "7", "3"]
    assert prior.dtype == np.int32
    assert prior.tolist() == [0, 5, 0, 3]  # overlaps_nucleus == 1 required


def test_assign_vendor_prior_string_ids_are_deterministic():
    ids = np.array(["UNASSIGNED", "bbb-1", "UNASSIGNED", "aaa-1", "bbb-1"], dtype=object)
    ov = np.array([0, 1, 1, 1, 0], dtype=np.uint8)
    vendor, prior = xc.assign_vendor_prior(ids, ov)
    assert vendor.tolist() == ["", "bbb-1", "", "aaa-1", "bbb-1"]
    # sorted unique nucleus-overlapping ids: aaa-1 -> 1, bbb-1 -> 2
    assert prior.tolist() == [0, 2, 0, 1, 0]


def test_assign_vendor_prior_row_order_independence():
    ids = np.array(["z-1", "a-1", "z-1"], dtype=object)
    ov = np.ones(3, dtype=np.uint8)
    _, p1 = xc.assign_vendor_prior(ids, ov)
    order = [1, 0, 2]
    _, p2 = xc.assign_vendor_prior(ids[order], ov[order])
    assert p1[order].tolist() == p2.tolist()


# ---------------------------------------------------------------------------
# molecule table
# ---------------------------------------------------------------------------


def _crop_table(n=10, seed=0):
    rng = np.random.default_rng(seed)
    return pa.table(
        {
            "x_location": pa.array(rng.uniform(0, 100, n), pa.float32()),
            "y_location": pa.array(rng.uniform(0, 100, n), pa.float32()),
            "qv": pa.array(rng.uniform(20, 40, n), pa.float32()),
            "feature_name": pa.array(rng.choice(["A", "B"], n).astype(str), pa.string()),
            "cell_id": pa.array(rng.integers(0, 5, n).astype(np.int32), pa.int32()),
            "overlaps_nucleus": pa.array(rng.integers(0, 2, n).astype(np.uint8), pa.uint8()),
        }
    )


def test_build_molecule_table_contract_and_sorting():
    t = xc.build_molecule_table(_crop_table(50))
    assert t.column_names == ["x", "y", "gene", "qv", "prior", "cell_vendor"]
    assert t.schema.field("x").type == pa.float64()
    assert t.schema.field("qv").type == pa.float32()
    assert t.schema.field("prior").type == pa.int32()
    y = np.asarray(t["y"].to_numpy(zero_copy_only=False))
    x = np.asarray(t["x"].to_numpy(zero_copy_only=False))
    assert np.all(np.diff(y) >= 0)
    # sorted lexicographically by (y, x)
    keys = list(zip(y.tolist(), x.tolist()))
    assert keys == sorted(keys)
    prior = np.asarray(t["prior"].to_numpy(zero_copy_only=False))
    import pyarrow.compute as pc

    src = _crop_table(50)
    idx = pc.sort_indices(
        src, sort_keys=[("y_location", "ascending"), ("x_location", "ascending")]
    )
    ov = np.asarray(src.take(idx)["overlaps_nucleus"].to_numpy(zero_copy_only=False))
    assert (prior[ov == 0] == 0).all()


# ---------------------------------------------------------------------------
# transcript reading
# ---------------------------------------------------------------------------


@pytest.fixture
def transcripts_v1(tmp_path):
    """A V1-style transcripts table (no is_gene/codeword_category)."""
    tbl = pa.table(
        {
            "transcript_id": pa.array(np.arange(12, dtype=np.uint64), pa.uint64()),
            "cell_id": pa.array(
                np.array([0, 1, 2, 0] * 3, dtype=np.int32), pa.int32()
            ),
            "overlaps_nucleus": pa.array(
                np.array([0, 1, 0, 1] * 3, dtype=np.uint8), pa.uint8()
            ),
            "feature_name": pa.array(
                ["ACTB", "BLANK_0001", "NegControlCodeword_0500", "GAPDH"] * 3,
                pa.string(),
            ),
            "x_location": pa.array(
                [5.0, 50.0, 5.0, 50.0, 150.0, 150.0, 5.0, 50.0, 150.0, 5.0, 50.0, 150.0],
                pa.float32(),
            ),
            "y_location": pa.array(
                [5.0, 5.0, 50.0, 50.0, 50.0, 150.0, 150.0, 150.0, 150.0, 5.0, 5.0, 5.0],
                pa.float32(),
            ),
            "qv": pa.array(
                [30.0, 30.0, 30.0, 15.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0],
                pa.float32(),
            ),
        }
    )
    path = tmp_path / "transcripts.parquet"
    pq.write_table(tbl, path, row_group_size=4)  # force multiple row groups
    return path


def test_read_transcript_crops_filters(transcripts_v1):
    (t,) = xc.read_transcript_crops(transcripts_v1, [(0, 0, 100, 100)], min_qv=20)
    genes = set(t["feature_name"].to_pylist())
    # BLANK/NegControl dropped, qv=15 row dropped, out-of-box rows dropped
    assert genes <= {"ACTB", "GAPDH"}
    x = np.asarray(t["x_location"].to_numpy(zero_copy_only=False))
    y = np.asarray(t["y_location"].to_numpy(zero_copy_only=False))
    assert ((x >= 0) & (x < 100) & (y >= 0) & (y < 100)).all()
    # qv filter applied
    assert min(t["qv"].to_pylist()) >= 20


def test_read_transcript_crops_multiple_boxes_single_pass(transcripts_v1):
    a, b = xc.read_transcript_crops(
        transcripts_v1, [(0, 0, 10, 10), (140, 140, 160, 160)], min_qv=20
    )
    assert a.num_rows >= 1
    assert b.num_rows >= 1
    xa = np.asarray(a["x_location"].to_numpy(zero_copy_only=False))
    assert (xa < 10).all()
    xb = np.asarray(b["x_location"].to_numpy(zero_copy_only=False))
    assert (xb >= 140).all()


def test_transcript_histogram_shape(transcripts_v1):
    hist = xc.transcript_histogram(transcripts_v1, min_qv=20, bounds=(0, 0, 200, 200), bin_um=50)
    assert hist.shape == (4, 4)
    assert hist.sum() > 0


def test_transcript_histogram_total_matches_reader(transcripts_v1):
    hist = xc.transcript_histogram(transcripts_v1, min_qv=20, bounds=(0, 0, 200, 200), bin_um=50)
    (t,) = xc.read_transcript_crops(transcripts_v1, [(0, 0, 200, 200)], min_qv=20)
    assert hist.sum() == t.num_rows


# ---------------------------------------------------------------------------
# crop selection
# ---------------------------------------------------------------------------


def _synthetic_hists():
    """100 x 100 bins (25 um), a dense patch in the middle and sparse rest."""
    mol = np.full((100, 100), 5, dtype=np.int64)
    cells = np.full((100, 100), 1, dtype=np.int64)
    mol[30:70, 30:70] = 40
    cells[30:70, 30:70] = 6
    mol[0:20, 0:20] = 0  # empty corner
    cells[0:20, 0:20] = 0
    return mol, cells


def test_pick_bbox_respects_limits():
    mol, cells = _synthetic_hists()
    res = xc.pick_bbox(
        mol, cells, 25.0, (0, 0, 2500, 2500),
        target=20000, min_mols=5000, max_mols=30000,
        min_cells=50, min_coverage=0.7,
    )
    x0, y0, x1, y1 = res["bbox_um"]
    assert 5000 <= res["mols"] <= 30000
    assert res["cells"] >= 50
    assert res["coverage"] >= 0.7
    assert x1 > x0 and y1 > y0


def test_pick_bbox_deterministic():
    mol, cells = _synthetic_hists()
    kw = dict(
        target=20000, min_mols=5000, max_mols=30000,
        min_cells=50, min_coverage=0.7,
    )
    r1 = xc.pick_bbox(mol, cells, 25.0, (0, 0, 2500, 2500), **kw, seed=0)
    r2 = xc.pick_bbox(mol, cells, 25.0, (0, 0, 2500, 2500), **kw, seed=0)
    assert r1 == r2


def test_pick_bbox_within_containment():
    mol, cells = _synthetic_hists()
    within = (750.0, 750.0, 1750.0, 1750.0)
    res = xc.pick_bbox(
        mol, cells, 25.0, (0, 0, 2500, 2500),
        target=10000, min_mols=1000, max_mols=30000,
        min_cells=10, min_coverage=0.5, within=within,
    )
    x0, y0, x1, y1 = res["bbox_um"]
    assert x0 >= within[0] and y0 >= within[1]
    assert x1 <= within[2] and y1 <= within[3]


def test_pick_bbox_density_hint_dense():
    mol, cells = _synthetic_hists()
    res = xc.pick_bbox(
        mol, cells, 25.0, (0, 0, 2500, 2500),
        target=20000, min_mols=5000, max_mols=30000,
        min_cells=50, min_coverage=0.7, density_hint="dense",
        dense_threshold=4000.0,
    )
    # cells/bin -> cells/mm^2: 6 per 625 um^2 = 9600/mm^2
    assert res["cells_per_mm2"] >= 4000.0


def test_pick_bbox_raises_when_impossible():
    mol, cells = _synthetic_hists()
    with pytest.raises(RuntimeError):
        xc.pick_bbox(
            mol, cells, 25.0, (0, 0, 2500, 2500),
            target=10**9, min_mols=10**8, max_mols=2 * 10**8,
            min_cells=10, min_coverage=0.5,
        )


# ---------------------------------------------------------------------------
# classes and scale
# ---------------------------------------------------------------------------


def test_density_class_boundaries():
    assert xc.density_class(2499.9) == "sparse"
    assert xc.density_class(2500) == "medium"
    assert xc.density_class(7000) == "medium"
    assert xc.density_class(7000.1) == "dense"


def test_gene_panel_class_boundaries():
    assert xc.gene_panel_class(49) == "tiny"
    assert xc.gene_panel_class(50) == "small"
    assert xc.gene_panel_class(249) == "small"
    assert xc.gene_panel_class(250) == "medium"
    assert xc.gene_panel_class(699) == "medium"
    assert xc.gene_panel_class(700) == "large"
    assert xc.gene_panel_class(1999) == "large"
    assert xc.gene_panel_class(2000) == "huge"


def test_estimate_scale_um_nucleus_method():
    import pandas as pd

    n = 100
    cells = pd.DataFrame(
        {
            "x_centroid": np.full(n, 50.0),
            "y_centroid": np.full(n, 50.0),
            "nucleus_area": np.full(n, math.pi * 9.0),  # radius 3 um
            "cell_area": np.full(n, 100.0),
        }
    )
    scale, method = xc.estimate_scale_um(cells, (0, 0, 100, 100))
    assert scale == pytest.approx(4.5, abs=0.01)  # 1.5 * 3
    assert "nucleus_area" in method


def test_estimate_scale_um_cell_fallback():
    import pandas as pd

    n = 100
    cells = pd.DataFrame(
        {
            "x_centroid": np.full(n, 50.0),
            "y_centroid": np.full(n, 50.0),
            "nucleus_area": np.full(n, np.nan),
            "cell_area": np.full(n, math.pi * 16.0),  # radius 4 um
        }
    )
    scale, method = xc.estimate_scale_um(cells, (0, 0, 100, 100))
    assert scale == pytest.approx(6.0, abs=0.01)
    assert "cell_area" in method


def test_panel_gene_count(tmp_path):
    panel = {
        "payload": {
            "targets": [
                {"type": {"descriptor": "gene"}},
                {"type": {"descriptor": "gene"}},
                {"type": {"descriptor": "negative_control"}},
            ]
        }
    }
    p = tmp_path / "gene_panel.json"
    p.write_text(json.dumps(panel))
    assert xc.panel_gene_count(p) == 2


# ---------------------------------------------------------------------------
# boundaries
# ---------------------------------------------------------------------------


def test_crop_boundaries_keeps_intersecting_objects(tmp_path):
    # object 1 fully inside, object 2 straddling, object 3 far away
    tbl = pa.table(
        {
            "cell_id": pa.array([1, 1, 1, 2, 2, 2, 3, 3], pa.int32()),
            "vertex_x": pa.array([10.0, 20.0, 15.0, 90.0, 110.0, 100.0, 500.0, 510.0]),
            "vertex_y": pa.array([10.0, 20.0, 5.0, 10.0, 10.0, 20.0, 500.0, 510.0]),
        }
    )
    p = tmp_path / "b.parquet"
    pq.write_table(tbl, p)
    out = xc.crop_boundaries(p, (0.0, 0.0, 100.0, 100.0))
    assert set(out["cell_id"].to_pylist()) == {1, 2}
    # untouched file keeps all objects
    assert set(pq.read_table(p)["cell_id"].to_pylist()) == {1, 2, 3}


def test_crop_boundaries_string_ids(tmp_path):
    tbl = pa.table(
        {
            "cell_id": pa.array(["a-1", "a-1", "b-1", "b-1"], pa.string()),
            "vertex_x": pa.array([10.0, 20.0, 400.0, 420.0]),
            "vertex_y": pa.array([10.0, 20.0, 400.0, 420.0]),
        }
    )
    p = tmp_path / "b.parquet"
    pq.write_table(tbl, p)
    out = xc.crop_boundaries(p, (0.0, 0.0, 100.0, 100.0))
    assert out["cell_id"].to_pylist() == ["a-1", "a-1"]


# ---------------------------------------------------------------------------
# images
# ---------------------------------------------------------------------------


def test_ome_info_and_window_read(tmp_path):
    import tifffile

    rng = np.random.default_rng(0)
    img = rng.integers(0, 1000, size=(2, 64, 128)).astype(np.uint16)
    p = tmp_path / "focus.ome.tif"
    tifffile.imwrite(
        str(p), img, ome=True,
        metadata={
            "axes": "CYX",
            "PhysicalSizeX": 0.5,
            "PhysicalSizeY": 0.5,
            "Channel": {"Name": ["DAPI", "Membrane"]},
        },
    )
    info = xc.ome_info(p)
    assert info["shape"] == (2, 64, 128)
    assert info["pixel_size_um"] == 0.5
    assert info["channels"][0] == "DAPI"
    # crop window (5..10 um x, 1..4 um y) -> px [10,20) x [2,8)
    win = xc.read_image_window(p, 0, (5.0, 1.0, 10.0, 4.0), 0.5)
    np.testing.assert_array_equal(win, img[0, 2:8, 10:20])


def test_write_tif_roundtrip(tmp_path):
    arr = np.arange(100, dtype=np.uint16).reshape(10, 10)
    p = tmp_path / "out" / "x.tif"
    xc.write_tif(p, arr)
    import tifffile

    np.testing.assert_array_equal(tifffile.imread(str(p)), arr)


# ---------------------------------------------------------------------------
# meta
# ---------------------------------------------------------------------------


def test_build_meta_contract_shape():
    stats = {
        "n_molecules": 1000,
        "n_genes": 50,
        "area_um2": 1e6,
        "molecules_per_um2": 0.001,
        "n_vendor_cells": 9000,
        "vendor_cells_per_mm2": 9000.0,
    }
    meta = xc.build_meta(
        dataset_id="x",
        tier="quick",
        source={"url": "u", "license": "L", "original_dataset": "O", "doi": None},
        bbox=(0, 0, 1000, 1000),
        note="n",
        stats=stats,
        panel_genes=3000,
        scale_um=4.5,
        scale_method="m",
        baysor={"config": "configs/xenium.toml", "prior": "column",
                "prior_confidence": 0.5, "scale_std": "25%",
                "min_molecules_per_cell": 50, "extra_args": []},
        images=[],
        retrieved="2026-09-29",
    )
    for key in ("id", "kind", "tier", "platform", "source", "crop", "stats",
                "difficulty", "baysor", "images", "truth"):
        assert key in meta
    assert meta["kind"] == "real"
    assert meta["truth"] is None
    assert meta["difficulty"] == {
        "cell_density": "dense",
        "gene_panel": "huge",
        "notes": "n",
    }
    assert meta["baysor"]["prior"] == "column"
    assert meta["baysor"]["scale_um"] == 4.5


def test_molecule_stats():
    t = xc.build_molecule_table(_crop_table(100))
    stats = xc.molecule_stats(t, (0, 0, 100, 100), n_vendor_cells=50)
    assert stats["n_molecules"] == 100
    assert stats["area_um2"] == 10000.0
    assert stats["n_genes"] == len(set(t["gene"].to_pylist()))
    assert stats["vendor_cells_per_mm2"] == pytest.approx(50 * 1e6 / 1e4)
