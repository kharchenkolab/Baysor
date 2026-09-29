"""Tests for benchmarks/fetch/other.py + other_utils.py (BENCH-REALO).

Run with:

    .deps/bench/bin/python -m pytest benchmarks/fetch/tests
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import other as O  # noqa: E402
import other_utils as U  # noqa: E402


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_solve_anubis_pow_matches_difficulty():
    nonce, digest = U.solve_anubis_pow("a" * 64, difficulty=4)
    assert digest.startswith("0000")
    assert digest == __import__("hashlib").sha256(
        ("a" * 64 + str(nonce)).encode()
    ).hexdigest()
    # deterministic
    assert U.solve_anubis_pow("a" * 64, difficulty=4) == (nonce, digest)


def test_solve_anubis_pow_odd_difficulty():
    _, digest = U.solve_anubis_pow("challenge", difficulty=3)
    # 1 zero byte + zero high nibble
    assert digest.startswith("000")
    assert int(digest[2], 16) < 16 or digest[2] == "0"


def test_pat_to_regex_and_gene_mask():
    assert U.pat_to_regex("Blank*") == "^Blank.*$"
    genes = pd.Series(["Blank1", "GeneA", "Negative2", "GeneB", "NEGR1"])
    keep = U.gene_mask(genes, ["Blank*", "Negative*"])
    assert list(keep) == [False, True, False, True, True]
    # no patterns -> everything kept
    assert U.gene_mask(genes, []).all()
    # NEGR1 must survive a NegControl* pattern (anchored match)
    assert list(U.gene_mask(pd.Series(["NEGR1", "NegControl1"]), ["NegControl*"])) == [True, False]


def test_strip_hybridization_suffix():
    s = pd.Series(["Acta2_Hybridization5", "Cnr1", "Vip_Hybridization2"])
    assert list(U.strip_hybridization_suffix(s)) == ["Acta2", "Cnr1", "Vip"]


def test_densest_window_is_deterministic_and_under_cap():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(0, 1, 40_000) + 50, rng.normal(100, 8, 400_000) + 50])
    y = np.concatenate([rng.normal(0, 1, 40_000) + 50, rng.normal(100, 8, 400_000) + 50])
    w1 = U.densest_window(x, y, target=150_000, cap=150_000, seed=7)
    w2 = U.densest_window(x, y, target=150_000, cap=150_000, seed=7)
    assert w1 == w2
    n = ((x >= w1[0]) & (x <= w1[2]) & (y >= w1[1]) & (y <= w1[3])).sum()
    assert 0.5 * 150_000 <= n <= 150_000
    # the window must sit on the dense blob (~150, 150), not the sparse one
    assert abs((w1[0] + w1[2]) / 2 - 150) < 20
    assert abs((w1[1] + w1[3]) / 2 - 150) < 20


def test_densest_window_tiny_window_still_under_cap():
    # pathological: very tight blob, window must shrink below the global side
    rng = np.random.default_rng(1)
    x = rng.normal(100, 0.5, 500_000)
    y = rng.normal(100, 0.5, 500_000)
    w = U.densest_window(x, y, target=150_000, cap=150_000, seed=3)
    n = ((x >= w[0]) & (x <= w[2]) & (y >= w[1]) & (y <= w[3])).sum()
    assert n <= 150_000
    assert n >= 50_000


def test_px_to_bbox_and_crop_px_are_half_open():
    assert O.px_to_bbox((1.2, 3.7, 10.01, 20.99)) == (1, 3, 11, 21)
    df = pd.DataFrame({"x": [0.0, 1.0, 9.9, 10.0, 11.0],
                       "y": [0.0, 5.0, 5.0, 5.0, 5.0]})
    out = O.crop_px(df, (1, 0, 10, 10), "x", "y")
    # keeps [1, 10): rows with x=1.0 and x=9.9
    assert list(out["x"]) == [1.0, 9.9]


def test_image_slice_bbox_clips():
    assert O.image_slice_bbox((0, 0, 10, 10), (5, 5)) == (0, 5, 0, 5)
    assert O.image_slice_bbox((-3, -4, 100, 100), (50, 40)) == (0, 50, 0, 40)


def test_contract_classes_boundaries():
    assert U.gene_panel_class(49) == "tiny"
    assert U.gene_panel_class(50) == "small"
    assert U.gene_panel_class(249) == "small"
    assert U.gene_panel_class(250) == "medium"
    assert U.gene_panel_class(699) == "medium"
    assert U.gene_panel_class(700) == "large"
    assert U.gene_panel_class(1999) == "large"
    assert U.gene_panel_class(2000) == "huge"
    assert U.cell_density_class(2499.9) == "sparse"
    assert U.cell_density_class(2500.0) == "medium"
    assert U.cell_density_class(7000.0) == "medium"
    assert U.cell_density_class(7000.1) == "dense"


def test_rasterize_polygons_and_lookup():
    polys = {1: [[0, 0], [10, 0], [10, 10], [0, 10]],
             7: [[5, 5], [20, 5], [20, 20], [5, 20]]}
    img = U.rasterize_polygons(polys, (0, 0, 25, 25))
    assert img.shape == (25, 25)
    assert img[2, 2] == 1 and img[15, 15] == 7 and img[0, 24] == 0
    labels = U.lookup_labels(polys, np.array([2.0, 17.0, 100.0]),
                             np.array([2.0, 17.0, 100.0]))
    assert list(labels) == [1, 7, 0]


def test_cc_label_binary():
    mask = np.zeros((20, 20), dtype=np.uint8)
    mask[2:5, 2:5] = 255
    mask[10:15, 10:15] = 255
    lab = U.cc_label_binary(mask)
    assert lab.dtype == np.uint16
    assert set(np.unique(lab)) == {0, 1, 2}


def test_sort_molecules_contract_order():
    df = pd.DataFrame({"y": [2, 1, 1], "x": [1, 3, 2], "gene": ["a", "b", "c"]})
    out = U.sort_molecules(df)
    assert list(out["x"]) == [2, 3, 1]


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def test_manifest_is_valid():
    m = O.load_manifest()
    ids = [d["id"] for d in m["datasets"]]
    assert len(ids) == len(set(ids))
    assert m["group"] == "real_other"
    tiers = set()
    for d in m["datasets"]:
        assert d["builder"] in O.BUILDERS, d["id"]
        assert d["tier"] in ("quick", "full")
        assert "source" in d and "url" in d["source"]
        assert "baysor" in d
        b = d["baysor"]
        assert (U.repo_root() / b["config"]).exists(), b["config"]
        assert isinstance(b["scale_um"], (int, float))
        assert b["scale_std"] is not None
        assert b["min_molecules_per_cell"] > 0
        cap = d["crop"]["cap"]
        assert cap <= (150_000 if d["tier"] == "quick" else 3_000_000)
        tiers.add(d["tier"])
    # at least one quick crop per builder, plus full crops for the big sources
    builders = {d["builder"] for d in m["datasets"]}
    assert {"iss", "osmfish", "starmap", "ileum", "cosmx_nsclc", "cosmx_wtx"} <= builders
    assert tiers == {"quick", "full"}
    # WTX: default mrf/ICA does not scale to ~19k genes -> louvain + recorded
    # default-timeout run (perf-stress: default ICA init)
    wtx = [d for d in m["datasets"] if d["builder"] == "cosmx_wtx"]
    for d in wtx:
        assert d["baysor"]["extra_args"] == ["--cluster-method", "louvain"], d["id"]
    quick_wtx = [d for d in wtx if d["tier"] == "quick"]
    assert len(quick_wtx) == 1
    assert quick_wtx[0].get("smoke_record_default_timeout") == 1200


# ---------------------------------------------------------------------------
# meta.json / smoke plumbing
# ---------------------------------------------------------------------------

def _dummy_spec():
    return {
        "id": "dummy_quick",
        "tier": "quick",
        "platform": "Test",
        "tissue": "test tissue",
        "seed": 1,
        "source": {"url": "http://example.com", "doi": None,
                   "license": "none", "original_dataset": "x"},
        "baysor": {"scale_um": 5.0, "scale_std": "25%",
                   "min_molecules_per_cell": 20, "config": "configs/example_config.toml",
                   "extra_args": [], "prior_confidence": 0.2},
    }


def test_build_meta_contract_shape():
    meta = O.build_meta(
        _dummy_spec(),
        bbox_um=[0, 0, 100, 200],
        z_range_um=None,
        n_molecules=1000,
        n_genes=42,
        n_cells=800,
        crop_note="note",
        images=[],
        prior="none",
        source_extra={"pixel_size_um": 0.5},
        difficulty_notes="dn",
    )
    assert set(meta) == {"id", "kind", "tier", "platform", "source", "crop",
                         "stats", "difficulty", "baysor", "images", "truth"}
    assert meta["kind"] == "real" and meta["truth"] is None
    assert meta["stats"]["area_um2"] == 20000.0
    assert meta["stats"]["molecules_per_um2"] == 0.05
    # 800 cells / 0.02 mm^2 = 40000 cells/mm^2 -> dense
    assert meta["difficulty"]["cell_density"] == "dense"
    assert meta["difficulty"]["gene_panel"] == "tiny"
    # prior=none -> prior_confidence null
    assert meta["baysor"]["prior"] == "none"
    assert meta["baysor"]["prior_confidence"] is None


def test_build_meta_prior_confidence_passthrough():
    spec = _dummy_spec()
    meta = O.build_meta(spec, bbox_um=[0, 0, 10, 10], z_range_um=None,
                        n_molecules=10, n_genes=5, n_cells=1, crop_note="",
                        images=[], prior="column", source_extra={},
                        difficulty_notes="")
    assert meta["baysor"]["prior_confidence"] == 0.2


def test_parse_time_v_output():
    import tempfile

    text = """\tCommand being timed: "baysor run"
\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:07.42
\tMaximum resident set size (kbytes): 2048000
"""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(text)
        p = Path(fh.name)
    wall, rss = O.parse_time_v(p)
    assert rss == 2048000
    assert abs(wall - 67.42) < 0.01


def test_baysor_command_construction(tmp_path):
    spec = _dummy_spec()
    ds = tmp_path / "ds"
    (ds / "images").mkdir(parents=True)
    meta = O.build_meta(spec, bbox_um=[0, 0, 1, 1], z_range_um=None, n_molecules=1,
                        n_genes=1, n_cells=1, crop_note="", images=[],
                        prior="image:images/mask.tif", source_extra={},
                        difficulty_notes="")
    meta["bayso r"] = None  # typo key must not break anything
    del meta["bayso r"]
    (ds / "meta.json").write_text(json.dumps(meta))
    (ds / "molecules.parquet").write_bytes(b"")
    cmd = O.build_baysor_command(spec, ds, Path("/bin/baysor"),
                                 tmp_path / "out", tmp_path / "time.txt")
    assert "run" in cmd and "-s" in cmd and "5.0" in cmd
    assert cmd[-1].endswith("images/mask.tif")
    assert (ds / "molecules.parquet") in [Path(c) for c in cmd]

    meta["baysor"]["prior"] = "column"
    (ds / "meta.json").write_text(json.dumps(meta))
    cmd = O.build_baysor_command(spec, ds, Path("/bin/baysor"),
                                 tmp_path / "out", tmp_path / "time.txt")
    assert cmd[-1] == ":prior"


def test_ensure_label_volume_detects_labels_and_binary():
    vol = np.arange(0, 1000, dtype=np.int32).reshape(10, 10, 10) % 60
    assert O._ensure_label_volume(vol) is vol
    binary = np.zeros((10, 10), dtype=np.uint8)
    binary[1:4, 1:4] = 255
    out = O._ensure_label_volume(binary)
    assert out.max() == 1 and out[2, 2] == 1
    with pytest.raises(RuntimeError):
        O._ensure_label_volume(np.array([[0, 3], [5, 0]]))  # 3 values: suspicious


def test_as_zyx_layout_normalisation():
    zyx = np.zeros((9, 50, 50))
    assert O._as_zyx(zyx) is zyx
    yxz = np.zeros((5000, 5000, 9))
    out = O._as_zyx(yxz)
    assert out.shape == (9, 5000, 5000)
