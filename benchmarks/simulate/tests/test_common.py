"""Unit tests for benchmarks/simulate/common.py."""
import numpy as np
import pandas as pd
import pytest

import common


def test_hex_centers_bounds_and_density():
    rng = np.random.default_rng(0)
    extent, spacing = 300.0, 15.0
    c = common.hex_centers(extent, spacing, 0.1, rng)
    assert c.shape[1] == 2
    assert (c >= 0).all() and (c <= extent + 1e-9).all()
    expected = extent ** 2 / (spacing ** 2 * np.sqrt(3) / 2)
    assert abs(len(c) - expected) / expected < 0.1

    rng = np.random.default_rng(0)
    c2 = common.hex_centers(extent, spacing, 0.1, rng, min_edge_dist_um=5.0)
    assert (c2 >= 5.0).all() and (c2 <= extent - 5.0).all()


def test_hex_centers_deterministic():
    a = common.hex_centers(200.0, 14.0, 0.1, np.random.default_rng(7))
    b = common.hex_centers(200.0, 14.0, 0.1, np.random.default_rng(7))
    assert np.array_equal(a, b)


def test_delaunay_adjacency_symmetric_and_connected():
    rng = np.random.default_rng(1)
    c = common.hex_centers(200.0, 15.0, 0.1, rng)
    adj = common.delaunay_adjacency(c)
    assert len(adj) == len(c)
    for i, nbrs in enumerate(adj):
        assert i not in nbrs
        for j in nbrs:
            assert i in adj[j]


def test_dsatur_proper_coloring_and_determinism():
    rng = np.random.default_rng(2)
    c = common.hex_centers(250.0, 15.0, 0.15, rng)
    adj = common.delaunay_adjacency(c)
    colors = common.dsatur_color(adj)
    assert (colors >= 0).all()
    for i, nbrs in enumerate(adj):
        for j in nbrs:
            assert colors[i] != colors[j], f"adjacent {i}-{j} share a colour"
    again = common.dsatur_color(adj)
    assert np.array_equal(colors, again)


def test_dsatur_small_graphs():
    # triangle needs 3 colours
    tri = [{1, 2}, {0, 2}, {0, 1}]
    assert len(set(common.dsatur_color(tri))) == 3
    # path of 4 needs 2
    path = [{1}, {0, 2}, {1, 3}, {2}]
    col = common.dsatur_color(path)
    assert len(set(col)) == 2
    assert col[0] != col[1] and col[2] != col[3]
    # empty graph
    assert common.dsatur_color([set(), set()]).tolist() == [0, 0]


def test_expression_profiles_normalised_and_markers():
    p = common.expression_profiles(4, 100, np.random.default_rng(3),
                                   markers_per_type=8, marker_mass=0.55)
    assert p.shape == (4, 100)
    assert np.allclose(p.sum(axis=1), 1.0)
    assert (p > 0).all()


def test_expression_profiles_disjoint_blocks():
    p = common.expression_profiles(5, 100, np.random.default_rng(4),
                                   markers_per_type=8, disjoint=True)
    assert np.allclose(p.sum(axis=1), 1.0)
    blocks = np.array_split(np.arange(100), 5)
    for t, block in enumerate(blocks):
        others = np.setdiff1d(np.arange(100), block)
        assert p[t, others].sum() == 0.0     # strictly disjoint gene sets
        assert abs(p[t, block].sum() - 1.0) < 1e-12
    # no gene used by two types
    used = (p > 0)
    assert (used.sum(axis=0) <= 1).all()


def test_sample_counts_bounds():
    n = common.sample_counts(1000, np.random.default_rng(5))
    assert n.min() >= 50 and n.max() <= 300
    assert 100 < n.mean() < 220


def test_nucleus_prior_synthetic():
    centers = np.array([[0.0, 0.0], [10.0, 0.0]])
    xy = np.array([[1.0, 0.0],     # inside cell 0
                   [9.0, 0.0],     # inside cell 1
                   [5.0, 0.0],     # between (3 um nucleus)
                   [2.5, 0.0]])    # inside cell 0
    label, idx = common.nucleus_prior(xy, centers, 3.0)
    assert label.tolist() == [1, 2, 0, 1]
    assert idx.tolist() == [0, 1, 0, 0]
    # 3D distances
    xyz = np.array([[0.0, 0.0, 2.5]])
    label3, _ = common.nucleus_prior(xyz, np.array([[0.0, 0.0, 0.0]]), 3.0)
    assert label3.tolist() == [1]


def test_class_boundaries():
    assert common.gene_panel_class(49) == "tiny"
    assert common.gene_panel_class(50) == "small"
    assert common.gene_panel_class(249) == "small"
    assert common.gene_panel_class(250) == "medium"
    assert common.gene_panel_class(700) == "large"
    assert common.gene_panel_class(2000) == "large"
    assert common.gene_panel_class(2001) == "huge"
    assert common.density_class(2499) == "sparse"
    assert common.density_class(2500) == "medium"
    assert common.density_class(7000) == "medium"
    assert common.density_class(7001) == "dense"


def test_make_molecules_sorts_and_casts():
    df = pd.DataFrame({
        "x": [2.0, 1.0, 1.0], "y": [1.0, 2.0, 1.0],
        "gene": ["g0", "g1", "g0"], "prior": [0, 1, 0],
        "cell": [1, 0, 2], "interior": [True, False, True],
        "celltype": ["a", "", "b"],
    })
    out = common.make_molecules(df)
    assert list(out["x"]) == [1.0, 2.0, 1.0]      # sorted by (y, x)
    assert list(out["y"]) == [1.0, 1.0, 2.0]
    assert out["prior"].dtype == np.int32
    assert out["cell"].dtype == np.int32
    assert out["interior"].dtype == bool
    assert out["x"].dtype == np.float64


def test_write_dataset_roundtrip(tmp_path):
    df = pd.DataFrame({
        "x": [1.0], "y": [2.0], "gene": ["g0"], "prior": [0],
        "cell": [1], "interior": [True], "celltype": ["ct"],
    })
    df = common.make_molecules(df)
    meta = common.make_meta(
        id="t", tier="quick", source={}, crop={}, 
        stats=common.make_stats(df, area_um2=100.0),
        difficulty={"cell_density": "sparse", "gene_panel": "tiny", "notes": ""},
        baysor=common.make_baysor(5.0), truth=None)
    hashes = common.write_dataset(tmp_path, df, meta)
    assert set(hashes) == {"molecules.parquet", "meta.json"}
    back = pd.read_parquet(tmp_path / "molecules.parquet")
    assert len(back) == 1
    import json
    assert json.loads((tmp_path / "meta.json").read_text())["id"] == "t"
    assert common.sha256(tmp_path / "molecules.parquet") == hashes["molecules.parquet"]


def test_make_meta_rejects_oversized_tier():
    df = pd.DataFrame({"cell": [1]})
    stats = {"n_molecules": 200_000, "n_genes": 1, "area_um2": 1.0,
             "molecules_per_um2": 1.0, "n_vendor_cells": 0,
             "vendor_cells_per_mm2": 0.0, "n_true_cells": 1,
             "true_cells_per_mm2": 1.0}
    with pytest.raises(ValueError, match="exceeds quick budget"):
        common.make_meta(id="t", tier="quick", source={}, crop={}, stats=stats,
                         difficulty={}, baysor={}, truth=None)


def test_child_rngs_independent_streams():
    a = common.child_rngs(42, ["x", "y"])
    b = common.child_rngs(42, ["x", "y"])
    assert np.array_equal(a["x"].random(10), b["x"].random(10))
    # streams are independent of each other
    c = common.child_rngs(42, ["x", "y"])
    x = c["x"].random(5)
    y = c["y"].random(5)
    assert not np.array_equal(x, y)
