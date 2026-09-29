"""Integration tests for the trivial scenarios: contract, exact ground truth,
geometry invariants and determinism.

The tests reconstruct the generator's internal state (cell centres, type
assignment) from ``meta.truth.params`` + the seed, so they independently
verify that the stored molecule table really matches the documented geometry.
"""
import numpy as np
import pandas as pd
import pytest
from scipy.spatial import cKDTree

import trivial

RNG_STREAMS = trivial._RNG_STREAMS

# tiny parameter sets: fast, but exercising every scenario
TINY = {
    "circles_gaps": dict(extent_um=160.0, spacing_um=16.0, radius_um=6.1),
    "tiled_distinct": dict(extent_um=160.0, spacing_um=14.5),
    "tiled_same": dict(extent_um=160.0, spacing_um=14.5),
    "mixed_sizes": dict(extent_um=180.0, spacing_um=20.5),
    "sparse_noisy": dict(extent_um=220.0, spacing_um=40.0),
    "circles_gaps_3d": dict(extent_um=160.0, spacing_um=16.0, radius_um=6.1),
    "elongated_gaps": dict(extent_um=200.0, spacing_um=27.0),
}


def build(scenario: str, seed: int = 1201, **over):
    kwargs = dict(TINY[scenario])
    kwargs.update(over)
    return trivial.SCENARIOS[scenario](seed=seed, dataset_id=f"t_{scenario}",
                                       **kwargs)


# ---------------------------------------------------------------------------
# generic contract checks
# ---------------------------------------------------------------------------

def check_contract(df: pd.DataFrame, meta: dict, tier_budget: int = 150_000):
    cols = list(df.columns)
    assert cols[:2] == ["x", "y"]
    assert {"gene", "prior", "cell", "interior", "celltype"} <= set(cols)
    assert df["x"].dtype == np.float64 and df["y"].dtype == np.float64
    assert df["prior"].dtype == np.int32 and df["cell"].dtype == np.int32
    assert df["interior"].dtype == bool
    assert df["gene"].dtype == object
    # sorted by (y, x)
    assert df[["y", "x"]].equals(df[["y", "x"]].sort_values(["y", "x"],
                                                            kind="stable").reset_index(drop=True))
    # meta contract
    assert meta["kind"] == "sim" and meta["tier"] in ("quick", "full")
    for key in ("id", "source", "crop", "stats", "difficulty", "baysor", "truth"):
        assert key in meta
    b = meta["baysor"]
    assert b["prior"] == "column" and b["scale_um"] > 0
    assert b["scale_std"] == "25%"
    t = meta["truth"]
    assert t["seed"] == meta["source"]["seed"]
    assert t["params"]["n_genes"] if "n_genes" in t["params"] else True
    assert len(df) <= tier_budget
    assert meta["stats"]["n_molecules"] == len(df)


def reconstruct_centers(scenario: str, meta: dict):
    """Replay the geometry RNG stream to recover the cell centres."""
    p = meta["truth"]["params"]
    rngs = trivial.common.child_rngs(meta["truth"]["seed"], RNG_STREAMS)
    if scenario in ("circles_gaps", "circles_gaps_3d", "mixed_sizes"):
        centers = trivial._lattice_with_gaps(
            p["extent_um"], p["spacing_um"], p["jitter_frac"], rngs["geo"],
            min_edge_dist_um=p.get("radius_um", p.get("radius_large_um")),
            min_gap_um=p["min_gap_um"],
            max_radius_um=p.get("radius_um", p.get("radius_large_um")))
    elif scenario == "elongated_gaps":
        a_max = max(p["semi_major_fibro_um"], p["semi_major_neuron_um"])
        centers = trivial._lattice_with_gaps(
            p["extent_um"], p["spacing_um"], p["jitter_frac"], rngs["geo"],
            min_edge_dist_um=a_max, min_gap_um=p["min_gap_um"],
            max_radius_um=a_max)
    else:
        centers = trivial.common.hex_centers(
            p["extent_um"], p["spacing_um"], p["jitter_frac"], rngs["geo"])
    return centers, rngs, p


@pytest.mark.parametrize("scenario", sorted(TINY))
def test_scenario_contract(scenario):
    df, meta = build(scenario)
    check_contract(df, meta)
    # ground-truth cell ids: background is 0, cells are 1..N contiguous
    cells = np.unique(df["cell"])
    assert cells.min() == 0                       # every scenario has background
    assert cells.max() == meta["stats"]["n_true_cells"]
    assert (cells[1:] == np.arange(1, cells.max() + 1)).all()
    # background molecules carry empty celltype and no interior exclusion rule
    bg = df["cell"] == 0
    assert (df.loc[bg, "celltype"] == "").all()
    assert (df.loc[bg, "prior"] == 0).sum() <= bg.sum()
    # at least some molecules get a prior
    assert 0.01 < (df["prior"] > 0).mean() < 0.8


@pytest.mark.parametrize("scenario", sorted(TINY))
def test_scenario_deterministic(scenario):
    df1, meta1 = build(scenario, seed=7)
    df2, meta2 = build(scenario, seed=7)
    pd.testing.assert_frame_equal(df1, df2)
    assert meta1 == meta2


# ---------------------------------------------------------------------------
# geometry invariants (reconstructed state)
# ---------------------------------------------------------------------------

def test_circles_gaps_geometry_and_prior():
    df, meta = build("circles_gaps")
    centers, rngs, p = reconstruct_centers("circles_gaps", meta)
    r, r_nuc = p["radius_um"], p["r_nucleus_um"]

    # clear gaps between cells
    d, _ = cKDTree(centers).query(centers, k=2)
    assert d[:, 1].min() >= 2 * r + p["min_gap_um"] - 1e-6

    cell = df["cell"].to_numpy()
    mol = df[["x", "y"]].to_numpy()
    inside = cell > 0
    dist_own = np.linalg.norm(mol[inside] - centers[cell[inside] - 1], axis=1)
    assert dist_own.max() <= r + 1e-9             # molecules inside their disc

    # background: uniform over the domain
    bg = ~inside
    assert abs(bg.mean() - p["bg_frac"]) < 0.01

    # prior rule: nearest centre within r_nucleus
    prior, _ = trivial.common.nucleus_prior(mol, centers, r_nuc)
    assert np.array_equal(prior, df["prior"].to_numpy())

    # interior rule: source cell far enough from the domain edge
    margin = r + 5.0
    c = centers
    interior_cell = ((c[:, 0] >= margin) & (c[:, 0] <= p["extent_um"] - margin)
                     & (c[:, 1] >= margin) & (c[:, 1] <= p["extent_um"] - margin))
    expected = np.zeros(len(df), bool)
    expected[inside] = interior_cell[cell[inside] - 1]
    bxy = mol[bg]
    expected[bg] = ((bxy[:, 0] >= margin) & (bxy[:, 0] <= p["extent_um"] - margin)
                    & (bxy[:, 1] >= margin) & (bxy[:, 1] <= p["extent_um"] - margin))
    assert np.array_equal(expected, df["interior"].to_numpy())


def test_tiled_distinct_neighbours_have_different_types():
    df, meta = build("tiled_distinct")
    centers, rngs, p = reconstruct_centers("tiled_distinct", meta)
    # one celltype per cell id, contiguous cell ids matching centre order
    per_cell = df[df.cell > 0].groupby("cell")["celltype"].nunique()
    assert (per_cell == 1).all()
    n_cells = meta["stats"]["n_true_cells"]
    assert n_cells == len(centers)
    ctype = np.empty(n_cells, object)
    sub = df[df.cell > 0]
    ctype[sub["cell"].to_numpy() - 1] = sub["celltype"].to_numpy()

    # Voronoi neighbours = Delaunay edges; all must have different types
    adj = trivial.common.delaunay_adjacency(centers)
    for i, nbrs in enumerate(adj):
        for j in nbrs:
            assert ctype[i] != ctype[j], f"neighbours {i}-{j} share a type"

    # disjoint gene sets: no gene observed for two different cell types
    genes_per_type = (sub.groupby("celltype")["gene"].unique())
    types = list(genes_per_type.index)
    for a in range(len(types)):
        for b in range(a + 1, len(types)):
            inter = set(genes_per_type[types[a]]) & set(genes_per_type[types[b]])
            assert not inter, f"gene sets of {types[a]}/{types[b]} overlap: {inter}"

    # Voronoi membership: every molecule's own centre is its nearest centre
    mol = sub[["x", "y"]].to_numpy()
    own = centers[sub["cell"].to_numpy() - 1]
    d_own = np.linalg.norm(mol - own, axis=1)
    d_all, _ = cKDTree(centers).query(mol, k=2)
    assert (d_own <= d_all[:, 0] + 1e-6).all()


def test_tiled_same_single_type():
    df, meta = build("tiled_same")
    sub = df[df.cell > 0]
    assert sub["celltype"].nunique() == 1
    # same geometry as tiled_distinct with the same seed
    dfd, _ = build("tiled_distinct")
    pd.testing.assert_series_equal(df["x"], dfd["x"])
    pd.testing.assert_series_equal(df["y"], dfd["y"])
    pd.testing.assert_series_equal(df["cell"], dfd["cell"])


def test_mixed_sizes_radii_and_composition():
    df, meta = build("mixed_sizes")
    centers, rngs, p = reconstruct_centers("mixed_sizes", meta)
    n_cells = len(centers)
    is_small = rngs["types"].random(n_cells) < p["p_small"]
    radius = np.where(is_small, p["radius_small_um"], p["radius_large_um"])

    # gaps: no pair closer than the two radii + min gap
    tree = cKDTree(centers)
    pairs = tree.query_pairs(2 * p["radius_large_um"] + p["min_gap_um"],
                             output_type="ndarray")
    if len(pairs):
        d = np.linalg.norm(centers[pairs[:, 0]] - centers[pairs[:, 1]], axis=1)
        need = radius[pairs[:, 0]] + radius[pairs[:, 1]] + p["min_gap_um"]
        assert (d >= need - 1e-6).all()

    sub = df[df.cell > 0]
    mol = sub[["x", "y"]].to_numpy()
    cid = sub["cell"].to_numpy() - 1
    d = np.linalg.norm(mol - centers[cid], axis=1)
    assert (d <= radius[cid] + 1e-9).all()        # inside own disc
    ct = sub["celltype"].to_numpy()
    assert (ct == np.where(is_small[cid], "immune", "tumour")).all()
    # large cells carry more molecules than small ones
    per_cell = sub.groupby("cell")["cell"].size().to_numpy()
    assert per_cell[~is_small].mean() > per_cell[is_small].mean()


def test_sparse_noisy_background_fraction():
    df, meta = build("sparse_noisy")
    p = meta["truth"]["params"]
    frac = (df["cell"] == 0).mean()
    assert 0.20 <= frac <= 0.30
    # sparse: density class must be sparse
    assert meta["difficulty"]["cell_density"] == "sparse"


def test_circles_gaps_3d_slab_and_prior():
    df, meta = build("circles_gaps_3d")
    centers, rngs, p = reconstruct_centers("circles_gaps_3d", meta)
    slab, r, r_nuc = p["slab_um"], p["radius_um"], p["r_nucleus_um"]
    # centre z drawn after the lattice on the same stream
    cz = rngs["geo"].uniform(0.0, slab, size=len(centers))
    centers3 = np.column_stack([centers, cz])
    assert "z" in df.columns and df["z"].dtype == np.float64
    assert df["z"].min() >= 0.0 and df["z"].max() <= slab

    cell = df["cell"].to_numpy()
    mol = df[["x", "y", "z"]].to_numpy()
    inside = cell > 0
    d = np.linalg.norm(mol[inside] - centers3[cell[inside] - 1], axis=1)
    assert d.max() <= r + 1e-9                     # inside the (cut) sphere
    prior, _ = trivial.common.nucleus_prior(mol, centers3, r_nuc)
    assert np.array_equal(prior, df["prior"].to_numpy())


# ---------------------------------------------------------------------------
# gene-panel ablation
# ---------------------------------------------------------------------------

def test_panel_variants_share_geometry():
    df100, m100 = build("circles_gaps", n_genes=100)
    df5000, m5000 = build("circles_gaps", n_genes=5000)
    pd.testing.assert_series_equal(df100["x"], df5000["x"])
    pd.testing.assert_series_equal(df100["y"], df5000["y"])
    pd.testing.assert_series_equal(df100["cell"], df5000["cell"])
    pd.testing.assert_series_equal(df100["prior"], df5000["prior"])
    assert df100["gene"].nunique() == 100
    assert 2000 < df5000["gene"].nunique() <= 5000
    assert m100["difficulty"]["gene_panel"] == "small"
    assert m5000["difficulty"]["gene_panel"] == "huge"


def test_expression_profiles_few_markers_per_type():
    """Molecule counts concentrate on a few markers plus shared background."""
    df, meta = build("circles_gaps", n_genes=100)
    sub = df[df.cell > 0]
    # per type, a handful of genes carry most of the mass
    for ct, g in sub.groupby("celltype")["gene"]:
        counts = g.value_counts().to_numpy()
        top = counts[:8].sum() / counts.sum()
        assert 0.35 < top < 0.85, f"{ct}: markers carry {top:.2f}"


def test_tier_budget_enforced():
    with pytest.raises(ValueError, match="exceeds quick budget"):
        build("circles_gaps", tier="quick", extent_um=1400.0)


# ---------------------------------------------------------------------------
# elongated (irregular-shape) cells and imperfect-prior variants
# ---------------------------------------------------------------------------

def test_elongated_geometry_gaps_and_ellipse_containment():
    df, meta = build("elongated_gaps")
    centers, rngs, p = reconstruct_centers("elongated_gaps", meta)
    a_max = max(p["semi_major_fibro_um"], p["semi_major_neuron_um"])

    # clear gaps: conservative circular envelopes never overlap
    d, _ = cKDTree(centers).query(centers, k=2)
    assert d[:, 1].min() >= 2 * a_max + p["min_gap_um"] - 1e-6

    sub = df[df.cell > 0]
    cid = sub["cell"].to_numpy() - 1
    ct = sub["celltype"].to_numpy()
    assert set(ct) <= {"fibro", "neuron"}
    a = np.where(ct == "fibro", p["semi_major_fibro_um"],
                 p["semi_major_neuron_um"])
    b = np.where(ct == "fibro", p["semi_minor_fibro_um"],
                 p["semi_minor_neuron_um"])

    # replay the per-cell orientation (types stream: shape draw, then theta;
    # proportions are fixed, no dirichlet draw in this scenario)
    rngs2 = trivial.common.child_rngs(meta["truth"]["seed"], RNG_STREAMS)
    n_cells = len(centers)
    rngs2["types"].random(n_cells)          # p_fibro draw
    theta = rngs2["types"].uniform(0.0, np.pi, size=n_cells)
    th = theta[cid]
    d = sub[["x", "y"]].to_numpy() - centers[cid]
    u = d[:, 0] * np.cos(th) + d[:, 1] * np.sin(th)
    v = -d[:, 0] * np.sin(th) + d[:, 1] * np.cos(th)
    inside = (u / a) ** 2 + (v / b) ** 2 <= 1.0 + 1e-9
    assert inside.all(), "molecules outside their ellipse"
    # genuinely elongated
    assert p["semi_major_fibro_um"] / p["semi_minor_fibro_um"] >= 2.5
    assert p["semi_major_neuron_um"] / p["semi_minor_neuron_um"] >= 4.0


def test_imprior_keeps_truth_and_degrades_prior():
    opts = {"kind": "imperfect", "seed": 7701, "base": "t_base",
            "miss_frac": 0.2, "shift_um": [1.0, 2.0], "merge_frac": 0.05}
    df0, meta0 = build("circles_gaps", prior_opts=None)
    df1, meta1 = build("circles_gaps", prior_opts=opts)
    # truth columns byte-identical
    for col in ("x", "y", "gene", "cell", "interior", "celltype"):
        assert (df0[col].to_numpy() == df1[col].to_numpy()).all(), col
    # prior degraded but still only nucleus labels or 0
    info = meta1["truth"]["prior"]
    n_cells = meta1["stats"]["n_true_cells"]
    assert info["n_missed"] == round(0.2 * n_cells)
    assert info["n_merge_sources"] == round(0.05 * n_cells)
    assert set(np.unique(df1["prior"])) <= set(range(0, n_cells + 1))
    assert (df1["prior"] > 0).mean() < (df0["prior"] > 0).mean()
    # perfect-prior dataset records no such block (schema stability)
    assert "prior" not in meta0["truth"]
    # determinism
    df2, meta2 = build("circles_gaps", prior_opts=opts)
    pd.testing.assert_frame_equal(df1, df2)
    assert meta1 == meta2
