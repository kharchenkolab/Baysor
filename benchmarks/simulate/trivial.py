"""Trivial simulated datasets with exact ground truth.

Six scenarios (see ``benchmarks/datasets/sim.yaml`` for the instances):

1. ``circles_gaps``    round cells of equal radius on a jittered hex lattice
                       with clear gaps; several cell types, ~1% uniform
                       background (``cell = 0``).  Easy: Baysor should be near
                       perfect.
2. ``tiled_distinct``  Voronoi tiling of a jittered hex lattice, no gaps; the
                       adjacency graph is coloured so neighbouring cells have
                       different types and the types have *disjoint* gene
                       sets, so boundaries are resolvable by composition only.
3. ``tiled_same``      the same geometry with every cell of one identical
                       type: boundaries are unresolvable by composition (hard
                       control).
4. ``mixed_sizes``     small (immune-like, r ~ 3 um) and large (tumour-like,
                       r ~ 9 um) round cells with small gaps and different
                       compositions.
5. ``sparse_noisy``    sparse cells in a large area with 20-30% background
                       noise.
6. ``circles_gaps_3d`` 3D spheres in a thin slab, ``z`` column present, cells
                       cut by the section.

Every scenario is a function of explicit parameters plus a seed and returns
``(molecules_df, meta)`` following ``benchmarks/README.md``.  Molecule counts
are lognormal in ``[50, 300]`` (per-cell medians vary by scenario); genes are
multinomial draws from per-type profiles with a few markers per type over a
shared low-level lognormal background.  The ``prior`` column labels molecules
within ``r_nucleus`` of the true centre (0 = none).

CLI (one dataset)::

    python benchmarks/simulate/trivial.py circles_gaps \\
        --id sim_circles_gaps_g100 --seed 1201 --n-genes 100 \\
        -o $BAYSOR_BENCH_DATA/sim/sim_circles_gaps_g100
"""
from __future__ import annotations

import argparse
import inspect
import json
import os
import sys

import numpy as np
import pandas as pd

try:
    from . import common
except ImportError:  # running as a plain script
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import common

GENERATOR_VERSION = 1

# RNG stream names shared by all scenarios (see common.child_rngs).
_RNG_STREAMS = ["geo", "types", "profiles", "counts", "pos", "genes", "bg"]


# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------

def _enforce_min_separation(centers: np.ndarray, threshold: float, lo: float,
                            hi: float, max_iter: int = 400) -> np.ndarray | None:
    """Iteratively push centres apart until every pair is ``>= threshold``.

    Violating pairs are separated halfway along their axis (averaged over all
    pairs a point takes part in, relaxed), then re-clipped to ``[lo, hi]``.
    Converges when the largest deficit drops below 1e-9 um (the dynamics
    approach the threshold from below, so exact emptiness is unreachable);
    returns ``None`` on stall, and the caller then falls back to a weaker
    jitter.
    """
    from scipy.spatial import cKDTree

    centers = np.asarray(centers, dtype=np.float64).copy()
    best_state, best_viol, best_it = None, np.inf, -1
    for it in range(max_iter):
        pairs = cKDTree(centers).query_pairs(threshold, output_type="ndarray")
        if len(pairs) == 0:
            return centers
        delta = centers[pairs[:, 1]] - centers[pairs[:, 0]]
        d = np.linalg.norm(delta, axis=1)
        deficit = threshold - d
        viol = float(deficit.max())
        if viol < 1e-9:
            return centers
        if viol < best_viol - 1e-12:
            best_viol, best_it, best_state = viol, it, centers.copy()
        elif it - best_it > 60:  # stalled
            break
        unit = delta / np.where(d > 0, d, 1.0)[:, None]
        move = (0.4 * deficit)[:, None] * unit
        disp = np.zeros_like(centers)
        cnt = np.zeros(len(centers))
        # push each end away from the other along the pair axis
        for a, b, m in ((pairs[:, 0], pairs[:, 1], -move),
                        (pairs[:, 1], pairs[:, 0], move)):
            np.add.at(disp, a, m)
            np.add.at(cnt, a, 1.0)
        cnt = np.maximum(cnt, 1.0)
        centers += disp / cnt[:, None]
        centers = np.clip(centers, lo, hi)
    # stalled: accept only a practically converged state
    if best_state is not None and best_viol < 1e-6:
        return best_state
    return None


def _lattice_with_gaps(extent_um: float, spacing_um: float, jitter_frac: float,
                       rng: np.random.Generator, min_edge_dist_um: float,
                       min_gap_um: float, max_radius_um: float) -> np.ndarray:
    """Jittered hex lattice whose discs of ``max_radius_um`` keep a gap.

    Draws a jittered lattice, then repairs violating pairs (see
    ``_enforce_min_separation``) so every centre pair is at least
    ``2 * radius + min_gap`` apart.  If repair stalls the jitter is halved and
    a fresh lattice is drawn; the regular lattice always satisfies the
    constraint, so this terminates.
    """
    threshold = 2.0 * max_radius_um + min_gap_um
    lo = min_edge_dist_um
    hi = extent_um - min_edge_dist_um
    scale, tries = 1.0, 0
    while True:
        centers = common.hex_centers(extent_um, spacing_um, jitter_frac * scale,
                                     rng, min_edge_dist_um=min_edge_dist_um)
        fixed = _enforce_min_separation(centers, threshold, lo, hi)
        if fixed is not None:
            return fixed
        tries += 1
        if tries > 50:
            raise RuntimeError("could not create a gap-preserving lattice")
        if tries % 5 == 0:
            scale *= 0.5


def _gene_names(n_genes: int) -> list[str]:
    w = len(str(n_genes - 1))
    return [f"g{i:0{w}d}" for i in range(n_genes)]


def _emit_genes(mol_cell: np.ndarray, type_of_cell: np.ndarray,
                profiles: np.ndarray, n_genes: int,
                rng: np.random.Generator) -> np.ndarray:
    """Per-cell multinomial gene draws expanded to one label per molecule.

    ``mol_cell`` must be sorted ascending and grouped by cell.
    """
    names = np.asarray(_gene_names(n_genes), dtype=object)
    genes = np.empty(len(mol_cell), dtype=object)
    starts = np.searchsorted(mol_cell, np.arange(type_of_cell.size))
    ends = np.searchsorted(mol_cell, np.arange(type_of_cell.size), side="right")
    for c in range(type_of_cell.size):
        n = ends[c] - starts[c]
        if n == 0:
            continue
        counts = rng.multinomial(int(n), profiles[type_of_cell[c]])
        genes[starts[c]:ends[c]] = np.repeat(names, counts)
    return genes


def _bg_molecules(n_bg: int, extent_um: float, slab_um: float | None,
                  pooled: np.ndarray, gene_names: list[str],
                  rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Uniform ambient molecules (``cell = 0``) and their genes."""
    xy = rng.uniform(0.0, extent_um, size=(n_bg, 2))
    if slab_um is not None:
        z = rng.uniform(0.0, slab_um, size=n_bg)
        xy = np.column_stack([xy, z])
    gi = rng.choice(len(gene_names), size=n_bg, p=pooled)
    return xy, np.asarray(gene_names, dtype=object)[gi]


def _capture_params(fn, env: dict, seed: int, dataset_id: str, tier: str) -> dict:
    """Actual runtime parameter values of a scenario (for ``truth.params``).

    ``env`` is the scenario's ``locals()``: every signature parameter is bound
    there, so the recorded values are what was really used, not the defaults.
    """
    names = inspect.signature(fn).parameters
    params = {k: env[k] for k in names}
    params.update(seed=seed, dataset_id=dataset_id, tier=tier)
    return params


def _assemble(*, dataset_id: str, tier: str, scenario: str, seed: int,
              params: dict, extent_um: float, slab_um: float | None,
              centers: np.ndarray, mol_xy: np.ndarray, mol_cell: np.ndarray,
              mol_genes: np.ndarray, type_of_cell: np.ndarray,
              type_names: list[str], interior_margin: float,
              r_nucleus: np.ndarray | float, bg_frac: float, bg_xy: np.ndarray,
              bg_genes: np.ndarray, scale_um: float, notes: str
              ) -> tuple[pd.DataFrame, dict]:
    """Build the contract molecule table + ``meta.json`` payload.

    ``centers`` is ``(N, 2)`` (2D scenarios) or ``(N, 3)`` (3D); ``mol_xy`` /
    ``bg_xy`` match that dimensionality.
    """
    n_cells = len(centers)
    interior_cell = (
        (centers[:, 0] >= interior_margin) & (centers[:, 0] <= extent_um - interior_margin)
        & (centers[:, 1] >= interior_margin) & (centers[:, 1] <= extent_um - interior_margin)
    )
    n_bg = len(bg_xy)

    obs_all = np.vstack([mol_xy, bg_xy])
    prior, _ = common.nucleus_prior(obs_all, centers, r_nucleus)

    cell = np.concatenate([mol_cell + 1, np.zeros(n_bg, dtype=np.int64)]).astype(np.int32)
    if n_bg:
        bg_interior = (
            (bg_xy[:, 0] >= interior_margin) & (bg_xy[:, 0] <= extent_um - interior_margin)
            & (bg_xy[:, 1] >= interior_margin) & (bg_xy[:, 1] <= extent_um - interior_margin)
        )
    else:
        bg_interior = np.zeros(0, dtype=bool)
    interior = np.concatenate([interior_cell[mol_cell], bg_interior]).astype(bool)
    celltype = np.concatenate([
        np.asarray(type_names, dtype=object)[type_of_cell[mol_cell]],
        np.full(n_bg, "", dtype=object),
    ])
    gene = np.concatenate([mol_genes, bg_genes])

    xy = obs_all
    df = pd.DataFrame({
        "x": xy[:, 0], "y": xy[:, 1], "gene": gene,
        "prior": prior, "cell": cell, "interior": interior, "celltype": celltype,
    })
    if xy.shape[1] == 3:
        df.insert(2, "z", xy[:, 2])
    df = common.make_molecules(df)

    stats = common.make_stats(df, area_um2=extent_um ** 2)
    if slab_um is not None:
        stats["volume_um3"] = float(extent_um ** 2 * slab_um)
    meta = common.make_meta(
        id=dataset_id, tier=tier,
        source={
            "generator": "benchmarks/simulate/trivial.py",
            "scenario": scenario,
            "generator_version": GENERATOR_VERSION,
            "seed": seed,
            "note": "synthetic data with exact ground truth; molecules are "
                    "uniform inside their cell (no displacement noise)",
        },
        crop={
            "bbox_um": [0.0, 0.0, extent_um, extent_um],
            "z_range_um": [0.0, float(slab_um)] if slab_um is not None else None,
            "note": "simulated domain",
        },
        stats=stats,
        difficulty={
            "cell_density": common.density_class(stats["true_cells_per_mm2"]),
            "gene_panel": common.gene_panel_class(stats["n_genes"]),
            "notes": notes,
        },
        baysor=common.make_baysor(scale_um),
        truth={
            "generator": "benchmarks/simulate/trivial.py",
            "generator_version": GENERATOR_VERSION,
            "scenario": scenario,
            "seed": seed,
            "params": params,
            "note": "truth exact by construction (no displacement); an oracle "
                    "accuracy is not computed for the trivial scenarios",
        },
    )
    return df, meta


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------

def circles_gaps(*, seed: int, dataset_id: str, tier: str = "quick",
                 n_genes: int = 100, extent_um: float = 380.0,
                 spacing_um: float = 16.0, radius_um: float = 6.1,
                 n_types: int = 4, jitter_frac: float = 0.12,
                 bg_frac: float = 0.01, r_nucleus_um: float = 3.0,
                 min_gap_um: float = 2.0, mean_molecules: float = 140.0,
                 molecules_sigma: float = 0.35, markers_per_type: int = 8,
                 marker_mass: float = 0.55) -> tuple[pd.DataFrame, dict]:
    """Round cells of equal radius on a jittered hex lattice with clear gaps."""
    rngs = common.child_rngs(seed, _RNG_STREAMS)
    centers = _lattice_with_gaps(extent_um, spacing_um, jitter_frac, rngs["geo"],
                                 min_edge_dist_um=radius_um, min_gap_um=min_gap_um,
                                 max_radius_um=radius_um)
    n_cells = len(centers)
    proportions = rngs["types"].dirichlet(np.full(n_types, 4.0))
    type_of_cell = rngs["types"].choice(n_types, size=n_cells, p=proportions).astype(np.int32)
    profiles = common.expression_profiles(
        n_types, n_genes, rngs["profiles"],
        markers_per_type=markers_per_type, marker_mass=marker_mass)
    n_c = common.sample_counts(n_cells, rngs["counts"], median=mean_molecules,
                               sigma=molecules_sigma, lo=50, hi=300)
    repeats = np.repeat(np.arange(n_cells), n_c)
    mol_xy = common.sample_in_discs_variable(
        centers[repeats], np.full(int(n_c.sum()), radius_um), rngs["pos"])
    mol_cell = repeats.astype(np.int32)
    mol_genes = _emit_genes(mol_cell, type_of_cell, profiles, n_genes, rngs["genes"])
    pooled = common.pooled_profile(profiles, proportions)
    n_bg = int(round(bg_frac / (1.0 - bg_frac) * n_c.sum()))
    bg_xy, bg_genes = _bg_molecules(n_bg, extent_um, None, pooled,
                                    _gene_names(n_genes), rngs["bg"])
    params = _capture_params(circles_gaps, locals(), seed, dataset_id, tier)
    return _assemble(
        dataset_id=dataset_id, tier=tier, scenario="circles_gaps", seed=seed,
        params=params, extent_um=extent_um, slab_um=None, centers=centers,
        mol_xy=mol_xy, mol_cell=mol_cell, mol_genes=mol_genes,
        type_of_cell=type_of_cell, type_names=[f"ct{t}" for t in range(n_types)],
        interior_margin=radius_um + 5.0, r_nucleus=r_nucleus_um,
        bg_frac=bg_frac, bg_xy=bg_xy, bg_genes=bg_genes, scale_um=radius_um,
        notes="round cells with clear gaps, ~1% uniform background")


def tiled_distinct(*, seed: int, dataset_id: str, tier: str = "quick",
                   n_genes: int = 100, extent_um: float = 380.0,
                   spacing_um: float = 14.5, jitter_frac: float = 0.15,
                   bg_frac: float = 0.01, r_nucleus_um: float = 3.0,
                   mean_molecules: float = 140.0, molecules_sigma: float = 0.35,
                   markers_per_type: int = 8) -> tuple[pd.DataFrame, dict]:
    """Gapless Voronoi tiling; neighbouring cells have different types with
    disjoint gene sets (boundaries resolvable by composition only)."""
    rngs = common.child_rngs(seed, _RNG_STREAMS)
    centers = common.hex_centers(extent_um, spacing_um, jitter_frac, rngs["geo"])
    n_cells = len(centers)
    colors = common.dsatur_color(common.delaunay_adjacency(centers))
    n_types = int(colors.max()) + 1
    type_of_cell = colors.astype(np.int32)
    profiles = common.expression_profiles(
        n_types, n_genes, rngs["profiles"], markers_per_type=markers_per_type,
        marker_mass=0.5, disjoint=True)
    n_c = common.sample_counts(n_cells, rngs["counts"], median=mean_molecules,
                               sigma=molecules_sigma, lo=50, hi=300)
    mol_xy, mol_cell = common.sample_in_voronoi(centers, n_c, extent_um, rngs["pos"])
    mol_cell = mol_cell.astype(np.int32)
    mol_genes = _emit_genes(mol_cell, type_of_cell, profiles, n_genes, rngs["genes"])
    proportions = np.bincount(type_of_cell, minlength=n_types) / n_cells
    pooled = common.pooled_profile(profiles, proportions)
    n_bg = int(round(bg_frac / (1.0 - bg_frac) * n_c.sum()))
    bg_xy, bg_genes = _bg_molecules(n_bg, extent_um, None, pooled,
                                    _gene_names(n_genes), rngs["bg"])
    params = _capture_params(tiled_distinct, locals(), seed, dataset_id, tier)
    return _assemble(
        dataset_id=dataset_id, tier=tier, scenario="tiled_distinct", seed=seed,
        params=params, extent_um=extent_um, slab_um=None, centers=centers,
        mol_xy=mol_xy, mol_cell=mol_cell, mol_genes=mol_genes,
        type_of_cell=type_of_cell, type_names=[f"ct{t}" for t in range(n_types)],
        interior_margin=spacing_um, r_nucleus=r_nucleus_um, bg_frac=bg_frac,
        bg_xy=bg_xy, bg_genes=bg_genes,
        scale_um=float(np.sqrt((extent_um ** 2 / n_cells) / np.pi)),
        notes=f"Voronoi tiling, {n_types} DSATUR colours, disjoint gene sets")


def tiled_same(*, seed: int, dataset_id: str, tier: str = "quick",
               n_genes: int = 100, extent_um: float = 380.0,
               spacing_um: float = 14.5, jitter_frac: float = 0.15,
               bg_frac: float = 0.01, r_nucleus_um: float = 3.0,
               mean_molecules: float = 140.0,
               molecules_sigma: float = 0.35) -> tuple[pd.DataFrame, dict]:
    """The ``tiled_distinct`` geometry with all cells of one identical type:
    boundaries are unresolvable by composition (hard control)."""
    rngs = common.child_rngs(seed, _RNG_STREAMS)
    centers = common.hex_centers(extent_um, spacing_um, jitter_frac, rngs["geo"])
    n_cells = len(centers)
    type_of_cell = np.zeros(n_cells, dtype=np.int32)
    profiles = common.expression_profiles(1, n_genes, rngs["profiles"],
                                          markers_per_type=8, marker_mass=0.5)
    n_c = common.sample_counts(n_cells, rngs["counts"], median=mean_molecules,
                               sigma=molecules_sigma, lo=50, hi=300)
    mol_xy, mol_cell = common.sample_in_voronoi(centers, n_c, extent_um, rngs["pos"])
    mol_cell = mol_cell.astype(np.int32)
    mol_genes = _emit_genes(mol_cell, type_of_cell, profiles, n_genes, rngs["genes"])
    pooled = common.pooled_profile(profiles, np.ones(1))
    n_bg = int(round(bg_frac / (1.0 - bg_frac) * n_c.sum()))
    bg_xy, bg_genes = _bg_molecules(n_bg, extent_um, None, pooled,
                                    _gene_names(n_genes), rngs["bg"])
    params = _capture_params(tiled_same, locals(), seed, dataset_id, tier)
    return _assemble(
        dataset_id=dataset_id, tier=tier, scenario="tiled_same", seed=seed,
        params=params, extent_um=extent_um, slab_um=None, centers=centers,
        mol_xy=mol_xy, mol_cell=mol_cell, mol_genes=mol_genes,
        type_of_cell=type_of_cell, type_names=["ct0"],
        interior_margin=spacing_um, r_nucleus=r_nucleus_um, bg_frac=bg_frac,
        bg_xy=bg_xy, bg_genes=bg_genes,
        scale_um=float(np.sqrt((extent_um ** 2 / n_cells) / np.pi)),
        notes="single cell type everywhere: composition gives no boundary signal")


def mixed_sizes(*, seed: int, dataset_id: str, tier: str = "quick",
                n_genes: int = 100, extent_um: float = 500.0,
                spacing_um: float = 20.5, radius_small_um: float = 3.0,
                radius_large_um: float = 8.8, p_small: float = 0.45,
                jitter_frac: float = 0.05, bg_frac: float = 0.01,
                r_nucleus_um: float = 2.5, min_gap_um: float = 1.8,
                median_small: float = 70.0, median_large: float = 190.0,
                molecules_sigma: float = 0.35,
                markers_per_type: int = 8) -> tuple[pd.DataFrame, dict]:
    """Small immune-like and large tumour-like round cells with small gaps."""
    rngs = common.child_rngs(seed, _RNG_STREAMS)
    centers = _lattice_with_gaps(extent_um, spacing_um, jitter_frac, rngs["geo"],
                                 min_edge_dist_um=radius_large_um,
                                 min_gap_um=min_gap_um,
                                 max_radius_um=radius_large_um)
    n_cells = len(centers)
    is_small = rngs["types"].random(n_cells) < p_small
    type_of_cell = np.where(is_small, 0, 1).astype(np.int32)  # 0 immune, 1 tumour
    radius = np.where(is_small, radius_small_um, radius_large_um)
    profiles = common.expression_profiles(
        2, n_genes, rngs["profiles"], markers_per_type=markers_per_type,
        marker_mass=0.6)
    n_small = int(is_small.sum())
    n_c = np.empty(n_cells, dtype=np.int64)
    n_c[is_small] = common.sample_counts(n_small, rngs["counts"],
                                         median=median_small, sigma=molecules_sigma,
                                         lo=40, hi=160)
    n_c[~is_small] = common.sample_counts(n_cells - n_small, rngs["counts"],
                                          median=median_large, sigma=molecules_sigma,
                                          lo=80, hi=300)
    repeats = np.repeat(np.arange(n_cells), n_c)
    mol_xy = common.sample_in_discs_variable(centers[repeats], radius[repeats],
                                             rngs["pos"])
    mol_cell = repeats.astype(np.int32)
    mol_genes = _emit_genes(mol_cell, type_of_cell, profiles, n_genes, rngs["genes"])
    proportions = np.array([p_small, 1.0 - p_small])
    pooled = common.pooled_profile(profiles, proportions)
    n_bg = int(round(bg_frac / (1.0 - bg_frac) * n_c.sum()))
    bg_xy, bg_genes = _bg_molecules(n_bg, extent_um, None, pooled,
                                    _gene_names(n_genes), rngs["bg"])
    params = _capture_params(mixed_sizes, locals(), seed, dataset_id, tier)
    return _assemble(
        dataset_id=dataset_id, tier=tier, scenario="mixed_sizes", seed=seed,
        params=params, extent_um=extent_um, slab_um=None, centers=centers,
        mol_xy=mol_xy, mol_cell=mol_cell, mol_genes=mol_genes,
        type_of_cell=type_of_cell, type_names=["immune", "tumour"],
        interior_margin=radius_large_um + 5.0, r_nucleus=r_nucleus_um,
        bg_frac=bg_frac, bg_xy=bg_xy, bg_genes=bg_genes,
        scale_um=float(np.sqrt(np.mean(radius ** 2))),
        notes="bimodal radii (3 um immune-like, 8.8 um tumour-like), small gaps")


def sparse_noisy(*, seed: int, dataset_id: str, tier: str = "quick",
                 n_genes: int = 100, extent_um: float = 600.0,
                 spacing_um: float = 40.0, radius_um: float = 8.0,
                 n_types: int = 3, jitter_frac: float = 0.1,
                 bg_frac: float = 0.25, r_nucleus_um: float = 3.0,
                 mean_molecules: float = 140.0, molecules_sigma: float = 0.35,
                 markers_per_type: int = 8,
                 marker_mass: float = 0.55) -> tuple[pd.DataFrame, dict]:
    """Sparse cells over a large area with 20-30% uniform background noise."""
    rngs = common.child_rngs(seed, _RNG_STREAMS)
    centers = common.hex_centers(extent_um, spacing_um, jitter_frac, rngs["geo"],
                                 min_edge_dist_um=radius_um)
    n_cells = len(centers)
    proportions = rngs["types"].dirichlet(np.full(n_types, 4.0))
    type_of_cell = rngs["types"].choice(n_types, size=n_cells, p=proportions).astype(np.int32)
    profiles = common.expression_profiles(
        n_types, n_genes, rngs["profiles"], markers_per_type=markers_per_type,
        marker_mass=marker_mass)
    n_c = common.sample_counts(n_cells, rngs["counts"], median=mean_molecules,
                               sigma=molecules_sigma, lo=50, hi=300)
    repeats = np.repeat(np.arange(n_cells), n_c)
    mol_xy = common.sample_in_discs_variable(
        centers[repeats], np.full(int(n_c.sum()), radius_um), rngs["pos"])
    mol_cell = repeats.astype(np.int32)
    mol_genes = _emit_genes(mol_cell, type_of_cell, profiles, n_genes, rngs["genes"])
    pooled = common.pooled_profile(profiles, proportions)
    n_bg = int(round(bg_frac / (1.0 - bg_frac) * n_c.sum()))
    bg_xy, bg_genes = _bg_molecules(n_bg, extent_um, None, pooled,
                                    _gene_names(n_genes), rngs["bg"])
    params = _capture_params(sparse_noisy, locals(), seed, dataset_id, tier)
    return _assemble(
        dataset_id=dataset_id, tier=tier, scenario="sparse_noisy", seed=seed,
        params=params, extent_um=extent_um, slab_um=None, centers=centers,
        mol_xy=mol_xy, mol_cell=mol_cell, mol_genes=mol_genes,
        type_of_cell=type_of_cell, type_names=[f"ct{t}" for t in range(n_types)],
        interior_margin=radius_um + 5.0, r_nucleus=r_nucleus_um, bg_frac=bg_frac,
        bg_xy=bg_xy, bg_genes=bg_genes, scale_um=radius_um,
        notes=f"sparse packing, {bg_frac:.0%} uniform background (cell = 0)")


def circles_gaps_3d(*, seed: int, dataset_id: str, tier: str = "quick",
                    n_genes: int = 100, extent_um: float = 380.0,
                    spacing_um: float = 16.0, radius_um: float = 6.1,
                    slab_um: float = 10.0, n_types: int = 4,
                    jitter_frac: float = 0.12, bg_frac: float = 0.01,
                    r_nucleus_um: float = 3.0, min_gap_um: float = 2.0,
                    mean_molecules: float = 140.0, molecules_sigma: float = 0.35,
                    markers_per_type: int = 8,
                    marker_mass: float = 0.55) -> tuple[pd.DataFrame, dict]:
    """3D variant: spheres in a thin slab; cells are cut by the section."""
    rngs = common.child_rngs(seed, _RNG_STREAMS)
    centers = _lattice_with_gaps(extent_um, spacing_um, jitter_frac, rngs["geo"],
                                 min_edge_dist_um=radius_um, min_gap_um=min_gap_um,
                                 max_radius_um=radius_um)
    n_cells = len(centers)
    center_z = rngs["geo"].uniform(0.0, slab_um, size=n_cells)
    centers3 = np.column_stack([centers, center_z])
    proportions = rngs["types"].dirichlet(np.full(n_types, 4.0))
    type_of_cell = rngs["types"].choice(n_types, size=n_cells, p=proportions).astype(np.int32)
    profiles = common.expression_profiles(
        n_types, n_genes, rngs["profiles"], markers_per_type=markers_per_type,
        marker_mass=marker_mass)
    n_c = common.sample_counts(n_cells, rngs["counts"], median=mean_molecules,
                               sigma=molecules_sigma, lo=50, hi=300)
    repeats = np.repeat(np.arange(n_cells), n_c)
    mol_xy = common.sample_in_slab_spheres(
        centers3[repeats], radius_um, slab_um, rngs["pos"])
    mol_cell = repeats.astype(np.int32)
    mol_genes = _emit_genes(mol_cell, type_of_cell, profiles, n_genes, rngs["genes"])
    pooled = common.pooled_profile(profiles, proportions)
    n_bg = int(round(bg_frac / (1.0 - bg_frac) * n_c.sum()))
    bg_xy, bg_genes = _bg_molecules(n_bg, extent_um, slab_um, pooled,
                                    _gene_names(n_genes), rngs["bg"])
    params = _capture_params(circles_gaps_3d, locals(), seed, dataset_id, tier)
    return _assemble(
        dataset_id=dataset_id, tier=tier, scenario="circles_gaps_3d", seed=seed,
        params=params, extent_um=extent_um, slab_um=slab_um, centers=centers3,
        mol_xy=mol_xy, mol_cell=mol_cell, mol_genes=mol_genes,
        type_of_cell=type_of_cell, type_names=[f"ct{t}" for t in range(n_types)],
        interior_margin=radius_um + 5.0, r_nucleus=r_nucleus_um, bg_frac=bg_frac,
        bg_xy=bg_xy, bg_genes=bg_genes, scale_um=radius_um,
        notes=f"spheres in a {slab_um:g} um slab; cells cut by the section")


SCENARIOS = {
    "circles_gaps": circles_gaps,
    "tiled_distinct": tiled_distinct,
    "tiled_same": tiled_same,
    "mixed_sizes": mixed_sizes,
    "sparse_noisy": sparse_noisy,
    "circles_gaps_3d": circles_gaps_3d,
}


def generate(scenario: str, *, seed: int, dataset_id: str, tier: str = "quick",
             params: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """Dispatch by scenario name; ``params`` are the scenario's keyword args."""
    if scenario not in SCENARIOS:
        raise KeyError(f"unknown scenario {scenario!r}; have {sorted(SCENARIOS)}")
    kwargs = dict(params or {})
    return SCENARIOS[scenario](seed=seed, dataset_id=dataset_id, tier=tier, **kwargs)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _str2bool(v: str) -> bool:
    if v.lower() in ("1", "true", "yes"):
        return True
    if v.lower() in ("0", "false", "no"):
        return False
    raise argparse.ArgumentTypeError(f"bad bool {v!r}")


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("scenario", choices=sorted(SCENARIOS))
    known, _ = pre.parse_known_args(argv)
    fn = SCENARIOS[known.scenario]

    parser = argparse.ArgumentParser(
        description="Generate one trivial simulated dataset "
                    "(contract: benchmarks/README.md).")
    parser.add_argument("scenario", choices=sorted(SCENARIOS))
    parser.add_argument("--id", required=True, help="dataset id (directory name)")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--tier", choices=["quick", "full"], default="quick")
    parser.add_argument("-o", "--out-dir", required=True,
                        help="output directory ($BAYSOR_BENCH_DATA/sim/<id>)")
    for name, prm in inspect.signature(fn).parameters.items():
        if prm.default is inspect.Parameter.empty:
            continue
        flag = "--" + name.replace("_", "-")
        if isinstance(prm.default, bool):
            parser.add_argument(flag, type=_str2bool, default=prm.default)
        else:
            parser.add_argument(flag, type=type(prm.default), default=prm.default)
    args = parser.parse_args(argv)

    kwargs = {name: getattr(args, name)
              for name in inspect.signature(fn).parameters if hasattr(args, name)}
    df, meta = generate(known.scenario, seed=args.seed, dataset_id=args.id,
                        tier=args.tier, params=kwargs)
    hashes = common.write_dataset(args.out_dir, df, meta)
    print(json.dumps({"id": args.id, "out": args.out_dir,
                      "n_molecules": int(len(df)),
                      "n_genes": int(df["gene"].nunique()),
                      "sha256": hashes}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
