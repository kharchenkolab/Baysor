"""Shared helpers for the simulated datasets of the Baysor benchmark suite.

Everything here is deterministic given a seed and follows the dataset contract
in ``benchmarks/README.md``: molecules live in ``molecules.parquet`` (sorted by
``(y, x)``), metadata in ``meta.json``.

The module provides:
  * geometry helpers (jittered hex lattices, Delaunay adjacency, DSATUR
    colouring of the adjacency graph);
  * expression-profile construction (shared background + per-type marker
    blocks, lognormal weights; optionally strictly disjoint gene sets);
  * molecule placement (uniform in a disc, in a Voronoi region, or in a
    sphere clipped to a slab);
  * the nucleus prior column (molecule within ``r_nucleus`` of a true cell
    centre);
  * contract assembly/writing (dtypes, sorting, stats, difficulty classes).
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Tier budgets (see benchmarks/README.md)
# ---------------------------------------------------------------------------
TIER_BUDGET = {"quick": 150_000, "full": 3_000_000}


def data_root() -> Path:
    """``$BAYSOR_BENCH_DATA`` (default ``<repo>/.bench-data``)."""
    env = os.environ.get("BAYSOR_BENCH_DATA")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[2] / ".bench-data"


def cache_dir() -> Path:
    """Simulation cache directory ``$BAYSOR_BENCH_DATA/cache/sim``."""
    return data_root() / "cache" / "sim"

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def hex_centers(extent_um: float, spacing_um: float, jitter_frac: float,
                rng: np.random.Generator, min_edge_dist_um: float = 0.0) -> np.ndarray:
    """Jittered hexagonal lattice of cell centres inside ``[0, extent_um] ^ 2``.

    Rows are ``spacing * sqrt(3)/2`` apart; odd rows are offset by half a
    spacing.  Centres are jittered by ``jitter_frac * spacing`` (Gaussian) and
    clipped to be at least ``min_edge_dist_um`` from every domain edge (used by
    the round-cell scenarios so whole cells lie inside the domain).

    Returns an ``(N, 2)`` float64 array ordered row-major (top row first).
    """
    dy = spacing_um * np.sqrt(3.0) / 2.0
    n_rows = max(int(np.floor(extent_um / dy)), 1)
    n_cols = max(int(np.floor(extent_um / spacing_um)), 1)
    pts = []
    for i in range(n_rows):
        y = (i + 0.5) * dy
        offset = 0.5 * spacing_um if i % 2 else 0.0
        for j in range(n_cols):
            x = (j + 0.5) * spacing_um + offset
            if x < extent_um:
                pts.append((x, y))
    centers = np.asarray(pts, dtype=np.float64)
    centers += rng.normal(0.0, jitter_frac * spacing_um, size=centers.shape)
    lo, hi = min_edge_dist_um, extent_um - min_edge_dist_um
    centers = np.clip(centers, lo, hi)
    return centers


def delaunay_adjacency(centers: np.ndarray) -> list[set[int]]:
    """Neighbour lists of the Delaunay triangulation (= Voronoi adjacency)."""
    from scipy.spatial import Delaunay

    tri = Delaunay(centers)
    n = len(centers)
    adj = [set() for _ in range(n)]
    for simplex in tri.simplices:
        for a in simplex:
            for b in simplex:
                if a != b:
                    adj[int(a)].add(int(b))
    return adj


def dsatur_color(adjacency: list[set[int]]) -> np.ndarray:
    """Deterministic DSATUR colouring of an undirected graph.

    Repeatedly colours the uncoloured vertex with the largest saturation
    degree (ties: larger degree, then smaller index).  Returns an ``int32``
    colour per vertex; adjacent vertices never share a colour.  Implemented
    with a lazy heap so large Delaunay graphs (full tier) stay fast.
    """
    import heapq

    n = len(adjacency)
    colour = np.full(n, -1, dtype=np.int32)
    sat: list[set[int]] = [set() for _ in range(n)]
    # heap entries: (-saturation, -degree, vertex); stale entries are skipped
    heap = [(0, -len(adjacency[v]), v) for v in range(n)]
    heapq.heapify(heap)
    while heap:
        neg_sat, _neg_deg, v = heapq.heappop(heap)
        if colour[v] >= 0:
            continue
        if -neg_sat != len(sat[v]):  # stale entry
            continue
        used = {colour[u] for u in adjacency[v] if colour[u] >= 0}
        c = 0
        while c in used:
            c += 1
        colour[v] = c
        for u in adjacency[v]:
            if colour[u] < 0:
                sat[u].add(c)
                heapq.heappush(heap, (-len(sat[u]), -len(adjacency[u]), u))
    if (colour < 0).any():
        raise RuntimeError("DSATUR left uncoloured vertices")
    return colour


# ---------------------------------------------------------------------------
# Expression profiles
# ---------------------------------------------------------------------------


def expression_profiles(n_types: int, n_genes: int, rng: np.random.Generator,
                        *, markers_per_type: int = 8, marker_mass: float = 0.55,
                        disjoint: bool = False) -> np.ndarray:
    """Per-type normalised expression compositions ``p`` of shape ``(K, G)``.

    Realistic shape: a few strongly expressed marker genes per type plus a
    shared low-level background over the rest of the panel; weights inside
    each block are lognormal (a Dirichlet-like draw), rows sum to 1.

    With ``disjoint=True`` each type owns a contiguous, strictly disjoint
    block of the panel (marker genes inside the block, low-level background
    over the remainder of the block) so the gene identifies the type exactly.
    """
    if not (0.0 < marker_mass < 1.0):
        raise ValueError("marker_mass must be in (0, 1)")
    p = np.zeros((n_types, n_genes), dtype=np.float64)
    if disjoint:
        if n_types > n_genes:
            raise ValueError("disjoint panels need at least one gene per type")
        blocks = np.array_split(np.arange(n_genes), n_types)
        for t, block in enumerate(blocks):
            m = int(min(markers_per_type, max(1, len(block) // 2)))
            w = np.zeros(len(block), dtype=np.float64)
            w[:m] = rng.lognormal(0.0, 0.5, size=m)
            w[m:] = rng.lognormal(0.0, 0.5, size=len(block) - m) * 0.05
            block_p = w / w.sum()
            p[t, block] = block_p
    else:
        background = rng.lognormal(0.0, 0.5, size=n_genes)
        background = background / background.sum()
        m = int(min(markers_per_type, n_genes // max(n_types, 1)))
        # distinct marker blocks across types (shuffled gene order, split)
        perm = rng.permutation(n_genes)
        for t in range(n_types):
            markers = perm[t * m:(t + 1) * m]
            marker_p = np.zeros(n_genes, dtype=np.float64)
            marker_p[markers] = rng.lognormal(0.0, 0.5, size=len(markers))
            marker_p[markers] /= marker_p[markers].sum()
            p[t] = marker_mass * marker_p + (1.0 - marker_mass) * background
            p[t] /= p[t].sum()
    if not np.allclose(p.sum(axis=1), 1.0):
        raise ValueError("profiles must sum to 1")
    return p


def pooled_profile(p: np.ndarray, proportions: np.ndarray) -> np.ndarray:
    """Tissue-average profile (mixture over types), used for ambient molecules."""
    pooled = proportions @ p
    return pooled / pooled.sum()


def sample_counts(n_cells: int, rng: np.random.Generator, *, median: float = 140.0,
                  sigma: float = 0.35, lo: int = 50, hi: int = 300) -> np.ndarray:
    """Lognormal molecules-per-cell counts clipped to ``[lo, hi]`` (Xenium-like)."""
    n = rng.lognormal(np.log(median), sigma, size=n_cells)
    return np.clip(np.round(n), lo, hi).astype(np.int64)


# ---------------------------------------------------------------------------
# Molecule placement
# ---------------------------------------------------------------------------


def sample_in_discs(centers: np.ndarray, radius_um: float,
                    rng: np.random.Generator) -> np.ndarray:
    """Uniform points, one per row position, inside discs of equal radius."""
    n = len(centers)
    r = radius_um * np.sqrt(rng.random(n))
    theta = rng.uniform(0.0, 2.0 * np.pi, size=n)
    return centers + np.column_stack([r * np.cos(theta), r * np.sin(theta)])


def sample_in_discs_variable(centers: np.ndarray, radius_um: np.ndarray,
                             rng: np.random.Generator) -> np.ndarray:
    """Uniform points inside discs of per-cell radius (one point per row)."""
    n = len(centers)
    r = radius_um * np.sqrt(rng.random(n))
    theta = rng.uniform(0.0, 2.0 * np.pi, size=n)
    return centers + np.column_stack([r * np.cos(theta), r * np.sin(theta)])


def sample_in_voronoi(centers: np.ndarray, n_per_cell: np.ndarray, extent_um: float,
                      rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Uniform-in-Voronoi-region points for every cell.

    Draws uniform points over the whole domain and keeps the first
    ``n_per_cell[c]`` points whose nearest centre is ``c`` (rejection in
    batches), which is uniform inside each region intersected with the
    domain.  Returns ``(xy, cell_index)``.
    """
    from scipy.spatial import cKDTree

    total = int(n_per_cell.sum())
    need = n_per_cell.astype(np.int64).copy()
    kept_xy: list[np.ndarray] = []
    kept_cell: list[np.ndarray] = []
    tree = cKDTree(centers)
    guard = 0
    while need.sum() > 0:
        guard += 1
        if guard > 200:
            raise RuntimeError("Voronoi rejection sampling failed to converge")
        batch = int(max(need.sum() * 1.4, 10_000))
        pts = rng.uniform(0.0, extent_um, size=(batch, 2))
        owner = tree.query(pts, k=1)[1].astype(np.int64)
        # take, per cell, at most what is still needed
        order = np.argsort(owner, kind="stable")
        owner_s = owner[order]
        pts_s = pts[order]
        starts = np.searchsorted(owner_s, np.arange(len(centers)), side="left")
        ends = np.searchsorted(owner_s, np.arange(len(centers)), side="right")
        for c in np.nonzero(need > 0)[0]:
            k = int(min(need[c], ends[c] - starts[c]))
            if k <= 0:
                continue
            sl = slice(starts[c], starts[c] + k)
            kept_xy.append(pts_s[sl])
            kept_cell.append(np.full(k, c, dtype=np.int64))
            need[c] -= k
    xy = np.concatenate(kept_xy, axis=0)
    cell = np.concatenate(kept_cell, axis=0)
    assert len(xy) == total
    # restore deterministic cell-major order
    order = np.argsort(cell, kind="stable")
    return xy[order], cell[order]


def sample_in_slab_spheres(centers: np.ndarray, radius_um: float, slab_um: float,
                           rng: np.random.Generator) -> np.ndarray:
    """Uniform points inside spheres of equal radius clipped to ``z in [0, slab]``.

    Centres are ``(N, 3)`` with ``z`` already inside the slab.  Rejection
    sampling keeps the conditional uniform distribution inside the clipped
    sphere.  Returns ``(N, 3)``.
    """
    from scipy.spatial import cKDTree  # noqa: F401  (kept for symmetry/API clarity)

    n = len(centers)
    out = np.empty((n, 3), dtype=np.float64)
    pending = np.arange(n)
    guard = 0
    while pending.size:
        guard += 1
        if guard > 500:
            raise RuntimeError("slab-sphere rejection sampling failed to converge")
        m = pending.size
        # uniform in ball: direction ~ N(0, I) normalised, radius ~ u^(1/3)
        dirs = rng.normal(size=(m, 3))
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        rad = radius_um * rng.random(m) ** (1.0 / 3.0)
        pts = centers[pending] + dirs * rad[:, None]
        ok = (pts[:, 2] >= 0.0) & (pts[:, 2] <= slab_um)
        out[pending[ok]] = pts[ok]
        pending = pending[~ok]
    return out


# ---------------------------------------------------------------------------
# Nucleus prior
# ---------------------------------------------------------------------------


def nucleus_prior(xy: np.ndarray, centers: np.ndarray,
                  radius_um: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """Prior label per molecule: 0 = no prior, else 1-based true centre index.

    A molecule gets the label of its nearest true cell centre when it lies
    within that centre's nucleus radius (3D when ``xy`` has three columns).
    ``radius_um`` may be a scalar or a per-centre array.
    """
    from scipy.spatial import cKDTree

    tree = cKDTree(centers)
    dist, idx = tree.query(xy, k=1)
    r = np.broadcast_to(np.asarray(radius_um, dtype=np.float64), (len(centers),))
    label = np.where(dist <= r[idx], idx + 1, 0).astype(np.int32)
    return label, idx.astype(np.int32)


def imperfect_nucleus_prior(xy: np.ndarray, centers: np.ndarray,
                            radius_um: np.ndarray | float, *, seed: int,
                            miss_frac: float = 0.2,
                            shift_um: tuple[float, float] = (1.0, 2.0),
                            merge_frac: float = 0.05
                            ) -> tuple[np.ndarray, dict]:
    """Prior label per molecule from an *imperfect* nucleus segmentation.

    Models a realistic vendor nucleus segmentation with three defects,
    applied in this order (all draws from an RNG seeded by ``seed``, so
    molecule placement elsewhere is unaffected):

    1. **missed nuclei** — ``miss_frac`` of the cells get no nucleus at all
       (no prior label anywhere in that cell);
    2. **shifted nuclei** — every remaining nucleus centre is displaced by a
       random direction and a magnitude drawn uniformly from ``shift_um``
       (default 1–2 µm);
    3. **merged nuclei** — ``merge_frac`` of the cells have their prior label
       remapped onto their nearest *unmerged* kept neighbour, i.e. one prior
       label covers two adjacent nuclei (sources are chosen first, targets
       only among non-sources, so the remap has no chains/cycles).

    Returns ``(label, info)``: the 1-based original cell id of the nucleus
    that claims the molecule (0 = no prior), and a dict describing the
    applied defects (recorded in ``meta.truth``).  Ground truth elsewhere is
    untouched: the returned labels are the *prior column only*.
    """
    from scipy.spatial import cKDTree

    if not (0.0 <= miss_frac < 1.0):
        raise ValueError("miss_frac must be in [0, 1)")
    if not (0.0 <= merge_frac < 1.0):
        raise ValueError("merge_frac must be in [0, 1)")
    if not (0.0 <= shift_um[0] <= shift_um[1]):
        raise ValueError("shift_um must be an ordered non-negative pair")
    n = len(centers)
    rng = np.random.default_rng(seed)

    # 1) missed nuclei
    kept = np.ones(n, dtype=bool)
    n_miss = int(round(miss_frac * n))
    n_miss = min(n_miss, max(n - 2, 0))  # keep >= 2 nuclei for merging
    if n_miss:
        kept[rng.choice(n, n_miss, replace=False)] = False
    idx_kept = np.nonzero(kept)[0]

    # 2) shifted nuclei (applied to every centre; missed ones are unused).
    #    The shift lies in the x-y plane (nuclei slide within the section);
    #    z components, if any, stay put.
    theta = rng.uniform(0.0, 2.0 * np.pi, size=n)
    mag = rng.uniform(shift_um[0], shift_um[1], size=n)
    dim = np.asarray(centers).shape[1]
    shift = np.zeros((n, dim), dtype=np.float64)
    shift[:, 0] = np.cos(theta) * mag
    shift[:, 1] = np.sin(theta) * mag
    shifted = np.asarray(centers, dtype=np.float64) + shift

    r = np.broadcast_to(np.asarray(radius_um, dtype=np.float64), (n,))
    tree = cKDTree(shifted[idx_kept])
    dist, j = tree.query(xy, k=1)
    label = np.where(dist <= r[idx_kept][j], idx_kept[j] + 1, 0).astype(np.int32)

    # 3) merged nuclei: source labels are remapped onto their nearest kept,
    #    unmerged neighbour (targets are never sources -> no chains)
    n_merge = int(round(merge_frac * n))
    n_merge = min(n_merge, max(len(idx_kept) - 2, 0))
    merge_src: np.ndarray = np.empty(0, dtype=np.int64)
    n_merged = 0
    if n_merge:
        merge_src = rng.choice(idx_kept, size=n_merge, replace=False)
        src_set = set(int(i) for i in merge_src)
        remap = np.arange(n + 1, dtype=np.int64)  # 1-based labels; 0 stays 0
        other_tree = cKDTree(shifted)
        for s in merge_src:
            _, order = other_tree.query(shifted[s], k=n)
            target = next((int(order[k]) for k in range(n)
                           if order[k] != s and int(order[k]) not in src_set
                           and kept[order[k]]), None)
            if target is None:  # no valid neighbour (tiny fields)
                continue
            remap[int(s) + 1] = target + 1
        label = remap[label].astype(np.int32)
        n_merged = int(np.sum(remap[1:] != np.arange(1, n + 1)))

    info = {
        "kind": "imperfect_nucleus",
        "seed": int(seed),
        "miss_frac": float(miss_frac),
        "n_nuclei": int(n),
        "n_missed": int(n_miss),
        "shift_um": [float(shift_um[0]), float(shift_um[1])],
        "merge_frac": float(merge_frac),
        "n_merge_sources": int(len(merge_src)),
        "n_merged": n_merged,
        "nucleus_radius_um": (float(radius_um) if np.isscalar(radius_um)
                              else np.asarray(radius_um, dtype=float).tolist()),
    }
    return label, info


# ---------------------------------------------------------------------------
# Contract assembly / IO
# ---------------------------------------------------------------------------


def gene_panel_class(n_genes: int) -> str:
    if n_genes < 50:
        return "tiny"
    if n_genes < 250:
        return "small"
    if n_genes < 700:
        return "medium"
    if n_genes <= 2000:
        return "large"
    return "huge"


def density_class(cells_per_mm2: float) -> str:
    if cells_per_mm2 < 2500:
        return "sparse"
    if cells_per_mm2 <= 7000:
        return "medium"
    return "dense"


def make_molecules(df: pd.DataFrame) -> pd.DataFrame:
    """Finalise the molecule table: dtypes, column order, sort by (y, x)."""
    cols = {"x": "float64", "y": "float64"}
    if "z" in df.columns:
        cols["z"] = "float64"
    cols.update({
        "gene": "object",
        "prior": "int32",
        "cell": "int32",
        "interior": "bool",
        "celltype": "object",
    })
    out = df.copy()
    for col, dtype in cols.items():
        if col not in out.columns:
            raise ValueError(f"missing required column {col!r}")
        out[col] = out[col].astype(dtype)
    out = out[[*cols]]
    out = out.sort_values(["y", "x"], kind="stable").reset_index(drop=True)
    return out


def make_meta(*, id: str, tier: str, source: dict, crop: dict, stats: dict,
              difficulty: dict, baysor: dict, truth: dict | None) -> dict:
    """Assemble the ``meta.json`` payload following the dataset contract."""
    if tier not in TIER_BUDGET:
        raise ValueError(f"tier must be one of {sorted(TIER_BUDGET)}")
    if stats["n_molecules"] > TIER_BUDGET[tier]:
        raise ValueError(
            f"{id}: {stats['n_molecules']} molecules exceeds {tier} budget "
            f"{TIER_BUDGET[tier]}")
    return {
        "id": id,
        "kind": "sim",
        "tier": tier,
        "platform": "simulated",
        "source": source,
        "crop": crop,
        "stats": stats,
        "difficulty": difficulty,
        "baysor": baysor,
        "images": [],
        "truth": truth,
    }


def make_stats(df: pd.DataFrame, area_um2: float,
               n_true_cells: int | None = None) -> dict:
    n_molecules = int(len(df))
    n_cells = int(df["cell"].max()) if n_true_cells is None else int(n_true_cells)
    return {
        "n_molecules": n_molecules,
        "n_genes": int(df["gene"].nunique()),
        "area_um2": float(area_um2),
        "molecules_per_um2": float(n_molecules / area_um2),
        "n_vendor_cells": 0,
        "vendor_cells_per_mm2": 0.0,
        "n_true_cells": n_cells,
        "true_cells_per_mm2": float(n_cells / area_um2 * 1e6),
    }


def make_baysor(scale_um: float, *, min_molecules_per_cell: int = 20,
                prior: str = "column", prior_confidence: float = 0.5,
                extra_args: list | None = None) -> dict:
    """Baysor invocation parameters stored in ``meta.json``.

    Simulated datasets carry no TOML config (``config`` is null): every
    parameter the runner needs lives in this block and is passed as explicit
    CLI flags.
    """
    return {
        "scale_um": float(scale_um),
        "scale_std": "25%",
        "min_molecules_per_cell": int(min_molecules_per_cell),
        "prior": prior,
        "prior_confidence": float(prior_confidence),
        "config": None,
        "extra_args": list(extra_args or []),
    }


def write_dataset(out_dir: str | os.PathLike, df: pd.DataFrame,
                  meta: dict) -> dict[str, str]:
    """Write ``molecules.parquet`` + ``meta.json``; returns their SHA-256 sums."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    parquet = out / "molecules.parquet"
    js = out / "meta.json"
    df.to_parquet(parquet, index=False, compression="snappy")
    with open(js, "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")
    return {
        "molecules.parquet": sha256(parquet),
        "meta.json": sha256(js),
    }


def sha256(path: str | os.PathLike, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def child_rngs(seed: int, names: list[str]) -> dict[str, np.random.Generator]:
    """Independent named RNG streams derived from one master seed.

    Deriving streams by name keeps downstream draws independent of call order:
    changing the panel size, for instance, does not perturb the geometry.
    """
    ss = np.random.SeedSequence(seed)
    children = ss.spawn(len(names))
    return {name: np.random.default_rng(c) for name, c in zip(names, children)}
