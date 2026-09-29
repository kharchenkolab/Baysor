#!/usr/bin/env python
"""Deliberately degrade a molecule-to-cell segmentation for audit validation.

Two degradation modes are provided, both operating on a tabular molecule table
with a vendor cell column (empty string = background / unassigned):

``border``
    Reassign a fraction of the most border-adjacent molecules to the nearest
    neighbouring cell. For every assigned molecule the distance to the nearest
    molecule of a *different* cell is computed; molecules are ranked by that
    distance and the top ``--fraction`` of all assigned molecules (the smallest
    cross-cell distances, i.e. the cell-border molecules) are re-assigned to the
    cell owning that nearest foreign molecule. The selection is a prefix of one
    deterministic ranking, so the 30% variant is a strict superset of the 10%
    variant and the degradation is monotone by construction.

``dilate``
    Dilate every cell by ``--distance`` um: each cell's molecules are enclosed
    in their convex hull and every molecule (assigned or background) within
    ``--distance`` um of a hull is claimed by the nearest such hull, ties broken
    by cell id. Molecules near a shared border flip to the neighbouring cell
    that reaches them first; background molecules within the band are absorbed.
    Molecules farther than the distance from every hull keep their label.

Output is an assignment parquet in the row order of the input molecules with
columns ``cell`` (int32, 0 = unassigned) and ``cell_label`` (the original cell
id string), which ``audit.py --assignment`` understands.

Usage:
    python degrade.py --molecules molecules.parquet --cell-column cell_vendor \
        --mode border --fraction 0.1 --out border10.parquet --report border10.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

# Same values the cellAdmix tabular loader treats as "no cell"; molecules with
# these ids are background, never a cell (imported from audit.py to keep one
# source of truth).
from audit import UNASSIGNED_TOKENS


def load_assignment_table(molecules_path: Path | str, cell_column: str) -> pd.DataFrame:
    """Read x, y and the vendor cell column of a molecule table.

    Loader-specific unassigned tokens (e.g. Xenium's "UNASSIGNED") are
    normalized to the empty string = background.
    """
    frame = pd.read_parquet(molecules_path, columns=["x", "y", cell_column])
    frame = frame.rename(columns={cell_column: "cell_vendor"})
    labels = frame["cell_vendor"].fillna("").astype(str)
    frame["cell_vendor"] = labels.mask(labels.isin(UNASSIGNED_TOKENS), "")
    return frame


def cell_codes(labels: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    """Map cell id strings to int codes (0 = unassigned); codes start at 1.

    Implemented with a single ``factorize`` pass — the obvious loop over the
    unique labels compares the full array once per label and is
    O(n_cells * n_molecules) on real crops.
    """
    labels = np.asarray(labels, dtype=object)
    order = sorted(pd.unique(labels[labels != ""]).tolist())
    raw_map = {label: i + 1 for i, label in enumerate(order)}
    mapping: dict[str, int] = {str(k): v for k, v in raw_map.items()}
    codes_u, uniques = pd.factorize(labels)
    remap = np.zeros(max(len(uniques), 1), dtype=np.int32)
    for i, u in enumerate(uniques):
        remap[i] = raw_map.get(u, 0)     # "" (and any oddity) -> unassigned
    codes = np.where(codes_u < 0, 0, remap[np.maximum(codes_u, 0)]).astype(np.int32)
    return codes, mapping


def border_reassign(xy: np.ndarray, labels: np.ndarray, fraction: float) -> tuple[np.ndarray, dict]:
    """Reassign the most border-adjacent ``fraction`` of molecules.

    ``xy`` is (n, 2) coordinates, ``labels`` cell-id strings ("" = background).
    Returns the new label array and a stats dict.
    """
    from scipy.spatial import cKDTree

    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"fraction must be in (0, 1], got {fraction}")
    assigned = labels != ""
    n_assigned = int(assigned.sum())
    if n_assigned == 0:
        return labels.copy(), {"n_reassigned": 0, "n_assigned": 0, "max_cross_cell_distance_um": None}

    idx_assigned = np.flatnonzero(assigned)
    coords = xy[idx_assigned]
    local_labels = labels[idx_assigned]

    k = min(10, n_assigned)
    tree = cKDTree(coords)
    dists, neigh = tree.query(coords, k=k)
    if k == 1:
        dists = dists[:, None]
        neigh = neigh[:, None]

    # For each molecule, the nearest neighbour owned by a different cell.
    # Vectorized over all molecules: `neigh` is distance-ordered per row, so
    # the first foreign label per row is the nearest foreign molecule (the
    # original per-molecule loop did exactly this scan).
    row_labels = local_labels[neigh]                      # (n, k)
    foreign = row_labels != local_labels[:, None]         # (n, k)
    has_foreign = foreign.any(axis=1)
    first = np.argmax(foreign, axis=1)
    rows = np.flatnonzero(has_foreign)
    cols = first[rows]
    best_d = np.full(n_assigned, np.inf)
    best_cell = np.full(n_assigned, "", dtype=object)
    best_d[rows] = dists[rows, cols]
    best_cell[rows] = row_labels[rows, cols]

    n_reassign = int(round(fraction * n_assigned))
    order = np.lexsort((idx_assigned, best_d))  # distance asc, then row index
    finite = order[np.isfinite(best_d[order])]
    take = finite[:n_reassign]

    out = labels.copy()
    for i in take:
        out[idx_assigned[i]] = best_cell[i]
    reassigned = np.zeros(len(labels), dtype=bool)
    reassigned[idx_assigned[take]] = True
    changed = reassigned & (out != labels)
    cut_d = float(np.max(best_d[take])) if len(take) else None
    stats = {
        "mode": "border",
        "fraction_requested": fraction,
        "n_assigned": n_assigned,
        "n_reassigned": int(changed.sum()),
        "n_unassigned": int((labels == "").sum()),
        "max_reassigned_cross_cell_distance_um": cut_d,
    }
    return out, stats


def dilate_cells(xy: np.ndarray, labels: np.ndarray, distance: float = 2.0) -> tuple[np.ndarray, dict]:
    """Dilate every cell by ``distance`` um, claiming neighbours and background."""
    import shapely

    if distance <= 0:
        raise ValueError(f"distance must be > 0, got {distance}")
    cell_values = sorted(pd.unique(labels[labels != ""]).tolist())
    cell_index = {c: i for i, c in enumerate(cell_values)}
    hulls = []
    if cell_values:
        # group row indices per label in one pass (the per-label boolean scan
        # over all molecules is what made this O(n_cells * n_molecules))
        codes, uniques = pd.factorize(labels)
        uniques = np.asarray(uniques, dtype=object)
        inv = np.argsort(codes, kind="stable")
        srt = codes[inv]
        starts = np.searchsorted(srt, np.arange(len(uniques)), side="left")
        ends = np.searchsorted(srt, np.arange(len(uniques)), side="right")
        grouped = {str(u): inv[s:e]
                   for u, s, e in zip(uniques, starts, ends)}
        hulls_arr = np.full(len(cell_values), None, dtype=object)
        for cell in cell_values:
            idx = grouped.get(str(cell))
            if idx is not None and len(idx):
                hulls_arr[cell_index[cell]] = shapely.MultiPoint(
                    xy[idx]).convex_hull
        hulls = list(hulls_arr)
    if not hulls:
        return labels.copy(), {"mode": "dilate", "n_reassigned": 0, "n_absorbed_background": 0}

    tree = shapely.STRtree(hulls)
    points = shapely.points(xy)
    # STRtree.query returns (input_geometry_indices, tree_geometry_indices).
    pt_idx, geom_idx = tree.query(points, predicate="dwithin", distance=float(distance))
    hull_arr = np.array(hulls, dtype=object)
    d = shapely.distance(points[pt_idx], hull_arr[geom_idx])
    # own[p] = cell index of the molecule's own cell (-1 = background/unknown)
    own_codes = np.array([-1 if str(lab) not in cell_index else cell_index[str(lab)]
                          for lab in labels], dtype=int)

    # Best *other* cell per molecule: smallest distance to its hull, ties
    # broken by cell index. Excluding the molecule's own cell is what makes
    # this a dilation: the band within `distance` inside a cell's border is
    # claimed by the neighbouring cell whose hull reaches it first.
    other = geom_idx != own_codes[pt_idx]
    p, g, dd = pt_idx[other], geom_idx[other], d[other]
    if len(p):
        # min per molecule by (distance, cell index), fully vectorized
        order = np.lexsort((g, dd, p))       # primary p, then dd, then g
        p_sorted = p[order]
        first = np.empty(len(order), dtype=bool)
        first[0] = True
        first[1:] = p_sorted[1:] != p_sorted[:-1]
        sel = order[first]
        best_p, best_g = p[sel], g[sel]
    else:
        best_p = np.empty(0, dtype=np.int64)
        best_g = np.empty(0, dtype=np.int64)

    out = labels.copy()
    n_flip = 0
    n_absorb = 0
    if len(best_p):
        target = np.asarray(cell_values, dtype=object)[best_g]
        current = out[best_p]
        absorb_mask = current == ""
        flip_mask = (~absorb_mask) & (current != target)
        n_absorb = int(absorb_mask.sum())
        n_flip = int(flip_mask.sum())
        out[best_p] = target
    stats = {
        "mode": "dilate",
        "distance_um": distance,
        "n_cells": len(cell_values),
        "n_reassigned_from_neighbours": n_flip,
        "n_absorbed_background": n_absorb,
        "n_assigned": int((labels != "").sum()),
        "n_unassigned": int((labels == "").sum()),
    }
    return out, stats


def write_assignment(path: Path, labels: np.ndarray, mapping: dict[str, int] | None = None) -> None:
    """Write the ``cell`` (int) + ``cell_label`` (str) assignment parquet."""
    if mapping is None:
        codes, mapping = cell_codes(labels)
    else:
        codes = np.zeros(len(labels), dtype=np.int32)
        for label, code in mapping.items():
            codes[labels == label] = code
    pd.DataFrame({"cell": codes, "cell_label": labels.astype(str)}).to_parquet(path, index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--molecules", type=Path, required=True)
    parser.add_argument("--cell-column", default="cell_vendor")
    parser.add_argument("--mode", choices=["border", "dilate"], required=True)
    parser.add_argument("--fraction", type=float, default=0.1,
                        help="border mode: fraction of assigned molecules reassigned")
    parser.add_argument("--distance", type=float, default=2.0,
                        help="dilate mode: dilation distance in um")
    parser.add_argument("--out", type=Path, required=True, help="output assignment parquet")
    parser.add_argument("--report", type=Path, help="output JSON stats")
    args = parser.parse_args(argv)

    frame = load_assignment_table(args.molecules, args.cell_column)
    xy = frame[["x", "y"]].to_numpy(dtype=float)
    labels = frame["cell_vendor"].fillna("").astype(str).to_numpy()

    if args.mode == "border":
        new_labels, stats = border_reassign(xy, labels, args.fraction)
    else:
        new_labels, stats = dilate_cells(xy, labels, args.distance)
    stats.update({
        "molecules": str(args.molecules),
        "cell_column": args.cell_column,
        "n_molecules": len(labels),
        "n_changed": int((new_labels != labels).sum()),
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_assignment(args.out, new_labels)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with open(args.report, "w") as fh:
            json.dump(stats, fh, indent=2)
            fh.write("\n")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
