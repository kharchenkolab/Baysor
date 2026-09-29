#!/usr/bin/env python
"""Transfer baseline cell types to a target segmentation by molecule overlap.

Cell typing for the cellAdmix admixture audit depends on the segmentation:
types assigned on the baseline (e.g. vendor) segmentation cannot be reused
verbatim for a different segmentation because the cell ids differ. The
recommended comparability-preserving scheme is to type the baseline once and
then, for every cell of the target segmentation, take the majority baseline
type among the molecules that the target cell contains. Cells whose molecules
carry no baseline type are dropped from the annotation (audit.py excludes
their molecules and records the filter).

Usage:
    python transfer.py --molecules molecules.parquet \
        --baseline-cell-column cell_vendor \
        --target-assignment baysor.parquet \
        --baseline-celltypes vendor_celltypes.parquet \
        --out baysor_celltypes.parquet --report transfer.json

When the baseline segmentation lives in a harness assignment table instead of
a molecules column, give ``--baseline-assignment`` (row-aligned per-molecule
``cell`` ints, 0 = unassigned) instead of ``--baseline-cell-column``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def read_celltypes(path: Path | str, cell_col: str | None = None, type_col: str | None = None) -> pd.Series:
    """Read a cell -> type table (parquet/CSV) as a str-indexed Series."""
    path = Path(path)
    frame = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    if cell_col is None:
        for cand in ("cell", "cell_id", "cell_vendor"):
            if cand in frame.columns:
                cell_col = cand
                break
        else:
            cell_col = frame.columns[0]
    if type_col is None:
        for cand in ("celltype", "cell_type", "type", "cluster", "label", "annotation"):
            if cand in frame.columns and cand != cell_col:
                type_col = cand
                break
        else:
            type_col = [c for c in frame.columns if c != cell_col][0]
    out = pd.Series(
        frame[type_col].astype(str).to_numpy(),
        index=frame[cell_col].astype(str).to_numpy(),
        name="celltype",
    )
    out = out[(out != "") & (out != "nan")]
    return out[~out.index.duplicated(keep="first")]


def read_target_cells(molecules_path: Path | str, *, assignment: Path | str | None,
                      cell_column: str | None, assignment_col: str = "cell") -> np.ndarray:
    """Per-molecule target cell labels (str, '' = unassigned), row-aligned."""
    if (assignment is None) == (cell_column is None):
        raise ValueError("exactly one of assignment / cell_column is required")
    if assignment is not None:
        import pyarrow.parquet as pq

        codes = pd.read_parquet(assignment, columns=[assignment_col])[assignment_col].to_numpy()
        schema = pq.read_schema(assignment)
        labels_col = None
        if "cell_label" in schema.names:
            labels_col = pd.read_parquet(assignment, columns=["cell_label"])["cell_label"].to_numpy()
        n_mol = len(pd.read_parquet(molecules_path, columns=["x"]))
        if len(codes) != n_mol:
            raise ValueError(f"assignment has {len(codes)} rows but molecules have {n_mol}")
        if labels_col is not None:
            out = labels_col.astype(str)
        else:
            out = np.where(codes > 0, codes.astype(str), "")
        return out
    col = pd.read_parquet(molecules_path, columns=[cell_column])[cell_column]
    return col.fillna("").astype(str).to_numpy()


def transfer_celltypes(
    target_cells: np.ndarray,
    baseline_cells: np.ndarray,
    baseline_types: pd.Series,
    *,
    min_votes: int = 1,
) -> tuple[pd.Series, dict]:
    """Majority baseline type per target cell, voted by shared molecules.

    Returns a Series indexed by target cell id and a stats dict. Target cells
    with fewer than ``min_votes`` molecules carrying a baseline type are absent
    from the result.
    """
    target_cells = np.asarray(target_cells, dtype=object)
    baseline_cells = np.asarray(baseline_cells, dtype=object)
    if len(target_cells) != len(baseline_cells):
        raise ValueError("target_cells and baseline_cells must be row-aligned")
    typed = baseline_cells != ""
    mapped = np.array([
        baseline_types.get(str(c), "") if c != "" else ""
        for c in baseline_cells
    ], dtype=object)
    usable = (target_cells != "") & (mapped != "")
    table = pd.crosstab(
        pd.Series(target_cells[usable], dtype=object),
        pd.Series(mapped[usable], dtype=object),
    )
    if table.empty:
        target_unique = set(target_cells[target_cells != ""])
        return pd.Series(dtype=object, name="celltype"), {
            "n_target_cells": len(target_unique),
            "n_typed_target_cells": 0,
            "n_dropped_untyped_target_cells": len(target_unique),
            "n_molecules_with_baseline_type": 0,
            "n_molecules_baseline_unassigned_or_untyped": int(((target_cells != "") & (mapped == "")).sum()),
            "min_votes": min_votes,
            "baseline_types": sorted(set(baseline_types.astype(str))),
        }
    counts = table.to_numpy(dtype=np.int64)
    types = table.columns.to_numpy()
    # Majority; ties broken by lexicographically smallest type (stable order).
    order = np.argsort(types, kind="stable")
    counts_sorted = counts[:, order]
    types_sorted = types[order]
    best = counts_sorted.argmax(axis=1)
    winners = types_sorted[best]
    votes = counts_sorted.max(axis=1)
    cell_ids = table.index.to_numpy().astype(str)
    keep = votes >= min_votes
    result = pd.Series(winners[keep], index=cell_ids[keep], name="celltype")

    target_unique = set(target_cells[target_cells != ""])
    stats = {
        "n_target_cells": len(target_unique),
        "n_typed_target_cells": int(keep.sum()),
        "n_dropped_untyped_target_cells": len(target_unique) - int(keep.sum()),
        "n_molecules_with_baseline_type": int(usable.sum()),
        "n_molecules_baseline_unassigned_or_untyped": int(((target_cells != "") & ~usable).sum()),
        "min_votes": min_votes,
        "baseline_types": sorted(set(baseline_types.astype(str))),
    }
    return result, stats


def write_celltypes(path: Path, types: pd.Series) -> None:
    out = types.rename_axis("cell").reset_index()
    out.columns = ["cell", "celltype"]
    out.to_parquet(path, index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--molecules", type=Path, required=True)
    parser.add_argument("--baseline-cell-column", default="cell_vendor")
    parser.add_argument("--baseline-assignment", type=Path,
                        help="baseline segmentation assignment parquet (row-aligned "
                             "int 'cell' column); alternative to --baseline-cell-column")
    parser.add_argument("--target-assignment", type=Path,
                        help="target segmentation assignment parquet (int 'cell' column)")
    parser.add_argument("--target-cell-column",
                        help="target cell column in the molecules table (alternative to --target-assignment)")
    parser.add_argument("--baseline-celltypes", type=Path, required=True)
    parser.add_argument("--celltype-cell-col")
    parser.add_argument("--celltype-type-col")
    parser.add_argument("--min-votes", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    if args.baseline_assignment is not None and args.baseline_cell_column != "cell_vendor":
        parser.error("--baseline-assignment and --baseline-cell-column are mutually exclusive")

    if args.baseline_assignment is not None:
        baseline_cells = read_target_cells(
            args.molecules, assignment=args.baseline_assignment, cell_column=None)
    else:
        baseline_cells = pd.read_parquet(
            args.molecules, columns=[args.baseline_cell_column]
        )[args.baseline_cell_column].fillna("").astype(str).to_numpy()
    baseline_types = read_celltypes(args.baseline_celltypes,
                                    cell_col=args.celltype_cell_col,
                                    type_col=args.celltype_type_col)
    target_cells = read_target_cells(
        args.molecules,
        assignment=args.target_assignment,
        cell_column=args.target_cell_column,
    )
    types, stats = transfer_celltypes(target_cells, baseline_cells, baseline_types,
                                      min_votes=args.min_votes)
    stats.update({
        "molecules": str(args.molecules),
        "baseline_cell_column": args.baseline_cell_column if args.baseline_assignment is None else None,
        "baseline_assignment": str(args.baseline_assignment) if args.baseline_assignment is not None else None,
        "baseline_celltypes": str(args.baseline_celltypes),
        "target": str(args.target_assignment or args.target_cell_column),
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_celltypes(args.out, types)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with open(args.report, "w") as fh:
            json.dump(stats, fh, indent=2)
            fh.write("\n")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
