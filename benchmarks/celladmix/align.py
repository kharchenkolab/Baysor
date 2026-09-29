#!/usr/bin/env python
"""Align a segmentation's per-molecule labels onto the dataset molecule rows.

Reads a segmentation output molecule table (e.g. Baysor's
``molecules.parquet`` with a string ``cell`` column), verifies that its rows
are in the same order as the contract-format dataset molecules, and writes the
``{cell: int32, cell_label: str}`` assignment parquet consumed by
``audit.py --assignment`` and ``transfer.py --target-assignment``.

Usage:
    python align.py --molecules molecules.parquet \
        --segmentation baysor_seg/molecules.parquet --cell-column cell \
        --out baysor.assignment.parquet --report baysor.align.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit import UNASSIGNED_TOKENS
from degrade import cell_codes, write_assignment


def align_segmentation(
    molecules_path: Path | str,
    segmentation_path: Path | str,
    cell_column: str,
    *,
    extra_unassigned: set[str] | None = None,
    atol: float = 1e-9,
) -> tuple[np.ndarray, dict]:
    """Return per-row cell labels ('' = unassigned) plus alignment stats.

    Raises ``ValueError`` when the segmentation table is not row-aligned with
    the dataset molecules (different length or different x/y/gene content).
    """
    base = pd.read_parquet(molecules_path, columns=["x", "y", "gene"])
    seg = pd.read_parquet(segmentation_path, columns=[cell_column])
    if len(base) != len(seg):
        raise ValueError(
            f"segmentation has {len(seg)} rows but dataset has {len(base)}; "
            "the segmentation must be in the same row order")
    # If the segmentation table carries coordinates, verify them; a bare label
    # column of equal length is trusted (that is the documented contract).
    import pyarrow.parquet as pq

    seg_columns = set(pq.read_schema(segmentation_path).names)
    if {"x", "y", "gene"} <= seg_columns:
        seg_full = pd.read_parquet(segmentation_path, columns=["x", "y", "gene"])
        same = (
            np.allclose(seg_full["x"].to_numpy(dtype=float),
                        base["x"].to_numpy(dtype=float), atol=atol, rtol=0)
            and np.allclose(seg_full["y"].to_numpy(dtype=float),
                            base["y"].to_numpy(dtype=float), atol=atol, rtol=0)
            and (seg_full["gene"].astype(str).to_numpy()
                 == base["gene"].astype(str).to_numpy()).all()
        )
        if not same:
            raise ValueError("segmentation x/y/gene do not match the dataset rows")

    labels = seg[cell_column].fillna("").astype(str).to_numpy()
    unassigned = set(UNASSIGNED_TOKENS) | set(extra_unassigned or set())
    labels = np.where(np.isin(labels, list(unassigned)), "", labels)
    stats = {
        "n_rows": int(len(labels)),
        "n_unassigned": int((labels == "").sum()),
        "n_cells": int(len(set(labels[labels != ""]))),
        "unassigned_labels": sorted(unassigned),
        "segmentation": str(segmentation_path),
        "cell_column": cell_column,
    }
    return labels, stats


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--molecules", type=Path, required=True,
                        help="dataset molecules parquet (row reference)")
    parser.add_argument("--segmentation", type=Path, required=True,
                        help="segmentation molecule table with per-molecule labels")
    parser.add_argument("--cell-column", default="cell")
    parser.add_argument("--unassigned", default="",
                        help="extra comma-separated label values to treat as unassigned")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)

    extra = {v for v in args.unassigned.split(",") if v}
    try:
        labels, stats = align_segmentation(
            args.molecules, args.segmentation, args.cell_column, extra_unassigned=extra)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_assignment(args.out, labels)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with open(args.report, "w") as fh:
            json.dump(stats, fh, indent=2)
            fh.write("\n")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
