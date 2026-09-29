#!/usr/bin/env python
"""Fetch the 10x Xenium pancreas FFPE bundle members and cut benchmark crops.

The cellAdmix audit validation needs a real Xenium dataset with a vendor
segmentation. This script downloads only the two bundle members that matter
(``transcripts.parquet`` and ``cells.parquet``, ~155 MB total) from the 10x
public zip with HTTP range requests (``remotezip``), then writes contract-format
dataset directories (``molecules.parquet`` + ``meta.json``) cropped to a fixed
square window.

Crop selection is deterministic: the window is the densest ``--side-um`` square
on a coarse grid over the vendor cell centroids (ties broken by lowest x, then
lowest y); if the window contains more than ``--max-molecules`` molecules after
gene/qv filtering, the side length is shrunk in ``--shrink-step-um`` increments
until the budget fits.

Filters applied (recorded in ``meta.json``):
  * control probes and blank codewords removed (feature names starting with
    NegControl, UnassignedCodeword, DeprecatedCodeword, Intergenic_Region);
  * ``qv < --qv-min`` molecules dropped;
  * rows outside the crop bbox dropped;
  * rows sorted by (y, x) as required by the dataset contract.

Usage:
    python fetch_pancreas.py --data-root $BAYSOR_BENCH_DATA \
        --crop-id pancreas_crop_quick --side-um 625 --max-molecules 150000
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ZIP_URL = (
    "https://cf.10xgenomics.com/samples/xenium/2.0.0/"
    "Xenium_V1_human_Pancreas_FFPE/Xenium_V1_human_Pancreas_FFPE_outs.zip"
)
DATASET_PAGE = (
    "https://www.10xgenomics.com/datasets/"
    "ffpe-human-pancreas-with-xenium-multimodal-cell-segmentation-1-standard"
)
CONTROL_PREFIXES = ("NegControl", "UnassignedCodeword", "DeprecatedCodeword", "Intergenic_Region")
MEMBERS = ("transcripts.parquet", "cells.parquet")
GRID_STEP_UM = 25.0  # grid step for the dense-window search


def default_data_root() -> Path:
    return Path(os.environ.get("BAYSOR_BENCH_DATA", "")) or Path(__file__).resolve().parents[2] / ".bench-data"


def fetch_members(raw_dir: Path, *, url: str = ZIP_URL, members: tuple[str, ...] = MEMBERS) -> dict[str, Path]:
    """Download the requested zip members with range requests, if missing."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    missing = [m for m in members if not (raw_dir / m).exists()]
    if missing:
        import remotezip

        print(f"fetching {len(missing)} member(s) from {url}", file=sys.stderr)
        with remotezip.RemoteZip(url) as zf:
            for name in missing:
                print(f"  {name} ...", file=sys.stderr, flush=True)
                zf.extract(name, path=raw_dir)
    for name in members:
        path = raw_dir / name
        if not path.exists():
            raise FileNotFoundError(f"missing zip member after fetch: {path}")
        out[name] = path
    return out


def pick_window(cell_xy: np.ndarray, side: float, grid_step: float = GRID_STEP_UM) -> tuple[float, float, float, float]:
    """Return (x0, y0, x1, y1) of the densest side x side window.

    Deterministic: candidates start at multiples of ``grid_step`` (the side
    length is rounded to a multiple of the grid step); ties are broken by
    lowest x, then lowest y. Uses a summed-area table over a coarse cell-count
    histogram, so the search is O(cells + grid).
    """
    if not len(cell_xy):
        raise ValueError("no cells to place a window")
    k = max(1, int(round(side / grid_step)))
    side = k * grid_step
    xmin, ymin = cell_xy.min(axis=0)
    xmax, ymax = cell_xy.max(axis=0)
    x0_start = np.floor(xmin / grid_step) * grid_step
    y0_start = np.floor(ymin / grid_step) * grid_step
    # Candidate starts run up to the largest cell so every cell is coverable;
    # pad by k so windows wider than the data extent stay valid (start at 0).
    nx = int(np.floor((xmax - x0_start) / grid_step)) + 1
    ny = int(np.floor((ymax - y0_start) / grid_step)) + 1
    ix = np.clip(((cell_xy[:, 0] - x0_start) / grid_step).astype(int), 0, nx - 1)
    iy = np.clip(((cell_xy[:, 1] - y0_start) / grid_step).astype(int), 0, ny - 1)
    hist = np.zeros((nx + k + 1, ny + k + 1), dtype=np.int64)
    np.add.at(hist, (ix + 1, iy + 1), 1)
    integral = hist.cumsum(axis=0).cumsum(axis=1)
    # Inclusive summed-area table: window over cell-bins [i, i+k) covers hist
    # rows [i+1, i+k], so its sum is I[i+k, j+k] - I[i, j+k] - I[i+k, j]
    # + I[i, j]. The k-bin padding keeps windows valid even when the side is
    # wider than the data extent (start i = 0 then covers every cell).
    counts = (
        integral[k:nx + k, k:ny + k]
        - integral[0:nx, k:ny + k]
        - integral[k:nx + k, 0:ny]
        + integral[0:nx, 0:ny]
    )
    # argmax over C-order picks the lowest i, then lowest j on ties.
    i, j = np.unravel_index(np.argmax(counts), counts.shape)
    x0 = x0_start + i * grid_step
    y0 = y0_start + j * grid_step
    return (float(x0), float(y0), float(x0 + side), float(y0 + side))


def build_crop(
    transcripts_path: Path,
    cells_path: Path,
    out_dir: Path,
    *,
    crop_id: str,
    side_um: float,
    max_molecules: int,
    qv_min: float,
    shrink_step_um: float = 25.0,
    min_side_um: float = 100.0,
) -> dict:
    """Cut one crop and write ``molecules.parquet`` + ``meta.json`` to ``out_dir``."""
    cells = pd.read_parquet(cells_path, columns=["cell_id", "x_centroid", "y_centroid"])
    cell_xy = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=float)

    tx = pd.read_parquet(
        transcripts_path,
        columns=["feature_name", "x_location", "y_location", "qv", "cell_id"],
    )
    gene_mask = ~tx["feature_name"].astype(str).str.startswith(CONTROL_PREFIXES)
    tx = tx[gene_mask & (tx["qv"].astype(float) >= qv_min)]
    n_control = int((~gene_mask).sum())

    side = float(side_um)
    while True:
        x0, y0, x1, y1 = pick_window(cell_xy, side)
        in_window = (
            (tx["x_location"] >= x0) & (tx["x_location"] < x1)
            & (tx["y_location"] >= y0) & (tx["y_location"] < y1)
        )
        n_mol = int(in_window.sum())
        if n_mol <= max_molecules or side <= min_side_um:
            break
        # Shrink toward the window center while staying on the same grid.
        side = max(min_side_um, side - shrink_step_um)
        print(f"  side={side:.0f} um -> {n_mol} molecules exceeds budget {max_molecules}; shrinking",
              file=sys.stderr)

    crop = tx[in_window].copy()
    crop = crop.rename(columns={
        "feature_name": "gene", "x_location": "x", "y_location": "y",
        "cell_id": "cell_vendor",
    })
    crop["cell_vendor"] = crop["cell_vendor"].fillna("").astype(str)
    crop["x"] = crop["x"].astype("float64")
    crop["y"] = crop["y"].astype("float64")
    crop["gene"] = crop["gene"].astype(str)
    crop["qv"] = crop["qv"].astype("float32")
    crop = crop[["x", "y", "gene", "qv", "cell_vendor"]]
    crop = crop.sort_values(["y", "x"], kind="mergesort").reset_index(drop=True)

    n_vendor_cells = int((crop["cell_vendor"] != "").sum())
    n_unique_vendor = int(crop.loc[crop["cell_vendor"] != "", "cell_vendor"].nunique())
    n_genes = int(crop["gene"].nunique())
    area = (x1 - x0) * (y1 - y0) / 1e6  # mm^2
    assigned = crop[crop["cell_vendor"] != ""]

    out_dir.mkdir(parents=True, exist_ok=True)
    crop.to_parquet(out_dir / "molecules.parquet", index=False)

    meta = {
        "id": crop_id,
        "kind": "real",
        "tier": "quick" if len(crop) <= 150_000 else "full",
        "platform": "Xenium",
        "location": "cache",  # benchmark cache dataset, not a suite real/ dataset
        "source": {
            "url": ZIP_URL,
            "dataset_page": DATASET_PAGE,
            "license": "10x Genomics public dataset terms",
            "original_dataset": "Xenium V1 human pancreas FFPE (multimodal cell segmentation)",
            "retrieved": pd.Timestamp.utcnow().date().isoformat(),
        },
        "crop": {
            "bbox_um": [round(float(x0), 3), round(float(y0), 3), round(float(x1), 3), round(float(y1), 3)],
            "z_range_um": None,
            "note": "densest vendor-cell window; shrunk to fit the molecule budget",
        },
        "stats": {
            "n_molecules": int(len(crop)),
            "n_genes": n_genes,
            "area_um2": round((x1 - x0) * (y1 - y0), 1),
            "molecules_per_um2": round(len(crop) / ((x1 - x0) * (y1 - y0)), 4),
            "n_vendor_cells": n_unique_vendor,
            "vendor_cells_per_mm2": round(n_unique_vendor / area, 1),
            "n_molecules_assigned_vendor": n_vendor_cells,
        },
        "difficulty": {
            "cell_density": "dense" if n_unique_vendor / area > 7000 else ("medium" if n_unique_vendor / area > 2500 else "sparse"),
            "gene_panel": ("tiny" if n_genes < 50 else "small" if n_genes < 250 else
                           "medium" if n_genes < 700 else "large" if n_genes < 2000 else "huge"),
            "notes": "cropped from the full Xenium pancreas FFPE bundle",
        },
        "filters": {
            "control_features_removed": n_control,
            "qv_min": qv_min,
            "control_feature_prefixes": list(CONTROL_PREFIXES),
        },
        "baysor": {
            "scale_um": 5.0,
            "scale_std": "25%",
            "min_molecules_per_cell": 20,
            "prior": "none",
            "prior_confidence": 0.5,
            # Contract-format tables already use x/y/gene columns, so no
            # configs/xenium.toml (that config maps raw Xenium column names).
            "config": None,
            "extra_args": ["--force-2d"],
        },
        "images": [],
        "truth": None,
    }
    with open(out_dir / "meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")
    return meta


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path, default=default_data_root(),
                        help="benchmark data root (default: $BAYSOR_BENCH_DATA)")
    parser.add_argument("--crop-id", default="pancreas_crop_quick", help="output dataset id")
    parser.add_argument("--side-um", type=float, default=625.0, help="crop window side length in um")
    parser.add_argument("--max-molecules", type=int, default=150_000, help="molecule budget for the crop")
    parser.add_argument("--qv-min", type=float, default=20.0, help="minimum qv kept")
    parser.add_argument("--shrink-step-um", type=float, default=25.0)
    parser.add_argument("--download-only", action="store_true", help="only fetch the raw zip members")
    args = parser.parse_args(argv)

    cache = args.data_root / "cache" / "celladmix"
    raw_dir = cache / "raw_pancreas"
    members = fetch_members(raw_dir)
    if args.download_only:
        return 0

    out_dir = cache / "datasets" / args.crop_id
    meta = build_crop(
        members["transcripts.parquet"], members["cells.parquet"], out_dir,
        crop_id=args.crop_id, side_um=args.side_um, max_molecules=args.max_molecules,
        qv_min=args.qv_min, shrink_step_um=args.shrink_step_um,
    )
    stats = meta["stats"]
    print(json.dumps({
        "out_dir": str(out_dir),
        "n_molecules": stats["n_molecules"],
        "n_genes": stats["n_genes"],
        "n_vendor_cells": stats["n_vendor_cells"],
        "bbox_um": meta["crop"]["bbox_um"],
        "tier": meta["tier"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
