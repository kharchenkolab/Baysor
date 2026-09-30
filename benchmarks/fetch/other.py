#!/usr/bin/env python3
"""BENCH-REALO: build cropped real benchmark datasets from non-Xenium platforms.

Reads ``benchmarks/datasets/real_other.yaml`` and materialises one directory
per dataset under ``$BAYSOR_BENCH_DATA/real/<id>/`` following the dataset
contract in ``benchmarks/README.md``:

    molecules.parquet   sorted by (y, x), coordinates in um, blanks removed
    meta.json           provenance, stats, difficulty classes, Baysor params
    images/*.tif        cropped stains (where the source has them)
    reference/*         vendor cell labels / polygons for the crop
    README.md           dataset-specific provenance notes

Sub-commands:

    build    download sources and (re)build datasets
    smoke    run the Release Baysor binary on every quick crop, recording
             wall time and peak RSS into benchmarks/baselines/
    report   print a markdown table with per-dataset stats

Sources (all downloadable without login; see the inventory in
``benchmarks/datasets/real_other_inventory.md``):

    ISS       pklab.med.harvard.edu (example files; fetched via headless
              Firefox because the host challenges non-browser clients) +
              figshare DAPI image and pciSeq spot-to-cell assignments
    osmFISH   pklab molecule csv + vitessce polygon mirror of polyT_seg
    STARmap   pklab molecules.csv + segmentation.tiff
    ileum     Dryad doi:10.5061/dryad.jm63xsjb2 (Anubis-protected, ranged zip)
    CosMx     nanostring-public-share S3 (NSCLC Lung5_Rep1) and
              objects.liquidweb.services smi-public bucket (WTX colon)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))

import other_utils as U  # noqa: E402

TODAY = date.today().isoformat()

# Performance finding (BENCH-REALO): the default mrf/ICA initialisation does not
# scale to the ~19k-gene CosMx WTX crops; recorded for the report, the dataset
# READMEs and meta.json difficulty.notes, and tagged for a future regression test.
WTX_PERF_NOTE = (
    "perf-stress: default ICA init — with the default --cluster-method mrf this "
    "~19k-gene crop never finished: >70 min stuck at 'Clustering molecules into 4 "
    "types (ICA init)' at 9.8 GB RSS (run killed by the coordinator; a 20-min capped "
    "rerun is recorded as a timeout in benchmarks/baselines/real_other_smoke.json). "
    "meta.baysor therefore uses --cluster-method louvain (113 s on the quick crop; "
    "leiden and none complete in ~111 s, all on 6 threads)."
)

# Documented pixel sizes (see real_other_inventory.md for sources).
ISS_PIXEL_UM = 1.0 / 3.0            # Qian et al.: 3 px = 1 um (top-hat filter radii)
OSMFISH_PIXEL_UM = 0.065            # linnarssonlab.org/osmFISH availability page
STARMAP_FRAME_PX = (17219, 3734)    # segmentation.tiff dimensions (x, y)
STARMAP_FRAME_UM = (1400.0, 300.0)  # Wang et al. Fig. 2B: "Full field: 1.4 by 0.3 mm"
STARMAP_F_UM = STARMAP_FRAME_UM[0] / STARMAP_FRAME_PX[0]
STARMAP_Z_PLANES = 36               # goodPoints z: 0..36
STARMAP_VOLUME_UM = 8.0             # Wang et al.: "Eight-um-thick volumes ... up to 1000 cells"
STARMAP_Z_STEP_UM = STARMAP_VOLUME_UM / STARMAP_Z_PLANES
ILEUM_ROOT = "data_release_baysor_merfish_gut"
NSCLC_PIXEL_UM = 0.18               # SMI Data File ReadMe (S3): "multiply by 0.18 um per pixel"
WTX_PIXEL_UM = 0.12028              # WTX README: "0.12028 um per pixel"

PKLAB = "http://pklab.med.harvard.edu/viktor/baysor"
VITESSE_MOLECULES = (
    "https://data-1.vitessce.io/0.0.33/main/codeluppi-2018/"
    "codeluppi_2018_nature_methods.molecules.csv"
)
VITESSE_SEGMENTS = (
    "https://data-1.vitessce.io/0.0.33/main/codeluppi-2018/"
    "codeluppi_2018_nature_methods.cells.segmentations.json"
)
DRYAD_FILE_URL = "https://datadryad.org/stash/downloads/file_stream/1019818"
NSCLC_BUCKET = "https://nanostring-public-share.s3.us-west-2.amazonaws.com"
WTX_BUCKET = "https://objects.liquidweb.services/smi-public/wtx_manuscript"

# Expected byte sizes of cached source files (guards against truncation).
SOURCE_FILES = {
    "iss_molecules": (f"{PKLAB}/iss/pciSeq_3-3.csv", 65058308),
    "iss_mask": (f"{PKLAB}/iss/DAPI_3-3_mask.tif", 653662138),
    "osmfish_molecules": (f"{PKLAB}/osm_fish/mRNA_coords_raw_counting.csv", 41573664),
    "starmap_molecules": (f"{PKLAB}/starmap/molecules.csv", 17985433),
    "starmap_segmentation": (f"{PKLAB}/starmap/segmentation.tiff", 128591750),
    "osmfish_segments": (VITESSE_SEGMENTS, 648485),
    "iss_dapi_jpg": ("https://ndownloader.figshare.com/files/13160405", 36081355),
    "iss_assign_left": ("https://ndownloader.figshare.com/files/18772673", 443716),
    "iss_assign_right": ("https://ndownloader.figshare.com/files/18772718", 574330),
    "nsclc_metadata": (f"{NSCLC_BUCKET}/SMI/Lung5_Rep1/Lung5_Rep1_metadata_file.csv", 11512298),
    "nsclc_fov": (f"{NSCLC_BUCKET}/SMI/Lung5_Rep1/Lung5_Rep1_fov_positions_file.csv", 1195),
    "nsclc_tx": (f"{NSCLC_BUCKET}/SMI/Lung5_Rep1/Lung5_Rep1_tx_file.csv", 3403685000),
    "wtx_metadata": (f"{WTX_BUCKET}/colon_discovery/S0_metadata_file.csv.gz", 39073339),
    "wtx_fov": (f"{WTX_BUCKET}/colon_discovery/S0_fov_positions_file.csv.gz", 4044),
    "wtx_readme": (f"{WTX_BUCKET}/colon_discovery/README_coloncancer.html", 690228),
    "wtx_tx": (f"{WTX_BUCKET}/colon_discovery/S0_tx_file.csv.gz", 11397927299),
}
DRYAD_MEMBERS = {
    "molecules": f"{ILEUM_ROOT}/raw_data/molecules.csv",
    "dapi": f"{ILEUM_ROOT}/raw_data/dapi_stack.tif",
    "membrane": f"{ILEUM_ROOT}/raw_data/membrane_stack.tif",
    "cellpose": f"{ILEUM_ROOT}/data_analysis/cellpose/cell_boundaries/results/cellpose_dapi.tif",
    "cellpose_membrane": f"{ILEUM_ROOT}/data_analysis/cellpose/cell_boundaries/results/cellpose_membrane.tif",
    "readme": f"{ILEUM_ROOT}/README.txt",
    "columns": f"{ILEUM_ROOT}/file_organization/molecules.txt",
}


# ---------------------------------------------------------------------------
# Manifest / dataset directory helpers
# ---------------------------------------------------------------------------

def load_manifest(path: Path | None = None) -> dict:
    path = path or (U.repo_root() / "benchmarks" / "datasets" / "real_other.yaml")
    with open(path) as fh:
        manifest = yaml.safe_load(fh)
    ids = [d["id"] for d in manifest["datasets"]]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate dataset ids in manifest")
    return manifest


def dataset_dir(dataset_id: str) -> Path:
    return U.real_dir() / dataset_id


def ensure_sources(keys: Sequence[str], *, browser_keys: Sequence[str] = ()) -> dict[str, Path]:
    """Download the manifest-listed source files into the cache (idempotent)."""
    out: dict[str, Path] = {}
    for key in keys:
        url, size = SOURCE_FILES[key]
        dest = U.cache_dir("real_other") / Path(url).name
        if key in browser_keys:
            out[key] = U.browser_download(url, dest)
        else:
            out[key] = U.ensure_file(url, dest, expected_size=size)
    return out


def px_to_bbox(bbox_px: tuple[float, float, float, float]) -> tuple[int, int, int, int]:
    """Floor/ceil a float pixel window to integer half-open bounds [lo, hi)."""
    x0, y0, x1, y1 = bbox_px
    return int(math.floor(x0)), int(math.floor(y0)), int(math.ceil(x1)), int(math.ceil(y1))


def crop_px(df: pd.DataFrame, bbox: tuple[int, int, int, int],
            xcol: str, ycol: str) -> pd.DataFrame:
    """Half-open pixel crop [x0, x1) x [y0, y1) matching image slice bounds."""
    x0, y0, x1, y1 = bbox
    m = (df[xcol] >= x0) & (df[xcol] < x1) & (df[ycol] >= y0) & (df[ycol] < y1)
    return df[m]


def image_slice_bbox(bbox: tuple[int, int, int, int], shape_yx: tuple[int, int]) -> tuple[int, int, int, int]:
    """Clip a (x0, y0, x1, y1) crop to the image bounds -> (y0, y1, x0, x1) slice."""
    x0, y0, x1, y1 = bbox
    h, w = shape_yx
    return max(y0, 0), min(y1, h), max(x0, 0), min(x1, w)


# ---------------------------------------------------------------------------
# meta.json construction (contract)
# ---------------------------------------------------------------------------

def build_meta(
    spec: dict,
    *,
    bbox_um: Sequence[float],
    z_range_um: Sequence[float] | None,
    n_molecules: int,
    n_genes: int,
    n_cells: int,
    crop_note: str,
    images: list[dict],
    prior: str,
    source_extra: dict,
    difficulty_notes: str,
) -> dict:
    bcfg = spec["baysor"]
    area = float((bbox_um[2] - bbox_um[0]) * (bbox_um[3] - bbox_um[1]))
    dens_cells = (n_cells / area) * 1e6 if area > 0 else 0.0
    prior_conf = bcfg.get("prior_confidence") if prior != "none" else None
    return {
        "id": spec["id"],
        "kind": "real",
        "tier": spec["tier"],
        "platform": spec["platform"],
        "source": {
            "url": spec["source"]["url"],
            "doi": spec["source"].get("doi"),
            "license": spec["source"].get("license"),
            "original_dataset": spec["source"]["original_dataset"],
            "retrieved": TODAY,
            **source_extra,
        },
        "crop": {
            "bbox_um": [float(v) for v in bbox_um],
            "z_range_um": [float(v) for v in z_range_um] if z_range_um else None,
            "note": crop_note,
        },
        "stats": {
            "n_molecules": int(n_molecules),
            "n_genes": int(n_genes),
            "area_um2": round(area, 3),
            "molecules_per_um2": round(n_molecules / area, 6) if area > 0 else 0.0,
            "n_vendor_cells": int(n_cells),
            "vendor_cells_per_mm2": round(dens_cells, 3),
        },
        "difficulty": {
            "cell_density": U.cell_density_class(dens_cells),
            "gene_panel": U.gene_panel_class(n_genes),
            "notes": difficulty_notes,
        },
        "baysor": {
            "scale_um": bcfg["scale_um"],
            "scale_std": bcfg["scale_std"],
            "min_molecules_per_cell": bcfg["min_molecules_per_cell"],
            "prior": prior,
            "prior_confidence": prior_conf,
            "config": bcfg["config"],
            "extra_args": list(bcfg.get("extra_args", [])),
        },
        "images": images,
        "truth": None,
    }


def write_dataset_readme(spec: dict, meta: dict, lines: Sequence[str]) -> None:
    ds = dataset_dir(spec["id"])
    body = [
        f"# {spec['id']}",
        "",
        f"{spec['platform']} — {spec['tissue']} ({spec['tier']} tier).",
        "",
        *lines,
        "",
        "## Rebuild",
        "",
        "```bash",
        f"python benchmarks/fetch/other.py build --only {spec['id']}",
        "```",
        "",
    ]
    (ds / "README.md").write_text("\n".join(body))


def finalize_dataset(
    spec: dict,
    df: pd.DataFrame,
    *,
    meta: dict,
    readme_lines: Sequence[str],
) -> Path:
    """Validate the contract and write molecules.parquet + meta.json + README."""
    ds = dataset_dir(spec["id"])
    ds.mkdir(parents=True, exist_ok=True)
    required = {"x", "y", "gene"}
    if not required.issubset(df.columns):
        raise RuntimeError(f"{spec['id']}: missing columns {required - set(df.columns)}")
    n = len(df)
    if spec["tier"] == "quick" and n > 150_000:
        raise RuntimeError(f"{spec['id']}: quick crop has {n} molecules (>150k)")
    if spec["tier"] == "full" and n > 3_000_000:
        raise RuntimeError(f"{spec['id']}: full crop has {n} molecules (>3M)")
    if meta["stats"]["n_molecules"] != n:
        raise RuntimeError(f"{spec['id']}: meta molecule count mismatch")
    keep = [c for c in ("x", "y", "z", "gene", "prior", "cell_vendor") if c in df.columns]
    df = U.sort_molecules(df[keep])
    U.write_molecules_parquet(df, ds / "molecules.parquet")
    U.json_dump(meta, ds / "meta.json")
    write_dataset_readme(spec, meta, readme_lines)
    return ds


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def _window_for(df: pd.DataFrame, spec: dict, xcol: str, ycol: str,
                *, cap: int | None = None, target: int | None = None,
                searchable: pd.DataFrame | None = None) -> tuple[int, int, int, int]:
    """Densest-window pixel crop for a crop spec (deterministic by seed).

    ``searchable`` optionally restricts where the window may be placed (used
    where only part of the field carries reference annotations).
    """
    crop = spec["crop"]
    target = target if target is not None else crop["target_molecules"]
    cap = cap if cap is not None else crop.get("cap", 150_000)
    src = searchable if searchable is not None else df
    bbox = U.densest_window(src[xcol].to_numpy(), src[ycol].to_numpy(),
                            target=target, cap=cap, seed=spec["seed"])
    x0, y0, x1, y1 = px_to_bbox(bbox)
    # px_to_bbox floors/ceils, which can widen the window past the cap;
    # shrink deterministically (longest side first) until under it.
    for _ in range(2000):
        m = ((df[xcol] >= x0) & (df[xcol] < x1)
             & (df[ycol] >= y0) & (df[ycol] < y1))
        n = int(m.sum())
        if n <= cap:
            break
        if (x1 - x0) >= (y1 - y0):
            x1 -= 1
        else:
            y1 -= 1
    else:
        raise RuntimeError("could not shrink window under the molecule cap")
    return x0, y0, x1, y1


def build_iss(spec: dict, *, download: bool = True) -> Path:
    """ISS mouse hippocampus (Qian et al. / pciSeq), section 3-3 right."""
    cache = U.cache_dir("real_other")
    if download:
        paths = ensure_sources(
            ["iss_molecules", "iss_mask", "iss_dapi_jpg", "iss_assign_left", "iss_assign_right"],
            browser_keys=["iss_molecules", "iss_mask", "iss_dapi_jpg",
                          "iss_assign_left", "iss_assign_right"],
        )
    else:
        paths = {k: cache / Path(SOURCE_FILES[k][0]).name for k in
                 ("iss_molecules", "iss_mask", "iss_dapi_jpg", "iss_assign_left", "iss_assign_right")}

    df = pd.read_csv(paths["iss_molecules"], usecols=["x", "y", "gene"])
    all_genes = set(df["gene"].unique())
    df = df[U.gene_mask(df["gene"], spec.get("exclude_genes", []))]
    excluded = sorted(all_genes - set(df["gene"].unique()))

    # pciSeq spot -> cell assignments (CA1 only, translated into the section
    # frame); resolved against the full section so the crop may fall anywhere.
    assign = _load_iss_assignments(paths["iss_assign_left"], paths["iss_assign_right"], df)
    if assign is not None:
        df = df.merge(assign, on=["gene", "x", "y"], how="left")
        df["cell_vendor"] = df["parent_id"].map(
            lambda v: "" if pd.isna(v) else str(int(v))
        )
        df = df.drop(columns=["parent_id"])
        # Keep the window inside the CA1 band so the crop actually carries
        # reference assignments.
        searchable = df[(df["x"] >= assign["x"].min()) & (df["x"] <= assign["x"].max())
                        & (df["y"] >= assign["y"].min()) & (df["y"] <= assign["y"].max())]
        if len(searchable) < 0.1 * len(df):
            searchable = df
        assign_note = (
            "cell_vendor: pciSeq spot-to-cell assignments (figshare article 10318610) "
            "joined by gene+coordinates after per-file translation inference; the "
            "window is restricted to the CA1 band where these assignments exist."
        )
    else:
        searchable = None
        df["cell_vendor"] = ""
        assign_note = "cell_vendor: not available (assignment coordinates did not align)."

    bbox = _window_for(df, spec, "x", "y", searchable=searchable)
    df = crop_px(df, bbox, "x", "y").copy()

    # Published watershed DAPI mask -> prior image (connected components = cells).
    import tifffile

    mask = tifffile.imread(paths["iss_mask"])
    sy0, sy1, sx0, sx1 = image_slice_bbox(bbox, mask.shape)
    if sy1 > mask.shape[0] or sx1 > mask.shape[1] or sy0 < 0 or sx0 < 0:
        raise RuntimeError("ISS crop falls outside the mask frame")
    labels = U.cc_label_binary(mask[sy0:sy1, sx0:sx1])
    del mask
    n_cells = int(labels.max())
    ds = dataset_dir(spec["id"])
    (ds / "images").mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(ds / "images" / "dapi_mask.tif", labels)

    # DAPI intensity image (grayscale crop of the figshare jpg).
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None
    with Image.open(paths["iss_dapi_jpg"]) as im:
        dapi = im.convert("L").crop((sx0, sy0, sx1, sy1))
        dapi.save(ds / "images" / "dapi.tif")

    f = ISS_PIXEL_UM
    bbox_um = [sx0 * f, sy0 * f, sx1 * f, sy1 * f]
    df["x"] = df["x"] * f
    df["y"] = df["y"] * f

    images = [
        {"name": "dapi", "file": "images/dapi.tif", "pixel_size_um": f,
         "origin_um": [bbox_um[0], bbox_um[1]]},
        {"name": "dapi_mask", "file": "images/dapi_mask.tif", "pixel_size_um": f,
         "origin_um": [bbox_um[0], bbox_um[1]]},
    ]
    prior = "image:images/dapi_mask.tif"
    meta = build_meta(
        spec,
        bbox_um=bbox_um,
        z_range_um=None,
        n_molecules=len(df),
        n_genes=df["gene"].nunique(),
        n_cells=n_cells,
        crop_note=(
            "Densest ~150k-molecule window of the section (seeded grid search); "
            "image crops cover the identical pixel window."
        ),
        images=images,
        prior=prior,
        source_extra={
            "pixel_size_um": f,
            "pixel_size_source": "Qian et al. 2019 (PMC6349128) Online Methods: "
                                  "3 px = 1 um (RCP top-hat radius), 24 px = 8 um (nuclei)",
            "files": {k: paths[k].stat().st_size for k in paths},
        },
        difficulty_notes=(
            f"Cell density from the {n_cells} connected components of the published "
            "watershed DAPI mask inside the crop. "
            f"cell_vendor coverage in the crop: {(df['cell_vendor'] != '').mean():.1%}. "
            + assign_note
        ),
    )
    notes = [
        f"Source molecules: pklab example file pciSeq_3-3.csv "
        f"({SOURCE_FILES['iss_molecules'][1] // 10**6} MB, {len(df)} rows kept).",
        f"Pixel -> um factor: {f:.6f} (1/3 um per pixel, documented in the pciSeq paper).",
        "Blanks/negative probes: none present in the source gene list "
        f"(genes kept: {df['gene'].nunique()}).",
        f"Prior: published watershed DAPI mask (DAPI_3-3_mask.tif, "
        f"{SOURCE_FILES['iss_mask'][1] // 10**6} MB), connected-component labelled in the crop.",
        assign_note,
        "Baysor parameters follow examples/iss/README.md: --scale 6.5 px -> "
        f"{6.5 * f:.4f} um, config configs/iss.toml (min_molecules_per_cell=2).",
    ]
    if excluded:
        notes.append(f"Excluded gene patterns matched: {excluded}")
    return finalize_dataset(spec, df, meta=meta, readme_lines=notes)


def _load_iss_assignments(left_path: Path, right_path: Path, iss_df: pd.DataFrame
                          ) -> pd.DataFrame | None:
    """Join figshare pciSeq spot->cell assignments to the section coordinates.

    Each of the two CA1 files lives in its own translated coordinate frame;
    the translation is recovered by voting over same-gene coordinate pairs
    (seeded, deterministic) and validated by the exact join rate.
    """
    parts = []
    for path in (left_path, right_path):
        a = pd.read_csv(path)
        off = _infer_translation(iss_df, a, seed=0)
        if off is None:
            return None
        a = a.assign(x=a["spotX"] + off[0], y=a["spotY"] + off[1])
        joined = a.merge(iss_df[["gene", "x", "y"]], on=["gene", "x", "y"], how="inner")
        if len(joined) < 0.95 * len(a):
            return None
        parts.append(a[["gene", "x", "y", "parent_id"]].drop_duplicates(["gene", "x", "y"]))
    return pd.concat(parts, ignore_index=True)


def _infer_translation(iss_df: pd.DataFrame, assign_df: pd.DataFrame, seed: int
                       ) -> tuple[int, int] | None:
    import collections

    rng = np.random.default_rng(seed)
    sample = assign_df.iloc[
        np.sort(rng.choice(len(assign_df), size=min(400, len(assign_df)), replace=False))
    ]
    by_gene = {g: d[["x", "y"]].to_numpy() for g, d in iss_df.groupby("gene")}
    votes: collections.Counter = collections.Counter()
    for row in sample.itertuples():
        pts = by_gene.get(row.gene)
        if pts is None or len(pts) == 0:
            continue
        off = pts[::7] - np.array([row.spotX, row.spotY])
        for dx, dy in off:
            votes[(int(dx), int(dy))] += 1
    if not votes:
        return None
    # Validate the top candidates by exact join rate and keep the best.
    best_off, best_rate = None, 0.0
    for (dx, dy), _ in votes.most_common(5):
        cand = assign_df.assign(x=assign_df["spotX"] + dx, y=assign_df["spotY"] + dy)
        merged = cand.merge(iss_df[["gene", "x", "y"]], on=["gene", "x", "y"], how="inner")
        rate = len(merged) / max(len(cand), 1)
        if rate > best_rate:
            best_off, best_rate = (dx, dy), rate
    return best_off if best_rate >= 0.95 else None


def build_osmfish(spec: dict, *, download: bool = True) -> Path:
    """osmFISH mouse somatosensory cortex (Codeluppi et al.)."""
    cache = U.cache_dir("real_other")
    if download:
        paths = ensure_sources(["osmfish_molecules", "osmfish_segments"])
    else:
        paths = {k: cache / Path(SOURCE_FILES[k][0]).name
                 for k in ("osmfish_molecules", "osmfish_segments")}

    df = pd.read_csv(paths["osmfish_molecules"])
    df["gene"] = U.strip_hybridization_suffix(df["gene"])
    df = df[U.gene_mask(df["gene"], spec.get("exclude_genes", []))]

    # Vendor polyT segmentation polygons (vitessce mirror of polyT_seg.pkl).
    import json as _json

    with open(paths["osmfish_segments"]) as fh:
        raw = _json.load(fh)
    polygons = {int(k): v for k, v in raw.items()}
    poly_pts = [p for pts in polygons.values() for p in pts]
    px_min = min(p[0] for p in poly_pts)
    py_min = min(p[1] for p in poly_pts)
    px_max = max(p[0] for p in poly_pts)
    py_max = max(p[1] for p in poly_pts)
    # The mirror only carries annotated cells; restrict the window search to
    # the polygon-covered part of the field so the prior has good coverage.
    searchable = df[(df["x"] >= px_min) & (df["x"] <= px_max)
                    & (df["y"] >= py_min) & (df["y"] <= py_max)]
    if len(searchable) < 0.5 * len(df):
        searchable = df
    bbox = _window_for(df, spec, "x", "y", searchable=searchable)
    df = crop_px(df, bbox, "x", "y").copy()

    x0, y0, x1, y1 = bbox
    # Keep polygons intersecting the crop (with a margin).
    keep_ids = [
        lab for lab, pts in polygons.items()
        if max(p[0] for p in pts) >= x0 and min(p[0] for p in pts) <= x1
        and max(p[1] for p in pts) >= y0 and min(p[1] for p in pts) <= y1
    ]
    crop_polys = {lab: polygons[lab] for lab in sorted(keep_ids)}

    labels = U.lookup_labels(
        crop_polys,
        df["x"].to_numpy(dtype=float),
        df["y"].to_numpy(dtype=float),
    )
    df["prior"] = labels.astype(np.int32)
    df["cell_vendor"] = np.where(labels > 0, labels.astype(str), "")
    prior_coverage = float((labels > 0).mean())

    n_cells = len(crop_polys)
    f = OSMFISH_PIXEL_UM
    bbox_um = [x0 * f, y0 * f, x1 * f, y1 * f]
    df["x"] = df["x"] * f
    df["y"] = df["y"] * f

    # Vendor polygons for the crop, in um.
    ds = dataset_dir(spec["id"])
    (ds / "reference").mkdir(parents=True, exist_ok=True)
    poly_rows = [
        {"cell_id": lab, "x": p[0] * f, "y": p[1] * f}
        for lab, pts in crop_polys.items() for p in pts
    ]
    pd.DataFrame(poly_rows).to_parquet(ds / "reference" / "cell_polygons.parquet", index=False)

    meta = build_meta(
        spec,
        bbox_um=bbox_um,
        z_range_um=None,
        n_molecules=len(df),
        n_genes=df["gene"].nunique(),
        n_cells=n_cells,
        crop_note="Densest ~150k-molecule window (seeded grid search) of the osmFISH field.",
        images=[],
        prior="column",
        source_extra={
            "pixel_size_um": f,
            "pixel_size_source": "linnarssonlab.org/osmFISH data availability page: "
                                  "1 pixel = 0.065 um",
            "files": {k: paths[k].stat().st_size for k in paths},
        },
        difficulty_notes=(
            f"Cell density from the {n_cells} vendor polyT-segmentation polygons "
            f"intersecting the crop; prior coverage (molecules inside a vendor cell): "
            f"{prior_coverage:.1%}. prior column and cell_vendor come from the same "
            "polygons (molecule-in-polygon containment)."
        ),
    )
    notes = [
        f"Source molecules: pklab example file mRNA_coords_raw_counting.csv "
        f"({SOURCE_FILES['osmfish_molecules'][1] // 10**6} MB).",
        f"Pixel -> um factor: {f} (documented on the official osmFISH availability page); "
        "the example's --scale 82 px becomes "
        f"{82 * f:.3f} um, --scale-std 48 px becomes {48 * f:.3f} um.",
        f"Genes: {df['gene'].nunique()} (gene names taken from the example csv; no blanks).",
        "prior/cell_vendor: vitessce mirror of polyT_seg.pkl polygons "
        "(cells.segmentations.json, annotated cells only; window search restricted to "
        f"the polygon-covered region, prior coverage {prior_coverage:.1%}); "
        "reference/cell_polygons.parquet holds the crop's vendor polygon vertices in um.",
        "No stained image is published for this dataset (the polyT/DAPI image stacks "
        "are only referenced by the defunct Google Storage bucket), so images/ is empty.",
        "Baysor parameters follow examples/osm-FISH/README.md (configs/osm_fish.toml, "
        "min_molecules_per_cell=30).",
    ]
    return finalize_dataset(spec, df, meta=meta, readme_lines=notes)


def build_starmap(spec: dict, *, download: bool = True) -> Path:
    """STARmap visual cortex, 1020 genes, 3D (Wang et al. 2018)."""
    cache = U.cache_dir("real_other")
    if download:
        paths = ensure_sources(["starmap_molecules", "starmap_segmentation"],
                               browser_keys=["starmap_molecules", "starmap_segmentation"])
    else:
        paths = {k: cache / Path(SOURCE_FILES[k][0]).name
                 for k in ("starmap_molecules", "starmap_segmentation")}

    df = pd.read_csv(paths["starmap_molecules"])
    df = df[U.gene_mask(df["gene"], spec.get("exclude_genes", []))]

    bbox = _window_for(df, spec, "x", "y")
    df = crop_px(df, bbox, "x", "y").copy()

    import tifffile

    seg = tifffile.imread(paths["starmap_segmentation"])
    if seg.shape != (STARMAP_FRAME_PX[1], STARMAP_FRAME_PX[0]):
        raise RuntimeError(f"unexpected segmentation shape {seg.shape}")
    sy0, sy1, sx0, sx1 = image_slice_bbox(bbox, seg.shape)
    seg_crop = seg[sy0:sy1, sx0:sx1]
    ds = dataset_dir(spec["id"])
    (ds / "images").mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(ds / "images" / "segmentation.tif", seg_crop)

    # cell_vendor from the vendor labels (2D projection, as in the example).
    ix = df["x"].to_numpy(dtype=np.int64) - sx0
    iy = df["y"].to_numpy(dtype=np.int64) - sy0
    labels = seg_crop[iy, ix].astype(np.int32)
    df["cell_vendor"] = np.where(labels > 0, labels.astype(str), "")
    n_cells = int(len(np.unique(seg_crop[seg_crop > 0])))

    f = STARMAP_F_UM
    df["x"] = df["x"] * f
    df["y"] = df["y"] * f
    df["z"] = df["z"].astype(float) * STARMAP_Z_STEP_UM
    bbox_um = [sx0 * f, sy0 * f, sx1 * f, sy1 * f]
    zr = [float(df["z"].min()), float(df["z"].max())]

    images = [
        {"name": "segmentation", "file": "images/segmentation.tif", "pixel_size_um": f,
         "origin_um": [bbox_um[0], bbox_um[1]]},
    ]
    prior = "image:images/segmentation.tif"
    meta = build_meta(
        spec,
        bbox_um=bbox_um,
        z_range_um=zr,
        n_molecules=len(df),
        n_genes=df["gene"].nunique(),
        n_cells=n_cells,
        crop_note=(
            "Densest ~150k-molecule window of the visual-cortex strip (seeded grid "
            "search); the vendor segmentation is cropped to the identical pixel window."
        ),
        images=images,
        prior=prior,
        source_extra={
            "pixel_size_um": f,
            "pixel_size_source": (
                "derived from Wang et al. Fig. 2B caption 'Full field: 1.4 by 0.3 mm': "
                f"1400 um / {STARMAP_FRAME_PX[0]} px = {f:.6f} um/px (cross-check: "
                f"{STARMAP_FRAME_PX[1]} px -> {STARMAP_FRAME_PX[1] * f:.1f} um ~ 300 um)"
            ),
            "z_step_um": STARMAP_Z_STEP_UM,
            "z_step_source": (
                "derived from Wang et al. 'Eight-um-thick volumes ... were imaged': "
                f"8 um / {STARMAP_Z_PLANES} z intervals = {STARMAP_Z_STEP_UM:.4f} um"
            ),
            "files": {k: paths[k].stat().st_size for k in paths},
        },
        difficulty_notes=(
            f"Cell density from the {n_cells} vendor segmentation labels inside the "
            "crop. 3D dataset: z is kept (planes scaled to um)."
        ),
    )
    scale = 89.995 * f
    scale_std = 16.398 * f
    spec["baysor"]["scale_um"] = round(scale, 4)
    spec["baysor"]["scale_std"] = round(scale_std, 4)
    notes = [
        f"Source: pklab example files molecules.csv "
        f"({SOURCE_FILES['starmap_molecules'][1] // 10**6} MB, {len(df)} rows kept) and "
        f"segmentation.tiff ({SOURCE_FILES['starmap_segmentation'][1] // 10**6} MB).",
        f"Pixel -> um factor: {f:.6f} (from the documented full-field size "
        f"1.4 mm / {STARMAP_FRAME_PX[0]} px); z planes -> um with step "
        f"{STARMAP_Z_STEP_UM:.4f} um (documented 8 um volumes over "
        f"{STARMAP_Z_PLANES} intervals).",
        f"Genes: {df['gene'].nunique()}; blanks: none in source.",
        "prior: provided vendor segmentation cropped to the window "
        "(examples/STARmap/README.md image-prior run); cell_vendor from the same labels.",
        "Baysor scale follows the README no-prior run "
        f"(89.995 px -> {scale:.3f} um, scale-std 16.398 px -> {scale_std:.3f} um), "
        "config configs/starmap.toml (min_molecules_per_cell=70).",
    ]
    return finalize_dataset(spec, df, meta=meta, readme_lines=notes)


def _dryad_member(member: str, dest_name: str, *, force: bool = False) -> Path:
    """Extract one member of the Dryad release zip (ranged reads) into the cache."""
    from remotezip import RemoteZip

    dest = U.cache_dir("real_other") / "ileum" / dest_name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force and dest.stat().st_size > 0:
        return dest
    s3 = U.dryad_zip_url(DRYAD_FILE_URL)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with RemoteZip(s3) as z:
        with z.open(member) as src, open(tmp, "wb") as out:
            shutil.copyfileobj(src, out, length=1 << 20)
    tmp.replace(dest)
    return dest


def _fit_pixel_size_um(x_px: np.ndarray, x_um: np.ndarray,
                       y_px: np.ndarray, y_um: np.ndarray) -> tuple[float, float, float]:
    """Least-squares slope/intercepts of the release's px<->um coordinate pairs.

    The axes are centred separately before pooling (their stage offsets differ),
    then a common slope is fitted; intercepts follow per axis.
    """
    x_px = np.asarray(x_px, dtype=float)
    y_px = np.asarray(y_px, dtype=float)
    x_um = np.asarray(x_um, dtype=float)
    y_um = np.asarray(y_um, dtype=float)
    px_c = np.concatenate([x_px - x_px.mean(), y_px - y_px.mean()])
    um_c = np.concatenate([x_um - x_um.mean(), y_um - y_um.mean()])
    f = float((px_c * um_c).sum() / (px_c * px_c).sum())
    bx = float(x_um.mean() - f * x_px.mean())
    by = float(y_um.mean() - f * y_px.mean())
    return f, bx, by


def _maybe_uint16(arr: np.ndarray) -> np.ndarray:
    """Down-cast integer label arrays to uint16 when values allow it."""
    if arr.dtype == np.uint16:
        return arr
    if arr.dtype.kind in "iu" and int(arr.max()) <= np.iinfo(np.uint16).max:
        return arr.astype(np.uint16)
    return arr


def _as_zyx(arr: np.ndarray) -> np.ndarray:
    """Normalise a tifffile volume to (z, y, x).

    tifffile stacks multi-page TIFFs on axis 0 already; this only corrects a
    (y, x, z) layout (small last axis) if some producer wrote one.
    """
    if arr.ndim == 3 and arr.shape[-1] <= 64 and arr.shape[0] > 1024 \
            and arr.shape[-1] < min(arr.shape[0], arr.shape[1]):
        return np.moveaxis(arr, -1, 0)
    return arr


def _ensure_label_volume(vol: np.ndarray) -> np.ndarray:
    """Validate that ``vol`` holds integer cell labels; CC-label binary masks.

    Samples up to ~2M elements to decide between a label image (many distinct
    positive values), a binary mask (values {0,1}/{0,255}) and something
    unexpected (raises).
    """
    step = max(vol.size // 2_000_000, 1)
    sample = np.asarray(vol).ravel()[::step]
    u = np.unique(sample)
    pos = u[u > 0]
    if len(pos) == 0:
        raise RuntimeError("label volume has no positive values")
    if len(pos) >= 40:
        return vol  # label image
    if set(pos.tolist()) <= {1, 255, 65535}:
        if vol.ndim == 2:
            return U.cc_label_binary(vol)
        raise RuntimeError("3D binary mask: unexpected for the ileum release")
    raise RuntimeError(f"unexpected value set in label volume: {pos[:10]}")


def build_ileum(spec: dict, *, download: bool = True) -> Path:
    """MERFISH mouse ileum (Petukhov et al. 2022, Dryad doi:10.5061/dryad.jm63xsjb2)."""
    tier = spec["tier"]
    if download:
        mol_p = _dryad_member(DRYAD_MEMBERS["molecules"], "molecules.csv")
        labels_p = _dryad_member(DRYAD_MEMBERS["cellpose"], "cellpose_dapi.tif")
        memb_labels_p = _dryad_member(DRYAD_MEMBERS["cellpose_membrane"],
                                      "cellpose_membrane.tif")
        dapi_p = _dryad_member(DRYAD_MEMBERS["dapi"], "dapi_stack.tif")
        memb_p = _dryad_member(DRYAD_MEMBERS["membrane"], "membrane_stack.tif")
    else:
        cache = U.cache_dir("real_other") / "ileum"
        mol_p, labels_p = cache / "molecules.csv", cache / "cellpose_dapi.tif"
        memb_labels_p = cache / "cellpose_membrane.tif"
        dapi_p, memb_p = cache / "dapi_stack.tif", cache / "membrane_stack.tif"

    usecols = ["gene", "x_pixel", "y_pixel", "z_pixel", "x_um", "y_um", "z_um"]
    df = pd.read_csv(mol_p, usecols=usecols)
    df = df[U.gene_mask(df["gene"], spec.get("exclude_genes", []))]
    # The release stores z_pixel as the *plane position* (0, 13.77, ...), not
    # as a plane index; keep the full-frame plane ladder for index mapping.
    z_planes = np.sort(df["z_pixel"].unique())

    f, bx, by = _fit_pixel_size_um(
        df["x_pixel"].to_numpy(), df["x_um"].to_numpy(),
        df["y_pixel"].to_numpy(), df["y_um"].to_numpy(),
    )
    if abs(f - 0.108) > 0.002:
        raise RuntimeError(f"ileum pixel size {f} deviates from expected ~0.108 um/px")

    if tier == "quick":
        bbox = _window_for(df, spec, "x_um", "y_um")
        # re-express the window in pixel bounds for exact image alignment
        px0 = int(math.floor((bbox[0] - bx) / f))
        py0 = int(math.floor((bbox[1] - by) / f))
        px1 = int(math.ceil((bbox[2] - bx) / f))
        py1 = int(math.ceil((bbox[3] - by) / f))
    else:
        px0 = int(math.floor(df["x_pixel"].min()))
        py0 = int(math.floor(df["y_pixel"].min()))
        px1 = int(math.ceil(df["x_pixel"].max())) + 1
        py1 = int(math.ceil(df["y_pixel"].max())) + 1

    m = ((df["x_pixel"] >= px0) & (df["x_pixel"] < px1)
         & (df["y_pixel"] >= py0) & (df["y_pixel"] < py1))
    df = df[m].copy()

    # Vendor Cellpose labels for the crop (3D) -> cell_vendor + reference.
    # The release ships two label volumes: cellpose_dapi (nuclei) and
    # cellpose_membrane (whole cells); molecules get the whole-cell label.
    import tifffile

    vol = _as_zyx(tifffile.imread(memb_labels_p))
    nuc_vol = _as_zyx(tifffile.imread(labels_p))
    if nuc_vol.shape != vol.shape:
        raise RuntimeError(f"cellpose volumes differ: {vol.shape} vs {nuc_vol.shape}")
    if vol.ndim == 3:
        sy0, sy1, sx0, sx1 = image_slice_bbox((px0, py0, px1, py1), vol.shape[1:])
        label_crop = vol[:, sy0:sy1, sx0:sx1]
        nuc_crop = nuc_vol[:, sy0:sy1, sx0:sx1]
    else:
        sy0, sy1, sx0, sx1 = image_slice_bbox((px0, py0, px1, py1), vol.shape)
        label_crop = vol[sy0:sy1, sx0:sx1]
        nuc_crop = nuc_vol[sy0:sy1, sx0:sx1]
    del vol, nuc_vol
    label_crop = _ensure_label_volume(label_crop)
    nuc_crop = _ensure_label_volume(nuc_crop)
    ix = (df["x_pixel"].to_numpy(dtype=np.int64) - sx0).astype(np.int64)
    iy = (df["y_pixel"].to_numpy(dtype=np.int64) - sy0).astype(np.int64)
    iz = np.searchsorted(z_planes, df["z_pixel"].to_numpy(dtype=float))
    if label_crop.ndim == 3:
        if len(z_planes) != label_crop.shape[0]:
            raise RuntimeError(
                f"z-plane mismatch: {len(z_planes)} molecule planes vs "
                f"{label_crop.shape[0]} label planes"
            )
        valid = ((ix >= 0) & (ix < label_crop.shape[2]) & (iy >= 0)
                 & (iy < label_crop.shape[1]) & (iz >= 0) & (iz < label_crop.shape[0]))
        labs = np.zeros(len(df), dtype=np.int64)
        labs[valid] = label_crop[iz[valid], iy[valid], ix[valid]]
    else:
        valid = ((ix >= 0) & (ix < label_crop.shape[1]) & (iy >= 0)
                 & (iy < label_crop.shape[0]))
        labs = np.zeros(len(df), dtype=np.int64)
        labs[valid] = label_crop[iy[valid], ix[valid]]
    if valid.mean() < 0.99:
        raise RuntimeError(f"only {valid.mean():.1%} of molecules inside the label frame")
    df["cell_vendor"] = np.where(labs > 0, labs.astype(str), "")
    n_cells = int(len(np.unique(label_crop[label_crop > 0])))

    # Scale estimated from the vendor (Cellpose) cell areas in the crop:
    # per-plane cross-sections, aggregated over planes.
    scale = spec["baysor"]["scale_um"]
    if n_cells > 0 and label_crop.ndim == 3:
        plane_medians = []
        for z in range(label_crop.shape[0]):
            b = np.bincount(label_crop[z].ravel())
            b = b[1:]
            b = b[b > 0]
            if len(b):
                plane_medians.append(float(np.median(b)))
        if plane_medians:
            scale = float(np.sqrt(np.median(plane_medians) / np.pi) * f)
    elif n_cells > 0:
        b = np.bincount(label_crop.ravel())
        b = b[1:]
        b = b[b > 0]
        if len(b):
            scale = float(np.sqrt(np.median(b) / np.pi) * f)
    scale = round(min(max(scale, 3.0), 40.0), 3)
    spec["baysor"]["scale_um"] = scale

    bbox_um = [sx0 * f + bx, sy0 * f + by, sx1 * f + bx, sy1 * f + by]
    zr = [float(df["z_um"].min()), float(df["z_um"].max())]

    df = df.rename(columns={"x_um": "x", "y_um": "y", "z_um": "z"})
    df = df[["x", "y", "z", "gene", "cell_vendor"]]

    ds = dataset_dir(spec["id"])
    (ds / "images").mkdir(parents=True, exist_ok=True)
    (ds / "reference").mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(ds / "reference" / "cellpose_cell_labels.tif",
                     _maybe_uint16(label_crop), photometric="minisblack")
    tifffile.imwrite(ds / "reference" / "cellpose_nuclei_labels.tif",
                     _maybe_uint16(nuc_crop), photometric="minisblack")
    del nuc_crop

    image_files = {}
    for name, src in (("dapi", dapi_p), ("membrane", memb_p)):
        stack = _as_zyx(tifffile.imread(src))
        if stack.ndim == 3:
            crop = stack[:, sy0:sy1, sx0:sx1]
        else:
            crop = stack[sy0:sy1, sx0:sx1]
        tifffile.imwrite(ds / "images" / f"{name}.tif", crop,
                         photometric="minisblack")
        del stack
        image_files[name] = src.stat().st_size

    images = [
        {"name": "dapi", "file": "images/dapi.tif", "pixel_size_um": f,
         "origin_um": [bbox_um[0], bbox_um[1]], "z_range_um": [2.5, 14.5]},
        {"name": "membrane", "file": "images/membrane.tif", "pixel_size_um": f,
         "origin_um": [bbox_um[0], bbox_um[1]], "z_range_um": [2.5, 14.5]},
    ]
    meta = build_meta(
        spec,
        bbox_um=bbox_um,
        z_range_um=zr,
        n_molecules=len(df),
        n_genes=df["gene"].nunique(),
        n_cells=n_cells,
        crop_note=(
            "Densest ~150k-molecule window (quick) of the released ileum region."
            if tier == "quick" else
            "Full released region (all molecules after blank removal)."
        ),
        images=images,
        prior="none",
        source_extra={
            "pixel_size_um": f,
            "pixel_size_source": (
                "least-squares fit of the release's own x_pixel/y_pixel vs x_um/y_um "
                f"columns: {f:.6f} um/px (x offset {bx:.3f} um, y offset {by:.3f} um)"
            ),
            "z_range_planes": [0, 8],
            "files": {
                "molecules.csv": mol_p.stat().st_size,
                "cellpose_dapi.tif": labels_p.stat().st_size,
                "cellpose_membrane.tif": memb_labels_p.stat().st_size,
                "dapi_stack.tif": dapi_p.stat().st_size,
                "membrane_stack.tif": memb_p.stat().st_size,
            },
        },
        difficulty_notes=(
            f"Cell density from the {n_cells} published Cellpose membrane labels "
            "(whole cells) inside the crop; cell_vendor is the Cellpose cell label at "
            "each molecule (3D lookup, z mapped via the release's z_pixel plane "
            "ladder); nuclei labels are kept in reference/. Prior: none "
            "(the release's primary Baysor run was mRNA-only; its membrane-prior run "
            "used a Cellpose prior). qc_score filtering (>0.8) was applied by the data "
            "release."
        ),
    )
    notes = [
        f"Source: Dryad data_release_baysor_merfish_gut.zip "
        f"(735,805,624 B total; molecules.csv {mol_p.stat().st_size} B), members fetched "
        "with ranged reads of the presigned S3 zip.",
        f"Coordinates: release x_um/y_um/z_um columns used directly; pixel factor "
        f"{f:.6f} um/px (fitted from the release's pixel/um column pairs) applies to the "
        "image and label crops only. Images: 9 z-planes, 1.5 um spacing, 2.5-14.5 um "
        "above the coverslip (release README).",
        f"Genes: {df['gene'].nunique()} after removing Blank* codewords "
        "(release-provided filter: qc_score > 0.8).",
        "cell_vendor: published Cellpose segmentation (release data_analysis/cellpose; "
        "membrane channel = whole cells), molecule-to-label containment; label crops "
        "kept in reference/cellpose_cell_labels.tif and reference/cellpose_nuclei_labels.tif.",
        "Images: cropped DAPI and Na+/K+-ATPase membrane stacks (images/dapi.tif, "
        "images/membrane.tif).",
        "Baysor scale estimated as the median equivalent radius of the vendor Cellpose "
        f"cells ({scale} um); min_molecules_per_cell=30 follows the release's published "
        "CLI parameters; config configs/example_config.toml.",
    ]
    return finalize_dataset(spec, df, meta=meta, readme_lines=notes)


def _choose_fov(metadata: pd.DataFrame, fov_col: str = "fov") -> int:
    """FOV with the most vendor cells (ties -> lowest FOV index)."""
    counts = metadata.groupby(fov_col).size()
    counts = counts.sort_index()
    return int(counts.idxmax())


def _stream_tx_fov(tx_path: Path, *, gzipped: bool, fov: int,
                   exclude_patterns: Sequence[str],
                   cache_p: Path | None = None) -> tuple[pd.DataFrame, dict]:
    """Read a cached CosMx transcript file and keep one FOV; rows + stats.

    The full tx file is downloaded once into the cache (with size check), then
    read in chunks; the kept rows (one FOV, genes filtered) are cached as
    parquet so the quick/full builds share a single pass.
    """
    stats_p = None
    if cache_p is not None:
        stats_p = cache_p.with_suffix(".stats.json")
        if cache_p.exists() and stats_p.exists():
            return pd.read_parquet(cache_p), json.loads(stats_p.read_text())
    usecols = (["fov", "cell_ID", "x_global_px", "y_global_px", "target"]
               if not gzipped else
               ["fov", "cell_ID", "cell", "x_global_px", "y_global_px", "target"])
    parts = []
    total = 0
    excluded_targets: set[str] = set()
    reader = pd.read_csv(tx_path, usecols=usecols, chunksize=500_000,
                         compression="gzip" if gzipped else None)
    for chunk in reader:
        total += len(chunk)
        keep = U.gene_mask(chunk["target"], exclude_patterns)
        excluded_targets.update(chunk.loc[~keep, "target"].astype(str).unique().tolist())
        part = chunk[(chunk["fov"] == fov) & keep]
        if len(part):
            parts.append(part)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=usecols)
    stats = {"total_transcripts": total, "excluded_targets": sorted(excluded_targets)}
    if cache_p is not None:
        cache_p.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache_p, index=False)
        stats_p.write_text(json.dumps(stats, indent=2))
    return df, stats


def _tx_cache_path(cache: Path, fov: int, patterns: Sequence[str]) -> Path:
    import hashlib as _hashlib

    key = _hashlib.sha1(",".join(patterns).encode()).hexdigest()[:8]
    return cache / f"fov{fov}_tx_{key}.parquet"


def _summarize_excluded(items: Sequence[str]) -> dict:
    """Compact record of removed control targets (full lists can be huge)."""
    return {"count": len(items), "examples": list(items[:12])}


def _excluded_note(items: Sequence[str]) -> str:
    if not items:
        return "none"
    head = ", ".join(items[:6])
    return f"{len(items)} control targets removed (e.g. {head}, ...)"


def build_cosmx_nsclc(spec: dict, *, download: bool = True) -> Path:
    """CosMx NSCLC Lung5_Rep1 (960-plex FFPE, NanoString/Bruker public data)."""
    cache = U.cache_dir("real_other") / "cosmx_nsclc"
    cache.mkdir(parents=True, exist_ok=True)
    if download:
        paths = ensure_sources(["nsclc_metadata", "nsclc_fov"])
    else:
        paths = {k: U.cache_dir("real_other") / Path(SOURCE_FILES[k][0]).name
                 for k in ("nsclc_metadata", "nsclc_fov")}

    metadata = pd.read_csv(paths["nsclc_metadata"])
    fov = _choose_fov(metadata)
    tx_path = ensure_sources(["nsclc_tx"])["nsclc_tx"]
    tx, stats = _stream_tx_fov(
        tx_path, gzipped=False, fov=fov,
        exclude_patterns=spec.get("exclude_genes", []),
        cache_p=_tx_cache_path(cache, fov, spec.get("exclude_genes", [])),
    )
    if tx.empty:
        raise RuntimeError("no transcripts for chosen FOV")

    tier = spec["tier"]
    if tier == "quick":
        bbox = _window_for(tx, spec, "x_global_px", "y_global_px")
        px0, py0, px1, py1 = px_to_bbox(bbox)
    else:
        px0 = int(math.floor(tx["x_global_px"].min()))
        py0 = int(math.floor(tx["y_global_px"].min()))
        px1 = int(math.ceil(tx["x_global_px"].max())) + 1
        py1 = int(math.ceil(tx["y_global_px"].max())) + 1
    tx = crop_px(tx, (px0, py0, px1, py1), "x_global_px", "y_global_px").copy()

    f = NSCLC_PIXEL_UM
    bbox_um = [px0 * f, py0 * f, px1 * f, py1 * f]
    tx["x"] = tx["x_global_px"] * f
    tx["y"] = tx["y_global_px"] * f
    tx["cell_vendor"] = np.where(
        tx["cell_ID"] > 0, f"{fov}_" + tx["cell_ID"].astype(str), ""
    )
    df = tx[["x", "y", "target", "cell_vendor"]].rename(columns={"target": "gene"})

    # Vendor cells inside the crop (from the cell metadata file).
    fov_meta = metadata[metadata["fov"] == fov]
    gx = fov_meta["CenterX_global_px"].to_numpy(dtype=float) * f
    gy = fov_meta["CenterY_global_px"].to_numpy(dtype=float) * f
    in_crop = ((gx >= bbox_um[0]) & (gx < bbox_um[2])
               & (gy >= bbox_um[1]) & (gy < bbox_um[3]))
    cells = fov_meta[in_crop].copy()
    cells["x_um"] = cells["CenterX_global_px"].astype(float) * f
    cells["y_um"] = cells["CenterY_global_px"].astype(float) * f
    cells["area_um2"] = cells["Area"].astype(float) * f * f
    n_cells = len(cells)
    scale = float(np.sqrt(np.median(cells["area_um2"]) / np.pi)) if n_cells else 10.0
    scale = round(min(max(scale, 4.0), 30.0), 3)
    spec["baysor"]["scale_um"] = scale

    # Vendor cell labels for the crop (per-FOV CellLabels tiff).
    celllabels_url = f"{NSCLC_BUCKET}/SMI/Lung5_Rep1/CellLabels/CellLabels_F{fov:03d}.tif"
    label_path = U.ensure_file(celllabels_url, cache / f"CellLabels_F{fov:03d}.tif")
    import tifffile

    labels_img = tifffile.imread(label_path)
    fov_positions = pd.read_csv(paths["nsclc_fov"])
    origin = fov_positions[fov_positions["fov"] == fov].iloc[0]
    ox, oy = float(origin["x_px"]), float(origin["y_px"])
    lx0, lx1 = int(math.floor(px0 - ox)), int(math.ceil(px1 - ox))
    ly0, ly1 = int(math.floor(py0 - oy)), int(math.ceil(py1 - oy))
    sy0, sy1, sx0, sx1 = image_slice_bbox((lx0, ly0, lx1, ly1), labels_img.shape)
    label_crop = labels_img[sy0:sy1, sx0:sx1]
    if label_crop.size == 0:
        raise RuntimeError("CellLabels crop is empty")

    ds = dataset_dir(spec["id"])
    (ds / "reference").mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(ds / "reference" / "cell_labels.tif", label_crop)

    meta = build_meta(
        spec,
        bbox_um=bbox_um,
        z_range_um=None,
        n_molecules=len(df),
        n_genes=df["gene"].nunique(),
        n_cells=n_cells,
        crop_note=(
            f"FOV {fov} (densest window, ~150k molecules) of CosMx NSCLC Lung5_Rep1."
            if tier == "quick" else
            f"Full FOV {fov} of CosMx NSCLC Lung5_Rep1."
        ),
        images=[],
        prior="none",
        source_extra={
            "pixel_size_um": f,
            "pixel_size_source": "SMI Data File ReadMe (S3): x_local_px * 0.18 um/px",
            "fov": fov,
            "files": {k: paths[k].stat().st_size for k in paths},
            "source_total_transcripts": stats["total_transcripts"],
            "excluded_targets": _summarize_excluded(stats["excluded_targets"]),
        },
        difficulty_notes=(
            f"Cell density from {n_cells} vendor cells (metadata) with centroids inside "
            f"the crop; cell_vendor = fov_cellID of the vendor transcript assignment. "
            f"Vendor cell labels kept in reference/cell_labels.tif. Excluded targets: "
            f"{_excluded_note(stats['excluded_targets'])}."
        ),
    )
    notes = [
        f"Source: nanostring-public-share S3, Lung5_Rep1_tx_file.csv "
        f"({SOURCE_FILES['nsclc_tx'][1]:,} B downloaded to cache; metadata "
        f"{SOURCE_FILES['nsclc_metadata'][1]} B).",
        f"Pixel -> um factor: {f} (documented in the SMI ReadMe). z-slices exist in the "
        "source but their physical spacing is not documented, so z is dropped (2D crop).",
        f"Genes: {df['gene'].nunique()} targets kept; excluded: "
        f"{_excluded_note(stats['excluded_targets'])}.",
        f"Chosen FOV: {fov} (most vendor cells); source has {len(metadata.groupby('fov'))} FOVs, "
        f"{stats['total_transcripts']:,} transcripts total.",
        "cell_vendor from the transcript file's cell_ID (0 = unassigned -> empty); "
        "vendor cell labels in reference/cell_labels.tif (per-FOV CellLabels tiff).",
        "No single-channel stain images are shipped in the public flat data (raw "
        "morphology images are a 41.7 GB tar; CellComposite jpgs are RGB composites), "
        "so images/ is empty.",
        f"Baysor scale estimated from vendor cell areas (median equivalent radius = "
        f"{scale} um); min_molecules_per_cell=20; config configs/example_config.toml.",
    ]
    return finalize_dataset(spec, df, meta=meta, readme_lines=notes)


def build_cosmx_wtx(spec: dict, *, download: bool = True) -> Path:
    """CosMx Whole Transcriptome (WTX) colon discovery dataset (~19k targets)."""
    cache = U.cache_dir("real_other") / "cosmx_wtx"
    cache.mkdir(parents=True, exist_ok=True)
    if download:
        paths = ensure_sources(["wtx_metadata", "wtx_fov", "wtx_readme"])
    else:
        paths = {k: U.cache_dir("real_other") / Path(SOURCE_FILES[k][0]).name
                 for k in ("wtx_metadata", "wtx_fov", "wtx_readme")}

    metadata = pd.read_csv(paths["wtx_metadata"], low_memory=False)
    fov = _choose_fov(metadata, "fov")
    tx_path = ensure_sources(["wtx_tx"])["wtx_tx"]
    tx, stats = _stream_tx_fov(
        tx_path, gzipped=True, fov=fov,
        exclude_patterns=spec.get("exclude_genes", []),
        cache_p=_tx_cache_path(cache, fov, spec.get("exclude_genes", [])),
    )
    if tx.empty:
        raise RuntimeError("no transcripts for chosen FOV")

    tier = spec["tier"]
    if tier == "quick":
        bbox = _window_for(tx, spec, "x_global_px", "y_global_px")
        px0, py0, px1, py1 = px_to_bbox(bbox)
    else:
        px0 = int(math.floor(tx["x_global_px"].min()))
        py0 = int(math.floor(tx["y_global_px"].min()))
        px1 = int(math.ceil(tx["x_global_px"].max())) + 1
        py1 = int(math.ceil(tx["y_global_px"].max())) + 1
    tx = crop_px(tx, (px0, py0, px1, py1), "x_global_px", "y_global_px").copy()

    f = WTX_PIXEL_UM
    bbox_um = [px0 * f, py0 * f, px1 * f, py1 * f]
    tx["x"] = tx["x_global_px"] * f
    tx["y"] = tx["y_global_px"] * f
    tx["cell_vendor"] = np.where(tx["cell_ID"] > 0, tx["cell"].astype(str), "")
    df = tx[["x", "y", "target", "cell_vendor"]].rename(columns={"target": "gene"})

    # Vendor cells inside the crop.
    fov_meta = metadata[metadata["fov"] == fov]
    ccols = [c for c in ("CenterX_global_px", "CenterY_global_px") if c in fov_meta.columns]
    if len(ccols) == 2:
        gx = fov_meta[ccols[0]].astype(float) * f
        gy = fov_meta[ccols[1]].astype(float) * f
    else:  # fall back to FOV-local centres + FOV origin
        origin = pd.read_csv(paths["wtx_fov"])
        row = origin[origin["FOV"] == fov].iloc[0]
        gx = fov_meta["CenterX_local_px"].astype(float) * f + float(row["x_global_px"]) * f
        gy = fov_meta["CenterY_local_px"].astype(float) * f + float(row["y_global_px"]) * f
    in_crop = ((gx >= bbox_um[0]) & (gx < bbox_um[2])
               & (gy >= bbox_um[1]) & (gy < bbox_um[3]))
    cells = fov_meta[in_crop].copy()
    cells["x_um"] = gx[in_crop].to_numpy()
    cells["y_um"] = gy[in_crop].to_numpy()
    if "Area" in cells.columns:
        cells["area_um2"] = cells["Area"].astype(float) * f * f
    n_cells = len(cells)
    if n_cells and "area_um2" in cells.columns:
        scale = float(np.sqrt(np.median(cells["area_um2"]) / np.pi))
    else:
        scale = 8.0
    scale = round(min(max(scale, 4.0), 30.0), 3)
    spec["baysor"]["scale_um"] = scale

    ds = dataset_dir(spec["id"])
    (ds / "reference").mkdir(parents=True, exist_ok=True)
    keep_cols = [c for c in cells.columns if c in (
        "fov", "cell_id", "cell_ID", "Area", "NucArea", "x_um", "y_um", "area_um2")]
    cells[keep_cols].to_parquet(ds / "reference" / "vendor_cells.parquet", index=False)

    meta = build_meta(
        spec,
        bbox_um=bbox_um,
        z_range_um=None,
        n_molecules=len(df),
        n_genes=df["gene"].nunique(),
        n_cells=n_cells,
        crop_note=(
            f"FOV {fov} (densest window, ~150k molecules) of the WTX colon dataset."
            if tier == "quick" else
            f"Full FOV {fov} of the WTX colon dataset."
        ),
        images=[],
        prior="none",
        source_extra={
            "pixel_size_um": f,
            "pixel_size_source": "README_coloncancer.html: pixel edge 120 nm, "
                                  "multiply px by 0.12028 um",
            "fov": fov,
            "files": {k: paths[k].stat().st_size for k in paths if paths[k].exists()},
            "source_total_transcripts": stats["total_transcripts"],
            "excluded_targets": _summarize_excluded(stats["excluded_targets"]),
        },
        difficulty_notes=(
            f"Cell density from {n_cells} vendor cells (metadata) with centroids inside "
            "the crop; cell_vendor is the study-wide cell id of the vendor transcript "
            "assignment (empty when cell_ID = 0). Excluded targets: "
            f"{_excluded_note(stats['excluded_targets'])}. " + WTX_PERF_NOTE
        ),
    )
    notes = [
        f"Source: objects.liquidweb.services smi-public bucket, colon_discovery "
        f"S0_tx_file.csv.gz ({SOURCE_FILES['wtx_tx'][1]:,} B downloaded to cache); "
        f"metadata ({SOURCE_FILES['wtx_metadata'][1]} B) and fov positions cached.",
        f"Pixel -> um factor: {f} (documented in README_coloncancer.html). z planes exist "
        "but the physical step is undocumented, so z is dropped (2D crop).",
        f"Genes: {df['gene'].nunique()} real targets kept; excluded control targets: "
        f"{_excluded_note(stats['excluded_targets'])}. The panel is the CosMx Human Whole "
        "Transcriptome panel (README: 'Genes Above LOD 7579' of the full panel).",
        f"Chosen FOV: {fov} (most vendor cells); source has "
        f"{len(metadata.groupby('fov'))} FOVs, {stats['total_transcripts']:,} transcripts.",
        "cell_vendor from the tx file (cell column; cell_ID=0 -> empty); vendor cell "
        "centroids/areas for the crop in reference/vendor_cells.parquet.",
        "No single-channel stains in the public files (Napari.zip is 68.8 GB), so "
        "images/ is empty.",
        f"Baysor scale estimated from vendor cell areas (median equivalent radius = "
        f"{scale} um); min_molecules_per_cell=20; config configs/example_config.toml; "
        "extra_args: --cluster-method louvain.",
        WTX_PERF_NOTE,
    ]
    return finalize_dataset(spec, df, meta=meta, readme_lines=notes)


BUILDERS: dict[str, Callable[..., Path]] = {
    "iss": build_iss,
    "osmfish": build_osmfish,
    "starmap": build_starmap,
    "ileum": build_ileum,
    "cosmx_nsclc": build_cosmx_nsclc,
    "cosmx_wtx": build_cosmx_wtx,
}


# ---------------------------------------------------------------------------
# Smoke runs (Release binary on every quick crop)
# ---------------------------------------------------------------------------

def parse_time_v(path: Path) -> tuple[float, int]:
    """Extract (wall seconds, max RSS KB) from ``/usr/bin/time -v`` output."""
    text = path.read_text(errors="replace")
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    rss = int(m.group(1)) if m else -1
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", text)
    wall = -1.0
    if m:
        parts = [float(p) for p in m.group(1).split(":")]
        while len(parts) < 3:
            parts.insert(0, 0.0)
        wall = parts[0] * 3600 + parts[1] * 60 + parts[2]
    return wall, rss


def build_baysor_command(spec: dict, ds: Path, binary: Path, out_dir: Path,
                         time_file: Path,
                         extra_args: Sequence[str] | None = None) -> list[str]:
    meta = json.loads((ds / "meta.json").read_text())
    b = meta["baysor"]
    cfg = U.repo_root() / b["config"]
    time_file.parent.mkdir(parents=True, exist_ok=True)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    args = b.get("extra_args", []) if extra_args is None else list(extra_args)
    cmd = [
        "/usr/bin/time", "-v", "-o", str(time_file),
        "taskset", "-c", "0-5", str(binary), "run",
        "-c", str(cfg),
        "-s", str(b["scale_um"]),
        "--scale-std", str(b["scale_std"]),
        "-m", str(b["min_molecules_per_cell"]),
    ]
    if b.get("prior_confidence") is not None:
        cmd += ["--prior-segmentation-confidence", str(b["prior_confidence"])]
    cmd += list(args)
    cmd += ["-o", str(out_dir), str(ds / "molecules.parquet")]
    prior = b["prior"]
    if prior == "column":
        cmd.append(":prior")
    elif prior.startswith("image:"):
        cmd.append(str(ds / prior[len("image:"):]))
    return cmd


def _run_baysor_capped(
    cmd: list[str], *, env: dict, timeout_s: float
) -> tuple[int, str, str, float, bool]:
    """Run Baysor with a wall-clock cap; kill the whole process group on timeout.

    Returns ``(returncode, stdout, stderr, wall_seconds, timed_out)``.
    """
    import signal

    t0 = time.time()
    proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True,
                            start_new_session=True)
    timed_out = False
    try:
        stdout, stderr = proc.communicate(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        timed_out = True
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            stdout, stderr = proc.communicate(timeout=15)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            stdout, stderr = proc.communicate()
    return proc.returncode, stdout or "", stderr or "", time.time() - t0, timed_out


def _last_log_stage(stdout: str, stderr: str) -> str:
    """Last log line of a (possibly killed) Baysor run: the stage it was at."""
    lines = [ln.strip() for ln in (stdout + "\n" + stderr).splitlines() if ln.strip()]
    info = [ln for ln in lines if "[info]" in ln]
    return (info or lines)[-1] if (info or lines) else ""


def smoke(manifest: dict, binary: Path, *, timeout_s: float = 1800.0,
          force: bool = False) -> Path:
    """Run the Release binary on every quick crop and record wall time + peak RSS.

    Each run is capped at ``timeout_s`` (default 30 min) and killed on timeout;
    timeouts are recorded as such instead of failing the pipeline.  Datasets
    with ``smoke_record_default_timeout`` get an extra run with the DEFAULT
    cluster method (extra_args stripped), capped at that many seconds, to
    document pathological default settings (see the WTX datasets).
    """
    runs_root = U.bench_data_root() / "runs" / "real_other_smoke"
    runs_root.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {
        "binary": str(binary),
        "generated": TODAY,
        "timeout_seconds": timeout_s,
        "runs": {},
        "recorded_default_timeouts": {},
    }
    env = {**os.environ, "OMP_NUM_THREADS": "6"}
    for spec in manifest["datasets"]:
        if spec["tier"] != "quick":
            continue
        ds = dataset_dir(spec["id"])
        meta = json.loads((ds / "meta.json").read_text())
        out_dir = runs_root / spec["id"]
        time_file = runs_root / f"{spec['id']}.time"
        if out_dir.exists():
            shutil.rmtree(out_dir)
        cmd = build_baysor_command(spec, ds, binary, out_dir, time_file)
        rc, stdout, stderr, wall_raw, timed_out = _run_baysor_capped(
            cmd, env=env, timeout_s=timeout_s
        )
        wall, rss = parse_time_v(time_file) if time_file.exists() else (-1.0, -1)
        if wall < 0:
            wall = wall_raw
        entry = {
            "command": cmd,
            "exit_code": rc,
            "status": "timeout" if timed_out else "ok",
            "wall_seconds": round(wall, 2),
            "max_rss_kb": rss,
            "n_molecules": meta["stats"]["n_molecules"],
            "n_genes": meta["stats"]["n_genes"],
            "prior": meta["baysor"]["prior"],
        }
        if timed_out:
            entry["stage"] = _last_log_stage(stdout, stderr)
            print(f"{spec['id']}: TIMEOUT after {wall:.0f}s at: {entry['stage']}")
        elif rc != 0:
            print(stdout[-3000:], file=sys.stderr)
            print(stderr[-3000:], file=sys.stderr)
            raise RuntimeError(f"baysor failed on {spec['id']} (exit {rc})")
        else:
            print(f"{spec['id']}: wall {wall:.1f}s, peak RSS {rss // 1024} MiB")
        results["runs"][spec["id"]] = entry

        # Optional: document how the DEFAULT cluster method behaves (capped).
        cap = spec.get("smoke_record_default_timeout")
        if cap:
            out_def = runs_root / f"{spec['id']}_default_mrf"
            tf_def = runs_root / f"{spec['id']}_default_mrf.time"
            if out_def.exists():
                shutil.rmtree(out_def)
            cmd_def = build_baysor_command(spec, ds, binary, out_def, tf_def,
                                           extra_args=[])
            rc2, so2, se2, wall2, to2 = _run_baysor_capped(
                cmd_def, env=env, timeout_s=float(cap)
            )
            wall_d, rss_d = parse_time_v(tf_def) if tf_def.exists() else (-1.0, -1)
            if wall_d < 0:
                wall_d = wall2
            rec = {
                "command": cmd_def,
                "cap_seconds": float(cap),
                "status": "timeout" if to2 else ("ok" if rc2 == 0 else "error"),
                "exit_code": rc2,
                "wall_seconds": round(wall_d, 2),
                "max_rss_kb": rss_d,
                "stage": _last_log_stage(so2, se2),
                "note": (
                    "default config (cluster-method mrf, ICA init) on a ~19k-gene crop; "
                    "recorded as a performance-stress datapoint (perf-stress: default "
                    "ICA init)"
                ),
            }
            results["recorded_default_timeouts"][spec["id"]] = rec
            print(f"{spec['id']} (default mrf): {rec['status']} after "
                  f"{wall_d:.0f}s at: {rec['stage']}")
    out = U.repo_root() / "benchmarks" / "baselines" / "real_other_smoke.json"
    U.json_dump(results, out)
    return out


def report(manifest: dict) -> str:
    rows = []
    for spec in manifest["datasets"]:
        ds = dataset_dir(spec["id"])
        meta = json.loads((ds / "meta.json").read_text())
        st = meta["stats"]
        size = sum(p.stat().st_size for p in ds.rglob("*") if p.is_file())
        rows.append(
            f"| {spec['id']} | {spec['tier']} | {meta['platform']} | {st['n_genes']} | "
            f"{st['n_molecules']:,} | {st['n_vendor_cells']:,} | "
            f"{meta['difficulty']['gene_panel']} | {meta['difficulty']['cell_density']} | "
            f"{size // 10**6} MB |"
        )
    header = ("| id | tier | platform | genes | molecules | vendor cells | panel | density | disk |\n"
              "|---|---|---|---|---|---|---|---|---|")
    out = "\n".join([header, *rows])
    print(out)
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build", help="build datasets from the manifest")
    b.add_argument("--manifest", type=Path, default=None)
    b.add_argument("--only", action="append", default=None,
                   help="dataset id (repeatable); default: all")
    b.add_argument("--no-download", action="store_true",
                   help="use only files already in the cache")
    b.add_argument("--force", action="store_true", help="rebuild even if outputs exist")
    s = sub.add_parser("smoke", help="run Release Baysor on all quick crops")
    s.add_argument("--manifest", type=Path, default=None)
    s.add_argument("--timeout", type=float, default=1800.0,
                   help="wall-clock cap per run in seconds (default 1800 = 30 min)")
    s.add_argument("--binary", type=Path,
                   default=Path("/home/vpetukhov/Projects/Baysor/.bench-data/binaries/"
                                "baysor-bugfixes-35e8a7e"))
    s.add_argument("--force", action="store_true")
    r = sub.add_parser("report", help="print per-dataset stats as markdown")
    r.add_argument("--manifest", type=Path, default=None)

    args = ap.parse_args(argv)
    manifest = load_manifest(args.manifest)

    if args.cmd == "build":
        specs = [d for d in manifest["datasets"]
                 if args.only is None or d["id"] in set(args.only)]
        if not specs:
            raise SystemExit("no matching datasets")
        for spec in specs:
            ds = dataset_dir(spec["id"])
            if args.force and ds.exists():
                shutil.rmtree(ds)
            builder = BUILDERS[spec["builder"]]
            print(f"== building {spec['id']} ({spec['builder']}, {spec['tier']})")
            t0 = time.time()
            out = builder(spec, download=not args.no_download)
            print(f"   -> {out} in {time.time() - t0:.1f}s")
        return 0
    if args.cmd == "smoke":
        out = smoke(manifest, args.binary.resolve(), timeout_s=args.timeout,
                    force=args.force)
        print(f"wrote {out}")
        return 0
    if args.cmd == "report":
        report(manifest)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
