"""Helpers for fetching and cropping public 10x Xenium datasets (BENCH-REALX).

The module contains the non-CLI parts of the pipeline:

* partial downloads of ``*_outs.zip`` bundles with :mod:`remotezip`
  (only the members listed in the manifest are ever fetched),
* streaming, column-projected reads of ``transcripts.parquet`` with gene and
  quality filtering while iterating row groups,
* deterministic crop-box selection on a molecule/cell density histogram,
* vendor ``cell_vendor``/``prior`` column construction per the dataset
  contract in ``benchmarks/README.md``,
* cropping of vendor boundary parquet files and OME-TIFF focus images,
* ``meta.json`` construction and the inventory report.

Data never lives in the repository; everything is written under
``$BAYSOR_BENCH_DATA`` (default ``<repo>/.bench-data``).
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Sequence

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

import download  # noqa: E402  (shared helper, same directory)

# ---------------------------------------------------------------------------
# Paths and manifest
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = REPO_ROOT / "benchmarks" / "datasets" / "real_xenium.yaml"

#: feature-name prefixes that mark control probes / blank codewords and are
#: dropped when the bundle has no ``is_gene``/``codeword_category`` column.
CONTROL_PREFIXES = (
    "negcontrol",
    "unassignedcodeword",
    "deprecatedcodeword",
    "intergenic_region",
    "blank",
)


def data_root() -> Path:
    """Root of the benchmark data tree (``$BAYSOR_BENCH_DATA``)."""
    env = os.environ.get("BAYSOR_BENCH_DATA")
    if env:
        return Path(env)
    return REPO_ROOT / ".bench-data"


def cache_root() -> Path:
    """Cache directory for raw Xenium bundle members of this task."""
    return data_root() / "cache" / "realx"


def load_manifest(path: Path | str | None = None) -> dict:
    import yaml

    path = Path(path) if path else DEFAULT_MANIFEST
    with open(path) as fh:
        return yaml.safe_load(fh)


def save_manifest(manifest: dict, path: Path | str | None = None) -> None:
    import yaml

    path = Path(path) if path else DEFAULT_MANIFEST
    with open(path, "w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, default_flow_style=False, allow_unicode=True)


def unique_sources(manifest: dict) -> list[dict]:
    """Unique source bundles referenced by the manifest, in first-seen order."""
    seen: dict[str, dict] = {}
    for ds in manifest["datasets"]:
        src = ds["source"]
        url = src["url"]
        if url not in seen:
            entry = dict(src)
            entry["prefix"] = ds.get("prefix", "")
            entry["members"] = dict(ds["members"])
            if ds.get("images"):
                entry["members"] = dict(entry["members"])
                entry["image_members"] = list(ds["images"]["focus_members"])
            seen[url] = entry
        elif ds.get("images") and "image_members" not in seen[url]:
            seen[url]["image_members"] = list(ds["images"]["focus_members"])
    return list(seen.values())


# ---------------------------------------------------------------------------
# Remote zip access
# ---------------------------------------------------------------------------


def head_info(url: str) -> tuple[int, int | None]:
    """HEAD ``url``; return ``(status, content_length)``."""
    import requests

    resp = requests.head(url, allow_redirects=True, timeout=60)
    length = resp.headers.get("Content-Length")
    return resp.status_code, (int(length) if length else None)


def member_cache_path(url: str, member: str, cache: Path | None = None) -> Path:
    zip_name = Path(url).name[:-4]  # strip .zip
    return (cache or cache_root()) / zip_name / member


def fetch_members(
    url: str,
    members: Sequence[str],
    cache: Path | None = None,
    *,
    expected: dict | None = None,
    record: dict | None = None,
) -> dict[str, Path]:
    """Ensure every member of the remote zip exists in the local cache.

    Streaming goes through the shared :mod:`download` helper: members are
    fetched with range requests only (the zip itself is never downloaded
    whole), retried with exponential backoff (429/5xx honouring
    Retry-After), written to a ``.part`` file and renamed atomically after
    size + sha256 checks.  ``expected`` maps member -> ``{"bytes",
    "sha256"}`` (typically from the manifest) and is verified on reuse;
    ``record`` is filled with the same structure for the caller to persist.
    Returns member -> local path.
    """
    cache = cache or cache_root()
    wanted = list(dict.fromkeys(members))  # de-duplicate, keep order
    dsts = {m: member_cache_path(url, m, cache) for m in wanted}
    exp = {m: expected[m] for m in wanted if expected and m in expected}
    digests = download.fetch_zip_members(url, dsts, expected=exp)
    if record is not None:
        for m, digest in digests.items():
            record[m] = {"bytes": dsts[m].stat().st_size, "sha256": digest}
    return dsts


def cached_bytes(cache: Path | None = None) -> int:
    """Total size of the members cached for this task."""
    cache = cache or cache_root()
    if not cache.exists():
        return 0
    return sum(p.stat().st_size for p in cache.rglob("*") if p.is_file())


# ---------------------------------------------------------------------------
# Transcript filtering
# ---------------------------------------------------------------------------


def transcript_columns(schema_names: Sequence[str], with_z: bool = False) -> list[str]:
    """Columns to project from a Xenium transcripts table.

    ``with_z`` additionally projects ``z_location`` (3D crops).
    """
    cols = ["x_location", "y_location", "qv", "cell_id", "overlaps_nucleus", "feature_name"]
    if with_z:
        cols.append("z_location")
    if "is_gene" in schema_names:
        cols.append("is_gene")
    elif "codeword_category" in schema_names:
        cols.append("codeword_category")
    return [c for c in cols if c in schema_names]


def _as_str_array(col: pa.ChunkedArray) -> pa.Array:
    if pa.types.is_binary(col.type) or pa.types.is_large_binary(col.type):
        col = pc.cast(col, pa.string())
    return col.combine_chunks() if isinstance(col, pa.ChunkedArray) else col


def gene_mask(table: pa.Table) -> np.ndarray:
    """Boolean mask of real-gene molecules in a transcript table chunk.

    Uses ``is_gene`` when present, otherwise ``codeword_category`` (values
    ending in ``_gene``), otherwise drops the control feature-name prefixes in
    :data:`CONTROL_PREFIXES` (case-insensitive).
    """
    n = table.num_rows
    if "is_gene" in table.column_names:
        return np.asarray(table["is_gene"].to_numpy(zero_copy_only=False), dtype=bool)
    if "codeword_category" in table.column_names:
        cwc = _as_str_array(table["codeword_category"]).to_numpy(zero_copy_only=False)
        return np.char.endswith(cwc.astype(str), "_gene")
    names = _as_str_array(table["feature_name"]).to_numpy(zero_copy_only=False).astype(str)
    lowered = np.char.lower(names)
    keep = np.ones(n, dtype=bool)
    for prefix in CONTROL_PREFIXES:
        keep &= ~np.char.startswith(lowered, prefix)
    return keep


def filter_mask(table: pa.Table, min_qv: float) -> np.ndarray:
    """Combined ``qv >= min_qv`` + real-gene mask."""
    qv = np.asarray(table["qv"].to_numpy(zero_copy_only=False))
    return (qv >= min_qv) & gene_mask(table)


def _bbox_mask(x: np.ndarray, y: np.ndarray, bbox: Sequence[float]) -> np.ndarray:
    x0, y0, x1, y1 = bbox
    return (x >= x0) & (x < x1) & (y >= y0) & (y < y1)


def read_transcript_crops(
    path: Path | str,
    bboxes: Sequence[Sequence[float]],
    min_qv: float,
    columns: Sequence[str] | None = None,
    with_z: bool = False,
) -> list[pa.Table]:
    """Read several crop boxes out of a transcripts table in one pass.

    Row groups are iterated one at a time; only the projected columns are
    materialised and only rows passing ``qv``/gene/bbox filters are kept, so
    multi-GB tables (Xenium Prime 5K) never need to fit in memory.
    ``with_z`` keeps the ``z_location`` column for 3D crops.
    """
    pf = pq.ParquetFile(str(path))
    cols = transcript_columns(pf.schema_arrow.names, with_z=with_z)
    final_cols = list(columns) if columns else [
        "x_location", "y_location", "qv", "feature_name", "cell_id", "overlaps_nucleus",
    ]
    if with_z:
        final_cols.append("z_location")
    final_cols = [c for c in final_cols if c in pf.schema_arrow.names]
    parts: list[list[pa.Table]] = [[] for _ in bboxes]
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=cols)
        if chunk.num_rows == 0:
            continue
        keep = filter_mask(chunk, min_qv)
        if not keep.any():
            continue
        chunk = chunk.filter(pa.array(keep))
        x = np.asarray(chunk["x_location"].to_numpy(zero_copy_only=False), dtype=np.float64)
        y = np.asarray(chunk["y_location"].to_numpy(zero_copy_only=False), dtype=np.float64)
        for i, bbox in enumerate(bboxes):
            m = _bbox_mask(x, y, bbox)
            if m.any():
                parts[i].append(chunk.filter(pa.array(m)))
    out = []
    for i, ps in enumerate(parts):
        if ps:
            t = pa.concat_tables([p.select(final_cols) for p in ps])
        else:
            schema = pa.schema(
                [pf.schema_arrow.field(c) for c in final_cols]
            )
            t = pa.Table.from_pylist([], schema=schema)
        out.append(t)
    return out


def assign_vendor_prior(
    cell_ids: np.ndarray, overlaps_nucleus: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Build ``cell_vendor`` (string) and ``prior`` (int32) columns.

    * ``cell_vendor``: the vendor cell id as a string, ``""`` when unassigned
      (``0``/negative for numeric ids, ``UNASSIGNED``/empty/NA for strings).
    * ``prior``: numeric vendor cell id when ``overlaps_nucleus == 1``, else 0.
      String ids are mapped deterministically to ``1..N`` in sorted order of
      the assigned ids that appear with nucleus overlap.
    """
    overlaps = np.asarray(overlaps_nucleus).astype(np.int8) == 1
    ids = np.asarray(cell_ids)
    if ids.dtype.kind in "iu":  # numeric ids, 0 = unassigned
        ids_i = ids.astype(np.int64)
        assigned = ids_i > 0
        vendor = np.where(assigned, ids_i.astype(str), "")
        prior = np.where(assigned & overlaps, ids_i, 0).astype(np.int32)
        return vendor, prior
    ids_s = ids.astype(str)
    unassigned = (ids_s == "") | (ids_s == "UNASSIGNED") | (ids_s == "nan") | (ids_s == "None")
    assigned = ~unassigned
    vendor = np.where(assigned, ids_s, "")
    sel = assigned & overlaps
    keys = np.unique(ids_s[sel])
    prior = np.zeros(len(ids_s), dtype=np.int32)
    if len(keys):
        prior[sel] = np.searchsorted(keys, ids_s[sel]).astype(np.int32) + 1
    return vendor, prior


def build_molecule_table(
    crop: pa.Table, bbox: Sequence[float] | None = None, *, keep_z: bool = False
) -> pa.Table:
    """Turn a filtered transcript crop into the contract molecule table.

    Columns: ``x``, ``y`` (float64, µm), ``z`` (float64, µm; only when
    ``keep_z`` and the crop carries ``z_location``), ``gene`` (string),
    ``qv`` (float32), ``prior`` (int32), ``cell_vendor`` (string).
    Sorted by (y, x).
    """
    x = np.asarray(crop["x_location"].to_numpy(zero_copy_only=False), dtype=np.float64)
    y = np.asarray(crop["y_location"].to_numpy(zero_copy_only=False), dtype=np.float64)
    if bbox is not None:  # defensive: reader already filtered
        m = _bbox_mask(x, y, bbox)
        x, y = x[m], y[m]
        crop = crop.filter(pa.array(m)) if not m.all() else crop
    qv = np.asarray(crop["qv"].to_numpy(zero_copy_only=False), dtype=np.float32)
    gene = _as_str_array(crop["feature_name"]).to_numpy(zero_copy_only=False).astype(str)
    vendor, prior = assign_vendor_prior(
        crop["cell_id"].to_numpy(zero_copy_only=False),
        crop["overlaps_nucleus"].to_numpy(zero_copy_only=False),
    )
    data = {
        "x": pa.array(x, pa.float64()),
        "y": pa.array(y, pa.float64()),
    }
    if keep_z and "z_location" in crop.column_names:
        z = np.asarray(crop["z_location"].to_numpy(zero_copy_only=False),
                       dtype=np.float64)
        data["z"] = pa.array(z, pa.float64())
    data.update({
        "gene": pa.array(gene, pa.string()),
        "qv": pa.array(qv, pa.float32()),
        "prior": pa.array(prior, pa.int32()),
        "cell_vendor": pa.array(vendor, pa.string()),
    })
    table = pa.table(data)
    idx = pc.sort_indices(table, sort_keys=[("y", "ascending"), ("x", "ascending")])
    return table.take(idx)


# ---------------------------------------------------------------------------
# Density histogram and crop selection
# ---------------------------------------------------------------------------


def transcript_histogram(
    path: Path | str,
    min_qv: float,
    bounds: Sequence[float],
    bin_um: float = 25.0,
) -> np.ndarray:
    """2-D histogram of filtered molecules over ``bounds`` (x0, y0, x1, y1)."""
    x0, y0, x1, y1 = bounds
    nx = int(math.ceil((x1 - x0) / bin_um))
    ny = int(math.ceil((y1 - y0) / bin_um))
    hist = np.zeros((ny, nx), dtype=np.int64)
    pf = pq.ParquetFile(str(path))
    schema_names = pf.schema_arrow.names
    cols = transcript_columns(schema_names)
    # The gene column is only needed for filtering when no is_gene/cwc exists.
    if "is_gene" in schema_names or "codeword_category" in schema_names:
        cols = [c for c in cols if c != "feature_name"]
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=cols)
        if chunk.num_rows == 0:
            continue
        keep = filter_mask(chunk, min_qv)
        if not keep.any():
            continue
        x = np.asarray(chunk["x_location"].to_numpy(zero_copy_only=False), dtype=np.float64)[keep]
        y = np.asarray(chunk["y_location"].to_numpy(zero_copy_only=False), dtype=np.float64)[keep]
        ix = ((x - x0) / bin_um).astype(np.int64)
        iy = ((y - y0) / bin_um).astype(np.int64)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        hist += np.bincount(iy[ok] * nx + ix[ok], minlength=ny * nx).reshape(ny, nx)
    return hist


def cells_histogram(cells, bounds: Sequence[float], bin_um: float = 25.0,
                    mask: np.ndarray | None = None) -> np.ndarray:
    """2-D histogram of vendor cell centroids over ``bounds``.

    ``mask`` optionally restricts the counted cells (e.g. one vendor cluster).
    """
    x0, y0, x1, y1 = bounds
    nx = int(math.ceil((x1 - x0) / bin_um))
    ny = int(math.ceil((y1 - y0) / bin_um))
    x = np.asarray(cells["x_centroid"], dtype=np.float64)
    y = np.asarray(cells["y_centroid"], dtype=np.float64)
    if mask is not None:
        x, y = x[mask], y[mask]
    ix = ((x - x0) / bin_um).astype(np.int64)
    iy = ((y - y0) / bin_um).astype(np.int64)
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    return np.bincount(iy[ok] * nx + ix[ok], minlength=ny * nx).reshape(ny, nx)


def _integral(a: np.ndarray) -> np.ndarray:
    """Zero-padded 2-D integral image; box sums via 4 corner lookups."""
    return np.pad(np.cumsum(np.cumsum(a, axis=0), axis=1), ((1, 0), (1, 0)))


def _box_sums(ii: np.ndarray, h: int, w: int) -> np.ndarray:
    """Sums of all h×w windows of the original array (from its integral image)."""
    return ii[h:, w:] - ii[:-h, w:] - ii[h:, :-w] + ii[:-h, :-w]


def density_class(cells_per_mm2: float) -> str:
    """Contract density class: sparse < 2500 <= medium <= 7000 < dense."""
    if cells_per_mm2 < 2500:
        return "sparse"
    if cells_per_mm2 <= 7000:
        return "medium"
    return "dense"


def gene_panel_class(n_panel_genes: int) -> str:
    """Contract gene-panel class from the panel size."""
    if n_panel_genes < 50:
        return "tiny"
    if n_panel_genes < 250:
        return "small"
    if n_panel_genes < 700:
        return "medium"
    if n_panel_genes < 2000:
        return "large"
    return "huge"


def pick_bbox(
    mol_hist: np.ndarray,
    cell_hist: np.ndarray,
    bin_um: float,
    bounds: Sequence[float],
    *,
    target: float,
    min_mols: float,
    max_mols: float,
    min_cells: int = 150,
    min_coverage: float = 0.7,
    density_hint: str = "any",
    dense_threshold: float = 7000.0,
    stroma_range: tuple[float, float] = (1000.0, 7000.0),
    within: Sequence[float] | None = None,
    seed: int = 0,
    cluster_hists: Sequence[np.ndarray] | None = None,
    criteria: dict | None = None,
) -> dict:
    """Deterministically choose a square crop box on the density histograms.

    Candidate boxes are all axis-aligned squares (in histogram bins) whose
    molecule count lies in ``[min_mols, max_mols]``, that contain at least
    ``min_cells`` vendor cells and at least ``min_coverage`` occupied bins
    (this is how "no large empty areas" is enforced).  Among those the box
    whose molecule count is closest to ``target`` (log-scale) wins; exact ties
    are broken with a seeded RNG so runs are reproducible.

    ``density_hint`` restricts the vendor-cell density of the box:
    ``any``, ``dense`` (>= ``dense_threshold`` cells/mm², relaxed to the
    densest qualifying quartile if impossible) or ``stromal`` (inside
    ``stroma_range`` with high tissue coverage).

    ``criteria`` carries the explicit crop-selection criteria documented in
    ``benchmarks/fetch/README.md`` (applied to *new* crops only):

    * ``min_clusters``: at least this many distinct vendor clusters must be
      present in the box; needs ``cluster_hists`` - one cells histogram per
      vendor cluster (from the bundle's ``analysis/.../clusters.csv``).
    * ``max_coverage``: at most this fraction of the box's bins may be
      occupied - how "tissue edge or empty space" crops are selected.
    * ``min_empty_border``: at least this fraction of the box's 2-bin-wide
      border ring must be empty (the empty space must touch the crop edge,
      i.e. a real tissue edge, not an interior hole).

    Returns a dict with ``bbox_um`` and per-box pick statistics (including
    the measured ``n_clusters``/``empty_border`` when the criteria apply).
    """
    if density_hint not in ("any", "dense", "stromal"):
        raise ValueError(f"unknown density_hint {density_hint!r}")
    criteria = dict(criteria or {})
    min_clusters = int(criteria.get("min_clusters") or 0)
    max_coverage = criteria.get("max_coverage")
    min_empty_border = criteria.get("min_empty_border")
    if min_clusters:
        if not cluster_hists or len(cluster_hists) < min_clusters:
            have = 0 if not cluster_hists else len(cluster_hists)
            raise RuntimeError(
                f"min_clusters={min_clusters} requested but only {have} "
                f"vendor cluster histograms available"
            )
    x0, y0, x1, y1 = bounds
    ny, nx = mol_hist.shape
    cov_hist = (mol_hist > 0).astype(np.int64)
    ii_m, ii_c, ii_v = _integral(mol_hist), _integral(cell_hist), _integral(cov_hist)
    ii_k = [_integral(h) for h in cluster_hists] if min_clusters else []

    # Restrict candidate top-left corners when a containing box is requested.
    if within is not None:
        wx0, wy0, wx1, wy1 = within
        col_lo = max(0, int(math.ceil((wx0 - x0) / bin_um)))
        col_hi = min(nx, int(math.floor((wx1 - x0) / bin_um)))
        row_lo = max(0, int(math.ceil((wy0 - y0) / bin_um)))
        row_hi = min(ny, int(math.floor((wy1 - y0) / bin_um)))
    else:
        col_lo, col_hi, row_lo, row_hi = 0, nx, 0, ny

    total_area = mol_hist.size * bin_um * bin_um
    mean_density = max(mol_hist.sum() / total_area, 1e-6)
    side0 = math.sqrt(target / mean_density)  # µm
    sizes = np.unique(
        np.round(np.arange(0.35, 3.0, 0.05) * side0 / bin_um).astype(np.int64)
    )
    sizes = sizes[(sizes >= 2) & (sizes <= min(row_hi - row_lo, col_hi - col_lo))]

    def collect(hint: str) -> list[dict]:
        candidates: list[dict] = []
        for s0 in sizes:
            s = int(s0)
            sums_m = _box_sums(ii_m, s, s)
            sums_c = _box_sums(ii_c, s, s)
            sums_v = _box_sums(ii_v, s, s)
            # window top-left (r, c) needs r + s <= row_hi, c + s <= col_hi
            r_lo, r_hi = row_lo, min(row_hi - s, ny - s) + 1
            c_lo, c_hi = col_lo, min(col_hi - s, nx - s) + 1
            if r_lo >= r_hi or c_lo >= c_hi:
                continue
            sm = sums_m[r_lo:r_hi, c_lo:c_hi]
            sc = sums_c[r_lo:r_hi, c_lo:c_hi]
            sv = sums_v[r_lo:r_hi, c_lo:c_hi]
            area = (s * bin_um) ** 2
            density = sc * 1e6 / area
            coverage = sv / float(s * s)
            ok = (
                (sm >= min_mols)
                & (sm <= max_mols)
                & (sc >= min_cells)
                & (coverage >= min_coverage)
            )
            if hint == "dense":
                ok &= density >= dense_threshold
            elif hint == "stromal":
                ok &= (density >= stroma_range[0]) & (density <= stroma_range[1])
            # explicit crop-selection criteria (never relaxed)
            ncl = None
            if ii_k:
                ncl = np.zeros(ok.shape, dtype=np.int32)
                for ii_clus in ii_k:
                    ncl += (_box_sums(ii_clus, s, s)[r_lo:r_hi, c_lo:c_hi] > 0)
                ok &= ncl >= min_clusters
            if max_coverage is not None:
                ok &= coverage <= max_coverage
            empty_border = None
            if min_empty_border is not None:
                if s > 4:
                    inner = _box_sums(ii_v, s - 4, s - 4)
                    ring = sv - inner[r_lo + 2:r_hi + 2, c_lo + 2:c_hi + 2]
                    ring_bins = s * s - (s - 4) * (s - 4)
                    empty_border = 1.0 - ring / float(ring_bins)
                    ok &= empty_border >= min_empty_border
                else:
                    ok = np.zeros_like(ok)  # border ring needs s > 4 bins
            if not ok.any():
                continue
            rr, cc = np.nonzero(ok)
            for r, c in zip(rr.tolist(), cc.tolist()):
                candidates.append(
                    {
                        "bbox_idx": (r + r_lo, c + c_lo, s),
                        "mols": int(sm[r, c]),
                        "cells": int(sc[r, c]),
                        "coverage": float(coverage[r, c]),
                        "density": float(density[r, c]),
                        "size_um": float(s * bin_um),
                        "n_clusters": (int(ncl[r, c]) if ncl is not None
                                       else None),
                        "empty_border": (round(float(empty_border[r, c]), 4)
                                          if empty_border is not None else None),
                    }
                )
        return candidates

    rng = np.random.default_rng(seed)
    candidates = collect(density_hint)
    relaxed = False
    if not candidates and density_hint != "any":
        # relax the density constraint: keep boxes meeting everything else
        # and take the densest quartile (recorded in the pick stats).
        relaxed = True
        candidates = collect("any")
        if candidates:
            dens = np.array([c["density"] for c in candidates])
            thr = np.quantile(dens, 0.75)
            pool = [c for c in candidates if c["density"] >= thr]
            if pool:
                candidates = pool

    if not candidates:
        raise RuntimeError(
            f"no crop box found: target={target}, hint={density_hint}, "
            f"criteria={criteria}, bounds={bounds}; relax min_mols/min_cells/"
            f"min_coverage"
        )
    res = _select(candidates, target, bin_um, bounds, rng)
    res["density_hint"] = density_hint
    res["relaxed_density"] = relaxed
    return res


def _select(
    cands: list[dict], target: float, bin_um: float, bounds: Sequence[float], rng
) -> dict:
    scores = np.array([-abs(math.log(c["mols"] / target)) for c in cands])
    best = scores.max()
    tied = [i for i, s in enumerate(scores) if s >= best - 1e-12]
    chosen = cands[tied[int(rng.integers(len(tied)))]]
    x0, y0 = bounds[0], bounds[1]
    r, c, s = chosen["bbox_idx"]
    bbox = [
        round(x0 + c * bin_um, 3),
        round(y0 + r * bin_um, 3),
        round(x0 + (c + s) * bin_um, 3),
        round(y0 + (r + s) * bin_um, 3),
    ]
    return {
        "bbox_um": bbox,
        "mols": chosen["mols"],
        "cells": chosen["cells"],
        "coverage": round(chosen["coverage"], 4),
        "cells_per_mm2": round(chosen["density"], 1),
        "size_um": chosen["size_um"],
        **({"n_clusters": chosen["n_clusters"]}
           if chosen.get("n_clusters") is not None else {}),
        **({"empty_border": chosen["empty_border"]}
           if chosen.get("empty_border") is not None else {}),
    }


# ---------------------------------------------------------------------------
# Cells, scale, boundaries
# ---------------------------------------------------------------------------


def cells_in_bbox(cells, bbox: Sequence[float]):
    x0, y0, x1, y1 = bbox
    x = np.asarray(cells["x_centroid"])
    y = np.asarray(cells["y_centroid"])
    m = (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
    return cells[m] if hasattr(cells, "__getitem__") else m


def load_cell_clusters(clusters_path: Path | str, cells) -> np.ndarray:
    """Vendor cluster id per row of ``cells`` from an analysis clusters.csv.

    The 10x Xenium bundles ship ``analysis/clustering/*/clusters.csv`` with
    ``Barcode,Cluster`` columns, where the barcode is the ``cell_id`` of
    ``cells.parquet``.  Cells absent from the file get ``-1``.
    """
    import pandas as pd

    cl = pd.read_csv(clusters_path)
    if cl.shape[1] < 2:
        raise ValueError(f"{clusters_path}: expected Barcode,Cluster columns")
    barcodes, labels = cl.iloc[:, 0], cl.iloc[:, 1]
    ids_num = pd.to_numeric(cells["cell_id"], errors="coerce")
    keys_num = pd.to_numeric(barcodes, errors="coerce")
    if ids_num.notna().all() and keys_num.notna().all():
        mapping = dict(zip(keys_num.astype(np.int64), labels))
        out = ids_num.astype(np.int64).map(mapping)
    else:
        mapping = dict(zip(barcodes.astype(str), labels))
        out = cells["cell_id"].astype(str).map(mapping)
    return out.fillna(-1).to_numpy(dtype=np.int64)


def estimate_scale_um(cells, bbox: Sequence[float]) -> tuple[float, str]:
    """Baysor ``scale_um`` from vendor cell/nucleus areas.

    Primary method: 1.5 x the equivalent radius of the median vendor nucleus
    area (``cells.parquet.nucleus_area`` of cells whose centroid lies in the
    crop).  Falls back to the median cell radius when nucleus areas are
    missing.  The returned string records the method and the medians used.
    """
    sel = cells_in_bbox(cells, bbox)
    nuc = np.asarray(sel["nucleus_area"], dtype=np.float64)
    nuc = nuc[np.isfinite(nuc) & (nuc > 0)]
    if len(nuc) >= 30:
        med = float(np.median(nuc))
        radius = math.sqrt(med / math.pi)
        method = (
            f"1.5 * sqrt(median(nucleus_area)/pi); median nucleus_area={med:.1f} um^2 "
            f"from {len(nuc)} vendor cells in crop (cells.parquet)"
        )
    else:
        area = np.asarray(sel["cell_area"], dtype=np.float64)
        area = area[np.isfinite(area) & (area > 0)]
        med = float(np.median(area))
        radius = math.sqrt(med / math.pi)
        method = (
            f"median cell equivalent radius sqrt(median(cell_area)/pi); "
            f"median cell_area={med:.1f} um^2 from {len(area)} vendor cells in crop "
            f"(nucleus_area unavailable)"
        )
    return round(1.5 * radius, 2), method


def crop_boundaries(path: Path | str, bbox: Sequence[float]) -> pa.Table:
    """Crop a Xenium boundaries parquet (object_id, vertex_x, vertex_y) box.

    A boundary object is kept whole when its bounding box intersects the crop
    box; vertices outside the box are retained so the polygon stays closed.
    """
    x0, y0, x1, y1 = bbox
    t = pq.read_table(str(path))
    cid = t["cell_id"].to_numpy(zero_copy_only=False)
    vx = np.asarray(t["vertex_x"].to_numpy(zero_copy_only=False), dtype=np.float64)
    vy = np.asarray(t["vertex_y"].to_numpy(zero_copy_only=False), dtype=np.float64)
    if len(cid) == 0:
        return t
    starts = np.flatnonzero(np.concatenate(([True], cid[1:] != cid[:-1])))
    lens = np.diff(np.append(starts, len(cid)))
    xmin = np.minimum.reduceat(vx, starts)
    xmax = np.maximum.reduceat(vx, starts)
    ymin = np.minimum.reduceat(vy, starts)
    ymax = np.maximum.reduceat(vy, starts)
    keep_run = (xmin <= x1) & (xmax >= x0) & (ymin <= y1) & (ymax >= y0)
    row_keep = np.repeat(keep_run, lens)
    return t.filter(pa.array(row_keep))


def panel_gene_count(path: Path | str) -> int:
    """Number of real genes in a Xenium ``gene_panel.json`` (descriptor gene)."""
    with open(path) as fh:
        panel = json.load(fh)
    targets = panel.get("payload", {}).get("targets", [])
    genes = sum(1 for t in targets if t.get("type", {}).get("descriptor") == "gene")
    if genes:  # prefer descriptor-based count when available
        return genes
    # older/simple panels: count targets that are not control sources
    return sum(
        1
        for t in targets
        if t.get("source", {}).get("category", "current") in ("current", "base", "custom")
    )


# ---------------------------------------------------------------------------
# Focus images
# ---------------------------------------------------------------------------


def ome_info(path: Path | str) -> dict:
    """Pixel size, shape and channel names of an OME-TIFF focus image.

    Xenium "morphology_focus" outputs are *multi-file OME*: each channel lives
    in its own ``morphology_focus_NNNN.ome.tif`` (C=0/DAPI -> 0000,
    C=1/boundary stain -> 0001, ...) and every file's OME XML references all
    of them.  tifffile resolves the sibling files by relative path, so any
    member of the set can be opened to read every channel.
    """
    import tifffile

    with tifffile.TiffFile(str(path)) as tf:
        series = tf.series[0]
        md = tf.ome_metadata or ""
        px = None
        import re

        m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', md)
        if m:
            px = float(m.group(1))
        channels = re.findall(r'<Channel .*?Name="([^"]*)"', md)
    return {
        "pixel_size_um": px,
        "shape": tuple(series.shape),  # (C, Y, X)
        "axes": series.axes,
        "channels": channels,
    }


def _open_focus_array(path: Path | str):
    import tifffile
    import zarr

    tf = tifffile.TiffFile(str(path))
    obj = zarr.open(tf.series[0].aszarr(), mode="r")
    if isinstance(obj, zarr.Group):  # pyramidal OME-TIFF: full resolution is "0"
        obj = obj["0"]
    return tf, obj


def read_image_window(
    path: Path | str, channel: int, bbox: Sequence[float], pixel_size_um: float
) -> np.ndarray:
    """Read only the crop window of one channel of a focus OME-TIFF.

    The image origin is (0, 0) µm and pixel rows grow with y (verified against
    the vendor cell extents).  The window is clamped to the image; the caller
    should recompute the origin from the returned shape if clamping mattered.
    """
    tf, arr = _open_focus_array(path)
    try:
        if arr.ndim == 2:  # single-channel image
            height, width = arr.shape
            if channel != 0:
                raise IndexError(f"channel {channel} requested from 2-D image")
            plane = arr
        else:
            _, height, width = arr.shape
            plane = arr[channel]
        x0, y0, x1, y1 = bbox
        c0 = max(0, int(math.floor(x0 / pixel_size_um)))
        c1 = min(width, int(math.ceil(x1 / pixel_size_um)))
        r0 = max(0, int(math.floor(y0 / pixel_size_um)))
        r1 = min(height, int(math.ceil(y1 / pixel_size_um)))
        return np.asarray(plane[r0:r1, c0:c1])
    finally:
        tf.close()


def write_tif(path: Path | str, image: np.ndarray) -> None:
    """Write a single-channel crop image (uint16/uint8, zlib-compressed)."""
    import tifffile

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        str(path),
        image,
        photometric="minisblack",
        compression="zlib",
        metadata={"axes": "YX"},
    )


def rasterize_nucleus_labels(
    path: Path | str, bbox: Sequence[float], out_path: Path | str
) -> dict:
    """Rasterise vendor nucleus boundaries into a Baysor image-prior label TIFF.

    Baysor maps the molecule at ``(x, y)`` to pixel
    ``(round(x) - 1, round(y) - 1)`` (:func:`load_prior_from_image` in
    ``src/data_loading/prior_segmentation.cpp``), i.e. the label raster must
    be a **1 µm/pixel grid aligned with the absolute coordinate origin**:
    pixel ``(c, r)`` covers ``x in [c+0.5, c+1.5)``, ``y in [r+0.5, r+1.5)``,
    and the image spans ``[0, ceil(x1)) x [0, ceil(y1))`` pixels of the crop
    box.  Every polygon/pixel intersection is tested exactly (shapely), so
    any molecule inside a nucleus always lands on that nucleus's pixel.

    Labels are dense ``1..K`` in sorted order of the boundary ``cell_id``
    values intersecting the crop (``0`` = background, uint16), which keeps the
    TIFF small; the mapping is deterministic.  Returns the image spec that
    goes into ``meta.json``.
    """
    import math

    import shapely

    x0, y0, x1, y1 = (float(v) for v in bbox)
    if not (x1 > x0 and y1 > y0):
        raise ValueError(f"degenerate bbox {bbox}")
    width, height = math.ceil(x1), math.ceil(y1)

    t = crop_boundaries(path, bbox)
    if t.num_rows == 0:
        raise ValueError(f"no boundary objects intersect bbox {bbox} in {path}")
    cid = t["cell_id"].to_numpy(zero_copy_only=False)
    vx = np.asarray(t["vertex_x"].to_numpy(zero_copy_only=False), dtype=np.float64)
    vy = np.asarray(t["vertex_y"].to_numpy(zero_copy_only=False), dtype=np.float64)

    # vertices come grouped per object; fall back to a stable sort if not
    starts = np.flatnonzero(np.concatenate(([True], cid[1:] != cid[:-1])))
    if len(starts) != len(np.unique(cid)):
        order = np.argsort(cid, kind="stable")
        cid, vx, vy = cid[order], vx[order], vy[order]
        starts = np.flatnonzero(np.concatenate(([True], cid[1:] != cid[:-1])))
    ends = np.append(starts[1:], len(cid))
    run_ids = cid[starts]
    keys = np.unique(run_ids)
    if len(keys) > 65535:
        raise ValueError(
            f"{len(keys)} boundary objects in crop; label TIFF is uint16 "
            f"(max 65535) - shrink the crop"
        )

    labels = np.zeros((height, width), dtype=np.uint16)
    n_polygons = 0
    for run_id, s, e in zip(run_ids, starts, ends):
        if e - s < 3:
            continue
        pts = np.column_stack([vx[s:e], vy[s:e]])
        poly = shapely.Polygon(pts)
        if poly.is_empty:
            continue
        if not poly.is_valid:
            poly = poly.buffer(0)
            if poly.is_empty:
                continue
        n_polygons += 1
        label = int(np.searchsorted(keys, run_id)) + 1
        xmin, ymin, xmax, ymax = poly.bounds
        c0 = max(0, int(math.ceil(xmin - 1.5)))
        c1 = min(width - 1, int(math.floor(xmax - 0.5)))
        r0 = max(0, int(math.ceil(ymin - 1.5)))
        r1 = min(height - 1, int(math.floor(ymax - 0.5)))
        if c0 > c1 or r0 > r1:
            continue
        cols = np.arange(c0, c1 + 1, dtype=np.float64)
        rows = np.arange(r0, r1 + 1, dtype=np.float64)
        cc, rr = np.meshgrid(cols + 1.0, rows + 1.0)  # pixel centers
        boxes = shapely.box(cc - 0.5, rr - 0.5, cc + 0.5, rr + 0.5)
        hit = shapely.intersects(poly, boxes).reshape(len(rows), len(cols))
        region = labels[r0:r1 + 1, c0:c1 + 1]
        region[hit] = np.where(region[hit] == 0, label, region[hit])

    if not n_polygons:
        raise ValueError(f"no usable polygon in {path} for bbox {bbox}")
    # Uncompressed on purpose: Baysor reads the label TIFF row by row starting
    # at the first row of the molecule window (usually mid-strip), and this
    # libtiff build refuses sub-strip scanline reads of Deflate data with
    # "Compression algorithm does not support random access".  Uncompressed
    # strips seek arbitrarily; a full-section frame is still only ~30 MB.
    import tifffile

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        str(out_path),
        labels,
        photometric="minisblack",
        metadata={"axes": "YX"},
    )
    return {
        "file": str(out_path),
        "width": width,
        "height": height,
        "pixel_size_um": 1.0,
        "origin_um": [0.0, 0.0],
        "n_labels": int(len(keys)),
        "n_polygons": n_polygons,
        "n_labelled_pixels": int(np.count_nonzero(labels)),
        "label_source": "dense vendor cell_id mapping (sorted), 0 = background",
    }


# ---------------------------------------------------------------------------
# meta.json and reporting
# ---------------------------------------------------------------------------


def molecule_stats(table: pa.Table, bbox: Sequence[float], n_vendor_cells: int) -> dict:
    """Contract ``stats`` block from a built molecule table."""
    x0, y0, x1, y1 = bbox
    area = (x1 - x0) * (y1 - y0)
    n = table.num_rows
    return {
        "n_molecules": int(n),
        "n_genes": len(np.unique(table["gene"].to_numpy(zero_copy_only=False).astype(str))),
        "area_um2": round(area, 1),
        "molecules_per_um2": round(n / area, 4) if area else 0.0,
        "n_vendor_cells": int(n_vendor_cells),
        "vendor_cells_per_mm2": round(n_vendor_cells * 1e6 / area, 1) if area else 0.0,
    }


def build_meta(
    *,
    dataset_id: str,
    tier: str,
    source: dict,
    bbox: Sequence[float],
    note: str,
    stats: dict,
    panel_genes: int,
    scale_um: float,
    scale_method: str,
    baysor: dict,
    images: list[dict],
    retrieved: str,
    difficulty_notes: str = "",
    z_range: Sequence[float] | None = None,
) -> dict:
    """Assemble the contract ``meta.json`` structure.

    ``baysor`` carries the manifest defaults (config, prior settings, ...);
    ``scale_um``/``scale_method`` are filled in from the vendor cell areas.
    ``z_range`` (min, max in µm) is recorded for 3D crops that keep z.
    """
    density = density_class(stats["vendor_cells_per_mm2"])
    baysor_block = dict(baysor)
    baysor_block["scale_um"] = scale_um
    baysor_block["scale_um_method"] = scale_method
    return {
        "id": dataset_id,
        "kind": "real",
        "tier": tier,
        "platform": "Xenium",
        "source": {
            "url": source["url"],
            "doi": source.get("doi"),
            "license": source.get("license"),
            "original_dataset": source.get("original_dataset"),
            "retrieved": retrieved,
        },
        "crop": {
            "bbox_um": [round(float(v), 3) for v in bbox],
            "z_range_um": ([round(float(v), 3) for v in z_range]
                           if z_range is not None else None),
            "note": note,
        },
        "stats": stats,
        "difficulty": {
            "cell_density": density,
            "gene_panel": gene_panel_class(panel_genes),
            "notes": difficulty_notes or note,
        },
        "baysor": baysor_block,
        "images": images,
        "truth": None,
    }


def dataset_readme(meta: dict, extra_lines: Sequence[str] = ()) -> str:
    """Small provenance README written into each dataset directory."""
    src = meta["source"]
    st = meta["stats"]
    b = meta["baysor"]
    lines = [
        f"# {meta['id']}",
        "",
        f"Cropped from `{src['url']}`",
        f"(original dataset: {src['original_dataset']}; retrieved {src['retrieved']}).",
        "",
        f"* tier: {meta['tier']}; crop bbox (µm): {meta['crop']['bbox_um']}",
        f"* {st['n_molecules']} molecules, {st['n_genes']} genes, "
        f"{st['n_vendor_cells']} vendor cells "
        f"({st['vendor_cells_per_mm2']:.0f} cells/mm², "
        f"density class {meta['difficulty']['cell_density']})",
        f"* crop note: {meta['crop']['note']}",
        f"* scale_um: {b['scale_um']} ({b['scale_um_method']})",
        f"* prior: {b['prior']} (prior_confidence {b['prior_confidence']}), "
        f"config {b['config']}",
        "",
        "Regenerate with:",
        "",
        "```bash",
        "python benchmarks/fetch/xenium.py fetch " + meta["id"],
        "python benchmarks/fetch/xenium.py build " + meta["id"],
        "```",
        "",
    ]
    for extra in extra_lines:
        lines.append(extra)
    return "\n".join(lines) + "\n"
