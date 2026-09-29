#!/usr/bin/env python3
"""Contract validator for every benchmark dataset under ``sim/`` and ``real/``.

Checks the dataset contract defined in ``benchmarks/README.md``:

* **molecules.parquet** — required columns and dtypes, rows sorted by
  (y, x), finite coordinates, no control/blank gene names;
* **meta.json** — required fields, enums (kind, tier, density/panel classes,
  ``baysor.prior``, truth), tier molecule budgets, image/prior/config
  references, bbox sanity;
* **stats** — every recorded stat recomputed from the data (molecule and
  gene counts, area vs the crop bbox, densities, difficulty classes);
* **manifests** — dataset listed, tier/platform mirrored, ``sha256`` of
  ``molecules.parquet``/``meta.json`` verified where recorded, the manifest
  ``baysor`` block mirrored into ``meta.json``;
* **file references** — ``images[].file``, the ``baysor.prior`` image and
  ``baysor.config`` all exist; image frames cover the crop bbox.

Findings have two severities:

* ``error`` — a *metadata* problem (``meta.json``, the manifests, missing
  referenced files, inconsistent ``stats``): fix these;
* ``data`` — a *content* problem of ``molecules.parquet`` that can only be
  fixed by rebuilding the dataset (unsorted rows, wrong dtypes, control
  genes, vendor-label sentinels, molecules outside the recorded bbox...):
  these are reported, never silently rebuilt — baselines depend on the
  content hashes.

Exit code: 0 = no errors (``data`` findings allowed unless ``--strict``),
1 = errors found, 2 = usage error.

Usage::

    validate_datasets.py                       # all datasets
    validate_datasets.py --datasets quick,id   # tier spec / ids / globs
    validate_datasets.py --json report.json    # machine-readable report
    validate_datasets.py --strict              # data findings fail too
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

# --- contract enums / thresholds -----------------------------------------

KINDS = ("sim", "real")
TIERS = ("quick", "full")
QUICK_MAX_MOLECULES = 150_000        # "at most about 150k"
FULL_MAX_MOLECULES = 3_000_000       # "at most about 3M"
DENSITY_CLASSES = ("sparse", "medium", "dense")
PANEL_CLASSES = ("tiny", "small", "medium", "large", "huge")
META_REQUIRED = ("id", "kind", "tier", "platform", "source", "crop",
                 "stats", "difficulty", "baysor", "images", "truth")
REAL_STATS_REQUIRED = ("n_molecules", "n_genes", "area_um2",
                       "molecules_per_um2", "n_vendor_cells",
                       "vendor_cells_per_mm2")
SIM_STATS_REQUIRED = REAL_STATS_REQUIRED + ("n_true_cells",
                                            "true_cells_per_mm2")
CONTROL_GENE_PREFIXES = ("NegControl", "NegCtrl", "NegProbe", "NegPrb",
                         "Negative", "SystemControl", "System Control",
                         "UnassignedCodeword", "DeprecatedCodeword",
                         "Intergenic_Region", "BLANK")
VENDOR_SENTINELS = {"0", "na", "nan", "null", "none", "unassigned", "cell_0"}
OUTSIDE_BBOX_TOL_UM = 0.1     # coordinate rounding tolerance
Z_RANGE_TOL_UM = 0.1
IMAGE_FRAME_TOL_PX = 2        # origin may sit up to 2 pixels off the bbox

MANIFESTS = ("sim.yaml", "real_xenium.yaml", "real_other.yaml")


@dataclass
class Issue:
    severity: str              # "error" (metadata) | "data"
    dataset: str
    check: str
    message: str


@dataclass
class Report:
    datasets: int = 0
    errors: int = 0
    data_findings: int = 0
    issues: list[Issue] = field(default_factory=list)

    def add(self, severity: str, dataset: str, check: str, message: str) -> None:
        self.issues.append(Issue(severity, dataset, check, message))
        if severity == "error":
            self.errors += 1
        else:
            self.data_findings += 1


# --- helpers ---------------------------------------------------------------

def gene_panel_class(n_genes: int) -> str:
    """Contract gene-panel class for a gene count (boundary cases follow the
    documented ranges: tiny < 50, small 50-250, medium 250-700,
    large 700-2000, huge > 2000; a value on a shared boundary belongs to the
    lower class)."""
    if n_genes < 50:
        return "tiny"
    if n_genes <= 250:
        return "small"
    if n_genes <= 700:
        return "medium"
    if n_genes <= 2000:
        return "large"
    return "huge"


def cell_density_class(cells_per_mm2: float) -> str:
    """Contract density class: sparse < 2500, medium 2500-7000, dense > 7000."""
    if cells_per_mm2 < 2500:
        return "sparse"
    if cells_per_mm2 <= 7000:
        return "medium"
    return "dense"


def load_manifest_index(repo: Path) -> dict[str, dict]:
    """Index every dataset entry of ``benchmarks/datasets/*.yaml`` by id.

    Returns ``{id: {"group", "path", "entry"}}``.
    """
    idx: dict[str, dict] = {}
    for name in MANIFESTS:
        path = repo / "benchmarks" / "datasets" / name
        if not path.is_file():
            continue
        with open(path) as fh:
            doc = yaml.safe_load(fh) or {}
        defaults = (doc.get("defaults") or {}).get("baysor") or {}
        for entry in doc.get("datasets") or []:
            eid = entry.get("id")
            if not eid:
                continue
            # the baysor block a dataset's meta.json must mirror: the
            # manifest defaults merged with per-dataset overrides
            baysor = {**defaults, **(entry.get("baysor") or {})}
            idx[eid] = {"group": path.stem, "path": path, "entry": entry,
                        "baysor": baysor}
    return idx


def _is_str_dict(t) -> bool:
    return str(t).startswith("dictionary<values=string")


def _close(a: float, b: float, abs_tol: float, rel_tol: float) -> bool:
    return abs(a - b) <= max(abs_tol, rel_tol * max(abs(b), 1e-30))


def _sha256(path: Path) -> str:
    return common.sha256_file(path)


def _discover_dirs(root: Path) -> list[Path]:
    out: list[Path] = []
    for kind in ("sim", "real"):
        base = root / kind
        if base.is_dir():
            out += sorted(d for d in base.iterdir() if d.is_dir())
    return out


def _select(dirs: Iterable[Path], spec: Optional[str], root: Path) -> list[Path]:
    dirs = list(dirs)
    if not spec:
        return dirs
    spec = spec.strip()
    if spec in common.TIERS:
        tier = spec
        out = []
        for d in dirs:
            mp = d / "meta.json"
            if mp.is_file():
                try:
                    if common.read_json(mp).get("tier") == tier:
                        out.append(d)
                except (OSError, ValueError):
                    continue
        return out
    tokens = [t.strip() for t in spec.split(",") if t.strip()]
    out, seen = [], set()
    for tok in tokens:
        matched = [d for d in dirs if d.name == tok]
        if not matched:
            matched = [d for d in dirs if fnmatch.fnmatch(d.name, tok)]
        if not matched:
            raise ValueError(f"no dataset matches '{tok!r}' under {root}")
        for d in matched:
            if d not in seen:
                seen.add(d)
                out.append(d)
    return out


# --- per-dataset validation ------------------------------------------------

def validate_dataset(ds_dir: Path, repo: Path, manifests: dict[str, dict],
                     report: Report, check_hashes: bool = True) -> None:
    """Validate one dataset directory; appends findings to ``report``."""
    ds_id = ds_dir.name
    meta_path = ds_dir / "meta.json"
    mol_path = ds_dir / "molecules.parquet"

    if not meta_path.is_file():
        report.add("error", ds_id, "files", "meta.json missing")
        return
    if not mol_path.is_file():
        report.add("error", ds_id, "files", "molecules.parquet missing")
        return
    try:
        meta = common.read_json(meta_path)
    except (OSError, ValueError) as exc:
        report.add("error", ds_id, "meta", f"meta.json is not valid JSON: {exc}")
        return
    if not isinstance(meta, dict):
        report.add("error", ds_id, "meta", "meta.json is not a JSON object")
        return

    _validate_meta_shape(ds_id, meta, ds_dir, report)
    kind = meta.get("kind")

    if not pq.ParquetFile(mol_path).metadata:
        report.add("error", ds_id, "molecules", "unreadable parquet file")
        return
    schema = pq.read_schema(mol_path)
    _validate_columns(ds_id, kind, schema, report)

    try:
        df = pd.read_parquet(mol_path)
    except Exception as exc:                                    # noqa: BLE001
        report.add("data", ds_id, "molecules", f"cannot read table: {exc}")
        return

    _validate_content(ds_id, kind, df, report)
    if isinstance(meta.get("crop"), dict):
        _validate_geometry(ds_id, kind, meta, df, report)
    if isinstance(meta.get("stats"), dict):
        _validate_stats(ds_id, kind, meta, df, report)
    if isinstance(meta.get("difficulty"), dict) and isinstance(meta.get("stats"), dict):
        _validate_difficulty(ds_id, kind, meta, report)
    if isinstance(meta.get("baysor"), dict):
        _validate_baysor(ds_id, ds_dir, repo, kind, meta, schema, report)
    if isinstance(meta.get("images"), list):
        _validate_images(ds_id, ds_dir, meta, report)
    _validate_manifest(ds_id, meta, repo, manifests, mol_path, meta_path,
                       report, check_hashes=check_hashes)


def _validate_meta_shape(ds_id: str, meta: dict, ds_dir: Path,
                         report: Report) -> None:
    for key in META_REQUIRED:
        if key not in meta:
            report.add("error", ds_id, "meta", f"required field '{key}' missing")
    if meta.get("id") != ds_dir.name:
        report.add("error", ds_id, "meta",
                   f"id {meta.get('id')!r} does not match directory {ds_dir.name!r}")
    if meta.get("kind") not in KINDS:
        report.add("error", ds_id, "meta", f"kind {meta.get('kind')!r} not in {KINDS}")
    elif ds_dir.parent.name != meta["kind"]:
        report.add("error", ds_id, "meta",
                   f"kind {meta['kind']!r} but stored under {ds_dir.parent.name}/")
    tier = meta.get("tier")
    if tier not in TIERS:
        report.add("error", ds_id, "meta", f"tier {tier!r} not in {TIERS}")
    platform = meta.get("platform")
    if not isinstance(platform, str) or not platform.strip():
        report.add("error", ds_id, "meta", "platform must be a non-empty string")
    elif meta.get("kind") == "sim" and platform != "simulated":
        report.add("error", ds_id, "meta",
                   f"sim datasets must have platform 'simulated', got {platform!r}")

    source = meta.get("source")
    if not isinstance(source, dict):
        report.add("error", ds_id, "meta", "source must be an object")
    elif meta.get("kind") == "real":
        for key in ("url", "license", "original_dataset", "retrieved", "doi"):
            if key not in source:
                report.add("error", ds_id, "meta", f"source.{key} missing (real)")
    elif meta.get("kind") == "sim":
        for key in ("generator", "generator_version", "seed"):
            if key not in source:
                report.add("error", ds_id, "meta", f"source.{key} missing (sim)")

    crop = meta.get("crop")
    if not isinstance(crop, dict):
        report.add("error", ds_id, "meta", "crop must be an object")
    else:
        bb = crop.get("bbox_um")
        if (not isinstance(bb, list) or len(bb) != 4
                or not all(isinstance(v, (int, float)) for v in bb)
                or not (bb[0] < bb[2] and bb[1] < bb[3])):
            report.add("error", ds_id, "meta",
                       f"crop.bbox_um must be [x0, y0, x1, y1] with x0<x1, y0<y1; got {bb!r}")
        zr = crop.get("z_range_um")
        if "z_range_um" not in crop:
            report.add("error", ds_id, "meta", "crop.z_range_um missing (may be null)")
        elif zr is not None and not (isinstance(zr, list) and len(zr) == 2
                                     and zr[0] <= zr[1]):
            report.add("error", ds_id, "meta", f"crop.z_range_um invalid: {zr!r}")
        if "note" not in crop:
            report.add("error", ds_id, "meta", "crop.note missing")

    truth = meta.get("truth")
    if meta.get("kind") == "sim":
        if not isinstance(truth, dict):
            report.add("error", ds_id, "meta", "sim dataset must have a truth object")
        else:
            for key in ("generator", "generator_version", "seed", "params"):
                if key not in truth:
                    report.add("error", ds_id, "meta", f"truth.{key} missing (sim)")
            for key in ("oracle_accuracy", "naive_accuracy"):
                v = truth.get(key)
                if v is not None and not (isinstance(v, (int, float)) and 0 <= v <= 1):
                    report.add("error", ds_id, "meta",
                               f"truth.{key} must be in [0, 1] or null, got {v!r}")
    elif truth is not None:
        report.add("error", ds_id, "meta", "real dataset must have truth: null")

    difficulty = meta.get("difficulty")
    if not isinstance(difficulty, dict):
        report.add("error", ds_id, "meta", "difficulty must be an object")
    else:
        if difficulty.get("cell_density") not in DENSITY_CLASSES:
            report.add("error", ds_id, "meta",
                       f"difficulty.cell_density {difficulty.get('cell_density')!r} "
                       f"not in {DENSITY_CLASSES}")
        if difficulty.get("gene_panel") not in PANEL_CLASSES:
            report.add("error", ds_id, "meta",
                       f"difficulty.gene_panel {difficulty.get('gene_panel')!r} "
                       f"not in {PANEL_CLASSES}")
        if "notes" not in difficulty:
            report.add("error", ds_id, "meta", "difficulty.notes missing")

    bcfg = meta.get("baysor")
    if not isinstance(bcfg, dict):
        report.add("error", ds_id, "meta", "baysor must be an object")
    if "images" not in meta:
        report.add("error", ds_id, "meta", "images missing (may be [])")


def _validate_columns(ds_id: str, kind: Optional[str], schema,
                      report: Report) -> None:
    names = schema.names
    types = {n: str(schema.field(n).type) for n in names}
    required = ["x", "y", "gene"]
    if kind == "sim":
        required.append("cell")
    if kind == "real":
        required.append("cell_vendor")
    for col in required:
        if col not in types:
            report.add("data", ds_id, "columns", f"required column '{col}' missing")

    expect = {
        "x": ("double",), "y": ("double",), "z": ("double",),
        "qv": ("float",), "prior": ("int32",), "cell": ("int32",),
        "interior": ("bool",),
    }
    for col, allowed in expect.items():
        if col in types and types[col] not in allowed:
            report.add("data", ds_id, "dtypes",
                       f"column '{col}' has type {types[col]}, expected {'/'.join(allowed)}")
    for col in ("gene", "cell_vendor", "celltype"):
        if col in types and not (types[col] == "string" or _is_str_dict(types[col])):
            report.add("data", ds_id, "dtypes",
                       f"column '{col}' has type {types[col]}, expected string")
    for col in names:
        if col not in expect and col not in ("gene", "cell_vendor", "celltype"):
            if not col.startswith("aux_"):
                report.add("data", ds_id, "columns",
                           f"unexpected column '{col}' (only contract columns and "
                           "aux_* extras are allowed)")


def _validate_content(ds_id: str, kind: Optional[str], df: pd.DataFrame,
                      report: Report) -> None:
    cols = set(df.columns)
    if {"x", "y"} <= cols:
        xy = df[["x", "y"]].to_numpy()
        if not np.isfinite(xy).all():
            report.add("data", ds_id, "content", "non-finite x/y coordinates")
        y = df["y"].to_numpy()
        x = df["x"].to_numpy()
        if len(y) > 1:
            y_ok = y[:-1] <= y[1:]
            x_ok = (y[:-1] != y[1:]) | (x[:-1] <= x[1:])
            if not (y_ok.all() and x_ok.all()):
                bad = int((~(y_ok & x_ok)).argmax())
                report.add("data", ds_id, "sorted",
                           f"rows are not sorted by (y, x); first violation at row {bad}")
    if "gene" in cols:
        g = df["gene"].astype(str)
        if df["gene"].isna().any() or (g.str.strip() == "").any():
            report.add("data", ds_id, "content", "null/empty gene values present")
        for prefix in CONTROL_GENE_PREFIXES:
            n = int(g.str.startswith(prefix).sum())
            if n:
                report.add("data", ds_id, "content",
                           f"{n} molecules carry control/blank gene prefix "
                           f"{prefix!r} (must be removed by the fetch script)")
                break
    if kind == "sim" and "cell" in cols:
        cells = df["cell"].to_numpy()
        if (cells < 0).any():
            report.add("data", ds_id, "content", "negative sim cell ids")
        if "interior" in cols and not df["interior"].astype(bool).any():
            report.add("data", ds_id, "content", "interior is false for every molecule")
    if "prior" in cols and (df["prior"].to_numpy() < 0).any():
        report.add("data", ds_id, "content", "negative prior labels (0 = no prior)")
    if kind == "real" and "cell_vendor" in cols:
        raw = df["cell_vendor"].fillna("").astype(str).str.strip()
        bad = sorted({v for v in raw.str.lower().unique() if v in VENDOR_SENTINELS})
        if bad:
            n = int(raw.str.lower().isin(bad).sum())
            report.add("data", ds_id, "content",
                       f"cell_vendor uses unassigned sentinel(s) {bad} on {n} rows "
                       "(the contract requires an empty string for unassigned)")


def _validate_geometry(ds_id: str, kind: Optional[str], meta: dict,
                       df: pd.DataFrame, report: Report) -> None:
    crop = meta["crop"]
    bb = crop.get("bbox_um")
    if isinstance(bb, list) and len(bb) == 4 and {"x", "y"} <= set(df.columns):
        x0, y0, x1, y1 = map(float, bb)
        x = df["x"].to_numpy()
        y = df["y"].to_numpy()
        dx = np.maximum(np.maximum(x0 - x, x - x1), 0.0)
        dy = np.maximum(np.maximum(y0 - y, y - y1), 0.0)
        outside = (dx > OUTSIDE_BBOX_TOL_UM) | (dy > OUTSIDE_BBOX_TOL_UM)
        n = int(outside.sum())
        if n:
            over = float(max(dx.max(), dy.max()))
            report.add("data", ds_id, "bbox",
                       f"{n} molecules lie more than {OUTSIDE_BBOX_TOL_UM} um "
                       f"outside crop.bbox_um (max overshoot {over:.3f} um)")
    zr = crop.get("z_range_um")
    if isinstance(zr, list) and "z" in df.columns:
        lo, hi = map(float, zr)
        z = df["z"].to_numpy()
        if len(z) and (z.min() < lo - Z_RANGE_TOL_UM or z.max() > hi + Z_RANGE_TOL_UM):
            report.add("data", ds_id, "z_range",
                       f"observed z range [{z.min():.3f}, {z.max():.3f}] exceeds "
                       f"crop.z_range_um [{lo}, {hi}] (tol {Z_RANGE_TOL_UM} um)")
    if zr is None and "z" in df.columns:
        report.add("error", ds_id, "z_range",
                   "z column present but crop.z_range_um is null")
    if isinstance(zr, list) and "z" not in df.columns:
        report.add("error", ds_id, "z_range",
                   "crop.z_range_um set but the table has no z column")


def _validate_stats(ds_id: str, kind: Optional[str], meta: dict,
                    df: pd.DataFrame, report: Report) -> None:
    stats = meta["stats"]
    required = SIM_STATS_REQUIRED if kind == "sim" else REAL_STATS_REQUIRED
    for key in required:
        if key not in stats:
            report.add("error", ds_id, "stats", f"stats.{key} missing")
    if kind == "real" and "cell_vendor" not in df.columns:
        return                            # column issue already reported

    n_mol = len(df)
    if stats.get("n_molecules") != n_mol:
        report.add("error", ds_id, "stats",
                   f"stats.n_molecules {stats.get('n_molecules')} != {n_mol} rows")
    if "gene" in df.columns:
        n_genes = int(df["gene"].nunique())
        if stats.get("n_genes") != n_genes:
            report.add("error", ds_id, "stats",
                       f"stats.n_genes {stats.get('n_genes')} != {n_genes} observed")

    area = stats.get("area_um2")
    if not isinstance(area, (int, float)) or area <= 0:
        report.add("error", ds_id, "stats", f"stats.area_um2 invalid: {area!r}")
        return
    bb = (meta.get("crop") or {}).get("bbox_um")
    if isinstance(bb, list) and len(bb) == 4:
        bbox_area = float((bb[2] - bb[0]) * (bb[3] - bb[1]))
        if not _close(area, bbox_area, abs_tol=1.0, rel_tol=1e-4):
            report.add("error", ds_id, "stats",
                       f"stats.area_um2 {area} != bbox area {bbox_area}")
    if stats.get("molecules_per_um2") is not None and not _close(
            stats["molecules_per_um2"], n_mol / area, abs_tol=1e-6, rel_tol=1e-3):
        report.add("error", ds_id, "stats",
                   f"stats.molecules_per_um2 {stats['molecules_per_um2']} != "
                   f"{n_mol / area:.6g}")

    n_vendor = stats.get("n_vendor_cells")
    if not isinstance(n_vendor, int) or n_vendor < 0:
        report.add("error", ds_id, "stats",
                   f"stats.n_vendor_cells must be a non-negative int, got {n_vendor!r}")
    else:
        if kind == "real" and n_vendor == 0:
            report.add("error", ds_id, "stats", "real dataset with n_vendor_cells == 0")
        vrate = stats.get("vendor_cells_per_mm2")
        if vrate is not None and not _close(vrate, n_vendor / area * 1e6,
                                            abs_tol=1e-3, rel_tol=1e-3):
            report.add("error", ds_id, "stats",
                       f"stats.vendor_cells_per_mm2 {vrate} != "
                       f"{n_vendor / area * 1e6:.6g}")

    if kind == "sim":
        if "cell" not in df.columns:
            return
        n_true = int(df.loc[df["cell"] > 0, "cell"].nunique())
        if stats.get("n_true_cells") != n_true:
            report.add("error", ds_id, "stats",
                       f"stats.n_true_cells {stats.get('n_true_cells')} != {n_true}")
        trate = stats.get("true_cells_per_mm2")
        if trate is not None and not _close(trate, n_true / area * 1e6,
                                            abs_tol=1e-3, rel_tol=1e-3):
            report.add("error", ds_id, "stats",
                       f"stats.true_cells_per_mm2 {trate} != {n_true / area * 1e6:.6g}")


def _validate_difficulty(ds_id: str, kind: Optional[str], meta: dict,
                         report: Report) -> None:
    diff = meta["difficulty"]
    stats = meta["stats"]
    if diff.get("gene_panel") in PANEL_CLASSES and stats.get("n_genes") is not None:
        want = gene_panel_class(int(stats["n_genes"]))
        if diff["gene_panel"] != want:
            report.add("error", ds_id, "difficulty",
                       f"gene_panel {diff['gene_panel']!r} but {stats['n_genes']} "
                       f"genes is class {want!r}")
    density = diff.get("cell_density")
    if density in DENSITY_CLASSES and isinstance(stats.get("area_um2"), (int, float)):
        key = "n_true_cells" if kind == "sim" else "n_vendor_cells"
        n = stats.get(key)
        if isinstance(n, int) and stats["area_um2"] > 0:
            want = cell_density_class(n / stats["area_um2"] * 1e6)
            if density != want:
                report.add("error", ds_id, "difficulty",
                           f"cell_density {density!r} but {key}/area = "
                           f"{n / stats['area_um2'] * 1e6:.0f}/mm2 is class {want!r}")


def _validate_baysor(ds_id: str, ds_dir: Path, repo: Path, kind: Optional[str],
                     meta: dict, schema, report: Report) -> None:
    cfg = meta["baysor"]
    scale = cfg.get("scale_um")
    if not isinstance(scale, (int, float)) or scale <= 0:
        report.add("error", ds_id, "baysor", f"scale_um must be > 0, got {scale!r}")
    if "scale_std" not in cfg:
        report.add("error", ds_id, "baysor", "scale_std missing")
    mmpc = cfg.get("min_molecules_per_cell")
    if not isinstance(mmpc, int) or mmpc < 1:
        report.add("error", ds_id, "baysor",
                   f"min_molecules_per_cell must be an int >= 1, got {mmpc!r}")
    if not isinstance(cfg.get("extra_args", []), list) or not all(
            isinstance(t, str) for t in cfg.get("extra_args", [])):
        report.add("error", ds_id, "baysor", "extra_args must be a list of strings")

    prior = cfg.get("prior", "none")
    prior = "none" if prior in (None, "") else str(prior)
    if prior not in ("none", "column") and not prior.startswith("image:"):
        report.add("error", ds_id, "baysor",
                   f"prior {prior!r} must be 'none', 'column' or 'image:<path>'")
    elif prior.startswith("image:"):
        rel = prior[len("image:"):]
        if not rel or not (ds_dir / rel).is_file():
            report.add("error", ds_id, "files",
                       f"prior image '{rel}' referenced by baysor.prior does not exist")
    conf = cfg.get("prior_confidence")
    if prior != "none":
        if not isinstance(conf, (int, float)) or not (0 <= conf <= 1):
            report.add("error", ds_id, "baysor",
                       f"prior_confidence must be in [0, 1] when a prior is used, "
                       f"got {conf!r}")
        if prior == "column" and "prior" not in schema.names:
            report.add("data", ds_id, "columns",
                       "baysor.prior is 'column' but molecules.parquet has no "
                       "'prior' column")

    config = cfg.get("config")
    if config is not None:
        if not isinstance(config, str):
            report.add("error", ds_id, "baysor", f"config must be a path or null: {config!r}")
        elif not (repo / config).is_file():
            report.add("error", ds_id, "files",
                       f"baysor.config '{config}' does not exist relative to the "
                       "repository root")


def _validate_images(ds_id: str, ds_dir: Path, meta: dict,
                     report: Report) -> None:
    images = meta["images"]
    if not isinstance(images, list):
        report.add("error", ds_id, "images", "images must be a list")
        return
    bb = (meta.get("crop") or {}).get("bbox_um")
    prior_rel = str((meta.get("baysor") or {}).get("prior") or "")
    prior_rel = prior_rel[len("image:"):] if prior_rel.startswith("image:") else None
    for i, img in enumerate(images):
        if not isinstance(img, dict):
            report.add("error", ds_id, "images", f"images[{i}] must be an object")
            continue
        for key in ("name", "file", "pixel_size_um", "origin_um"):
            if key not in img:
                report.add("error", ds_id, "images", f"images[{i}].{key} missing")
        rel = img.get("file")
        if not isinstance(rel, str) or not rel:
            continue
        if rel.startswith("/") or ".." in Path(rel).parts:
            report.add("error", ds_id, "images",
                       f"images[{i}].file {rel!r} escapes the dataset directory")
            continue
        path = ds_dir / rel
        if not path.is_file():
            report.add("error", ds_id, "files",
                       f"images[{i}].file {rel!r} does not exist")
            continue
        px = img.get("pixel_size_um")
        org = img.get("origin_um")
        if not isinstance(px, (int, float)) or px <= 0:
            report.add("error", ds_id, "images",
                       f"images[{i}].pixel_size_um must be > 0, got {px!r}")
            continue
        if not (isinstance(org, list) and len(org) == 2
                and all(isinstance(v, (int, float)) for v in org)):
            report.add("error", ds_id, "images",
                       f"images[{i}].origin_um must be [x0, y0], got {org!r}")
            continue
        try:
            import tifffile
            with tifffile.TiffFile(path) as tf:
                shape = tf.series[0].shape
        except Exception as exc:                                # noqa: BLE001
            report.add("error", ds_id, "images",
                       f"images[{i}].file {rel!r} unreadable as TIFF: {exc}")
            continue
        h, w = shape[-2], shape[-1]
        # the image must cover the crop bbox (within one pixel)
        if isinstance(bb, list) and len(bb) == 4:
            covers = (org[0] <= bb[0] + px and org[1] <= bb[1] + px
                      and org[0] + w * px >= bb[2] - px
                      and org[1] + h * px >= bb[3] - px)
            if not covers:
                report.add("error", ds_id, "images",
                           f"images[{i}] ({rel}) does not cover the crop bbox "
                           f"(origin {org}, {w}x{h} @ {px} um/px vs bbox {bb})")
            # the origin should sit at the bbox corner, unless this is the
            # Baysor image-prior raster (1 um/px aligned to the origin)
            is_prior = prior_rel == rel
            aligned_prior = is_prior and px == 1.0 and org == [0, 0]
            if not aligned_prior:
                if (abs(org[0] - bb[0]) > IMAGE_FRAME_TOL_PX * px
                        or abs(org[1] - bb[1]) > IMAGE_FRAME_TOL_PX * px):
                    report.add("error", ds_id, "images",
                               f"images[{i}] ({rel}) origin {org} is not at the "
                               f"bbox corner {bb[:2]} (tol {IMAGE_FRAME_TOL_PX} px)")


def _validate_manifest(ds_id: str, meta: dict, repo: Path,
                       manifests: dict[str, dict], mol_path: Path,
                       meta_path: Path, report: Report,
                       check_hashes: bool = True) -> None:
    hit = manifests.get(ds_id)
    if hit is None:
        report.add("error", ds_id, "manifest",
                   "dataset is not listed in any of benchmarks/datasets/*.yaml")
        return
    entry = hit["entry"]
    group = hit["group"]
    if entry.get("tier") and entry["tier"] != meta.get("tier"):
        report.add("error", ds_id, "manifest",
                   f"tier {entry['tier']!r} in {group}.yaml != meta tier "
                   f"{meta.get('tier')!r}")
    if entry.get("platform") and entry["platform"] != meta.get("platform"):
        report.add("error", ds_id, "manifest",
                   f"platform {entry['platform']!r} in {group}.yaml != meta "
                   f"platform {meta.get('platform')!r}")

    recorded = dict(entry.get("sha256") or {})
    recorded.update(entry.get("outputs") or {})
    if check_hashes:
        for fname, path in (("molecules.parquet", mol_path), ("meta.json", meta_path)):
            want = recorded.get(fname)
            if want:
                got = _sha256(path)
                if got != want:
                    report.add("error", ds_id, "manifest",
                               f"{fname} sha256 mismatch vs {group}.yaml "
                               f"(recorded {want[:12]}..., on disk {got[:12]}...)")

    # the manifest's baysor block must be mirrored into meta.json
    mb = hit.get("baysor") or {}
    meta_b = meta.get("baysor") or {}
    for key, want in mb.items():
        if key not in meta_b:
            report.add("error", ds_id, "manifest",
                       f"baysor.{key} present in {group}.yaml but missing in meta.json")
        elif meta_b[key] != want:
            report.add("error", ds_id, "manifest",
                       f"baysor.{key}: {group}.yaml has {want!r} but meta.json "
                       f"has {meta_b[key]!r}")


def validate_all(root: Path, repo: Path, manifests: Optional[dict[str, dict]] = None,
                 spec: Optional[str] = None, kind: Optional[str] = None,
                 check_hashes: bool = True) -> Report:
    """Validate every selected dataset; returns the accumulated report."""
    if manifests is None:
        manifests = load_manifest_index(repo)
    report = Report()
    dirs = _select(_discover_dirs(root), spec, root)
    if kind:
        dirs = [d for d in dirs if d.parent.name == kind]
    for d in dirs:
        report.datasets += 1
        validate_dataset(d, repo, manifests, report, check_hashes=check_hashes)
    # manifest entries without a dataset on disk
    on_disk = {d.name for d in _discover_dirs(root)}
    for ds_id in sorted(manifests):
        if ds_id not in on_disk:
            report.add("error", ds_id, "manifest",
                       f"listed in {manifests[ds_id]['group']}.yaml but missing "
                       "under the data root")
    return report


# --- CLI -------------------------------------------------------------------

def format_report(report: Report) -> str:
    lines = []
    by_ds: dict[str, list[Issue]] = {}
    for issue in report.issues:
        by_ds.setdefault(issue.dataset, []).append(issue)
    for ds_id in sorted(by_ds):
        lines.append(f"{ds_id}:")
        for issue in by_ds[ds_id]:
            tag = "ERROR" if issue.severity == "error" else "data "
            lines.append(f"  [{tag}] {issue.check}: {issue.message}")
    lines.append(f"\nvalidated {report.datasets} dataset(s): "
                 f"{report.errors} metadata error(s), "
                 f"{report.data_findings} data finding(s)")
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=None,
                    help="data root (default $BAYSOR_BENCH_DATA or <repo>/.bench-data)")
    ap.add_argument("--repo", default=None, help="repository root (default: auto)")
    ap.add_argument("--datasets", default=None,
                    help="tier (quick|full) or comma-separated ids/globs")
    ap.add_argument("--kind", choices=list(KINDS), default=None)
    ap.add_argument("--json", metavar="PATH", default=None,
                    help="write the machine-readable report as JSON")
    ap.add_argument("--no-hashes", action="store_true",
                    help="skip sha256 verification against the manifests")
    ap.add_argument("--strict", action="store_true",
                    help="treat data findings as failures too")
    ap.add_argument("--quiet", "-q", action="store_true",
                    help="only print the summary line")
    args = ap.parse_args(argv)

    repo = Path(args.repo).resolve() if args.repo else common.repo_root()
    root = common.data_root(args.data_root)
    try:
        report = validate_all(root, repo, spec=args.datasets, kind=args.kind,
                              check_hashes=not args.no_hashes)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if not args.quiet and report.issues:
        print(format_report(report))
    else:
        print(f"validated {report.datasets} dataset(s): "
              f"{report.errors} metadata error(s), "
              f"{report.data_findings} data finding(s)")
    if args.json:
        payload = {
            "data_root": str(root),
            "datasets": report.datasets,
            "errors": report.errors,
            "data_findings": report.data_findings,
            "issues": [asdict(i) for i in report.issues],
        }
        common.write_json(Path(args.json), payload)
    if report.errors:
        return 1
    if args.strict and report.data_findings:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
