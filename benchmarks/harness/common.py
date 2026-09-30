"""Shared helpers for the Baysor benchmark harness.

Dataset discovery under ``$BAYSOR_BENCH_DATA``, JSON/parquet I/O, sha256
hashing, and robust parsing/normalisation of Baysor segmentation output into
the harness's per-molecule assignment table.
"""
from __future__ import annotations

import fnmatch
import hashlib
import json
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np
import pandas as pd

KIND_DIRS = {"sim": "sim", "real": "real"}
TIERS = ("quick", "full", "all")


def repo_root() -> Path:
    """Repository root (parent of ``benchmarks/``)."""
    return Path(__file__).resolve().parents[2]


def data_root(cli_value: Optional[str] = None) -> Path:
    """Resolve the shared data root: ``--data-root`` > ``$BAYSOR_BENCH_DATA``
    > ``<repo>/.bench-data``."""
    if cli_value:
        return Path(cli_value).expanduser().resolve()
    import os
    env = os.environ.get("BAYSOR_BENCH_DATA")
    if env:
        return Path(env).expanduser().resolve()
    return repo_root() / ".bench-data"


def baselines_root(cli_value: Optional[str] = None) -> Path:
    """Baseline store ``<data-root>/baselines``: per-dataset metric JSONs,
    ``resources.csv``, ``SUMMARY.md`` and the per-replicate assignment
    tables/cell types (local, never committed)."""
    return data_root(cli_value) / "baselines"


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=False, allow_nan=True)
        f.write("\n")


def read_json(path: Path) -> Any:
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Dataset discovery / selection
# ---------------------------------------------------------------------------

@dataclass
class Dataset:
    id: str
    kind: str                 # "sim" | "real"
    tier: Optional[str]       # "quick" | "full" | None
    path: Path
    meta: dict = field(repr=False, default_factory=dict)

    @property
    def molecules_path(self) -> Path:
        return self.path / "molecules.parquet"

    @property
    def baysor_cfg(self) -> dict:
        return self.meta.get("baysor") or {}

    @property
    def has_z(self) -> bool:
        try:
            import pyarrow.parquet as pq
            schema = pq.read_schema(self.molecules_path)
            return "z" in schema.names
        except Exception:
            return False


def load_dataset(ds_dir: Path) -> Optional[Dataset]:
    meta_path = ds_dir / "meta.json"
    if not meta_path.is_file() or not (ds_dir / "molecules.parquet").is_file():
        return None
    meta = read_json(meta_path)
    kind = meta.get("kind")
    if kind not in KIND_DIRS:
        return None
    tier = meta.get("tier")
    return Dataset(id=meta.get("id", ds_dir.name), kind=kind, tier=tier,
                   path=ds_dir, meta=meta)


def discover_datasets(root: Path, kind: Optional[str] = None) -> list[Dataset]:
    out: list[Dataset] = []
    kinds = [kind] if kind else list(KIND_DIRS)
    for k in kinds:
        base = root / KIND_DIRS[k]
        if not base.is_dir():
            continue
        for d in sorted(base.iterdir()):
            if d.is_dir():
                ds = load_dataset(d)
                if ds is not None:
                    out.append(ds)
    return out


def select_datasets(root: Path, spec: str, kind: Optional[str] = None) -> list[Dataset]:
    """Select datasets from a spec.

    The spec is a comma-separated list of dataset ids and/or shell globs, or
    exactly one manifest tier name (``quick``, ``full`` or ``all``).
    Tier membership comes from ``meta.json`` (the manifest mirrors it).
    """
    all_ds = discover_datasets(root, kind=kind)
    spec = (spec or "").strip()
    if not spec:
        raise ValueError("--datasets is required (ids, globs, or quick|full|all)")
    if spec in TIERS:
        if spec == "all":
            return all_ds
        return [d for d in all_ds if d.tier == spec]
    tokens = [t.strip() for t in spec.split(",") if t.strip()]
    by_id = {d.id: d for d in all_ds}
    out: list[Dataset] = []
    seen: set[str] = set()
    for tok in tokens:
        if tok in by_id:
            matched = [by_id[tok]]
        else:
            matched = [d for d in all_ds if fnmatch.fnmatch(d.id, tok)]
        if not matched:
            raise ValueError(f"no dataset matches '{tok}' under {root}")
        for d in matched:
            if d.id not in seen:
                seen.add(d.id)
                out.append(d)
    return out


# ---------------------------------------------------------------------------
# Baysor segmentation output -> normalized assignment table
# ---------------------------------------------------------------------------

@dataclass
class Segmentation:
    """Parsed Baysor segmentation table in input-molecule order."""
    cell: np.ndarray          # int64, 0 = unassigned/noise
    confidence: np.ndarray    # float64, NaN if the binary provided none
    n_cells: int
    source: str               # "molecules.parquet" | "segmentation.csv"


_CELL_PREFIX = re.compile(r"^cell_")


def cells_to_int(values: Iterable[Any]) -> np.ndarray:
    """Map Baysor cell labels (``0``, ``cell_12``, plain ints) to int64 ids.

    0 always means unassigned/noise. Non-numeric labels are factorized
    deterministically (sorted order) to consecutive ids starting at 1.
    """
    if isinstance(values, (np.ndarray, pd.Series)):
        raw = values.tolist()
    else:
        raw = list(values)
    uniq: dict[str, int] = {}
    out = np.zeros(len(raw), dtype=np.int64)
    next_id = 1
    for i, v in enumerate(raw):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            out[i] = 0
            continue
        if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
            out[i] = int(v)
            continue
        s = str(v).strip()
        if s == "" or s == "0":
            out[i] = 0
            continue
        s2 = _CELL_PREFIX.sub("", s)
        try:
            out[i] = int(s2)
            continue
        except ValueError:
            pass
        if s not in uniq:
            uniq[s] = next_id
            next_id += 1
        out[i] = uniq[s]
    return out


def _find_segmentation_file(seg_dir: Path) -> Path:
    for name in ("molecules.parquet", "segmentation.csv"):
        p = seg_dir / name
        if p.is_file():
            return p
    raise FileNotFoundError(
        f"no Baysor segmentation output (molecules.parquet / segmentation.csv) in {seg_dir}")


def read_segmentation(seg_dir: Path) -> Segmentation:
    """Read Baysor's per-molecule segmentation output (parquet or legacy CSV)."""
    path = _find_segmentation_file(seg_dir)
    if path.suffix == ".parquet":
        df = pd.read_parquet(path, columns=["cell", "assignment_confidence"])
        source = "molecules.parquet"
    else:
        df = pd.read_csv(path, usecols=["cell", "assignment_confidence"])
        source = "segmentation.csv"
    cell = cells_to_int(df["cell"].to_numpy())
    if "assignment_confidence" in df.columns:
        conf = pd.to_numeric(df["assignment_confidence"], errors="coerce").to_numpy(float)
    else:
        conf = np.full(len(df), np.nan)
    return Segmentation(cell=cell, confidence=conf,
                        n_cells=int(len(set(cell[cell > 0]))), source=source)


def _seg_coords_genes(seg_dir: Path) -> pd.DataFrame:
    path = _find_segmentation_file(seg_dir)
    cols = ["gene", "x", "y"]
    if path.suffix == ".parquet":
        return pd.read_parquet(path, columns=cols)
    return pd.read_csv(path, usecols=cols)


def _align_positions(input_df: pd.DataFrame, seg: pd.DataFrame) -> np.ndarray:
    """Return for every input row the index of its row in the segmentation table
    (-1 when the molecule is absent from the output).

    Baysor writes output rows in input order; that is the fast path and is
    verified against coordinates and gene names. The loader may legitimately
    drop molecules (``min_molecules_per_gene``, ``exclude_genes``, ``min_qv``,
    coordinate bounds), in which case lengths differ and we fall back to a
    join on (gene, rounded x, rounded y) with duplicate handling; unmatched
    input rows are reported as -1. Segmentation rows that cannot be matched
    back to any input row are an error (that would mean transformed
    coordinates or a mismatched input file).
    """
    n = len(input_df)
    og = input_df["gene"].astype(str).to_numpy()
    sg = seg["gene"].astype(str).to_numpy()
    ox = input_df["x"].to_numpy(float)
    oy = input_df["y"].to_numpy(float)
    sx = seg["x"].to_numpy(float)
    sy = seg["y"].to_numpy(float)
    if len(seg) == n:
        coords_ok = (np.allclose(ox, sx, rtol=1e-5, atol=1e-4)
                     and np.allclose(oy, sy, rtol=1e-5, atol=1e-4))
        if coords_ok and np.array_equal(og, sg):
            return np.arange(n, dtype=np.int64)

    # Key on gene + coordinates rounded to 1e-3 um (legacy CSV float
    # formatting keeps ~6 significant digits, so exact equality is unsafe).
    def keys(g, x, y):
        return list(zip(g, np.round(x, 3), np.round(y, 3)))

    buckets: dict[tuple, list[int]] = {}
    for i, k in enumerate(keys(sg, sx, sy)):
        buckets.setdefault(k, []).append(i)
    out = np.full(n, -1, dtype=np.int64)
    for i, k in enumerate(keys(og, ox, oy)):
        q = buckets.get(k)
        if q:
            out[i] = q.pop(0)
    leftover = sum(len(v) for v in buckets.values())
    if leftover:
        raise RuntimeError(
            f"{leftover} segmentation rows have no matching input molecule; "
            f"coordinates/genes do not correspond to the input file")
    return out


def normalize_assignment(dataset_dir: Path, seg_dir: Path) -> pd.DataFrame:
    """Build the normalized assignment table for one replicate.

    Returns a DataFrame with columns ``mol_index`` (input row order),
    ``cell`` (int64, 0 = unassigned/noise or filtered by the Baysor loader)
    and ``confidence`` (float64, NaN when the binary does not report an
    assignment for the molecule, e.g. it was dropped during input loading).
    """
    input_df = pd.read_parquet(dataset_dir / "molecules.parquet", columns=["x", "y", "gene"])
    seg = read_segmentation(seg_dir)
    coords = _seg_coords_genes(seg_dir)
    order = _align_positions(input_df, coords)
    n = len(input_df)
    cell = np.zeros(n, dtype=np.int64)
    conf = np.full(n, np.nan)
    matched = order >= 0
    cell[matched] = seg.cell[order[matched]]
    conf[matched] = seg.confidence[order[matched]]
    return pd.DataFrame({
        "mol_index": np.arange(n, dtype=np.int64),
        "cell": cell,
        "confidence": conf,
    })


def load_assignment(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    for col in ("mol_index", "cell"):
        if col not in df.columns:
            raise ValueError(f"assignment table {path} lacks column {col}")
    return df


def assignment_cells(path: Path) -> np.ndarray:
    """Cell vector (int64) of an assignment table, ordered by mol_index."""
    df = load_assignment(path)
    df = df.sort_values("mol_index")
    return df["cell"].to_numpy(np.int64)
