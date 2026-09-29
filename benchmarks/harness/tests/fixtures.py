"""Contract-conformant fixture datasets and synthetic run helpers for tests.

Everything here follows ``benchmarks/README.md``: ``molecules.parquet`` sorted
by (y, x) plus a ``meta.json`` with tier, baysor parameters and truth.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

BAYSOR_DEFAULTS = {
    "scale_um": 5.0,
    "scale_std": "25%",
    "min_molecules_per_cell": 10,
    "prior": "none",
    "prior_confidence": 0.5,
    "config": None,
    "extra_args": [],
}


def _write(path: Path, df: pd.DataFrame, meta: dict) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path / "molecules.parquet", index=False)
    with open(path / "meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    return path


def _gene_for(cell: int, rng, n_genes: int = 15) -> str:
    """Cell-specific gene pool (disjoint per cell) so segmenters can separate
    cells by expression; background noise draws from all genes."""
    if cell > 0:
        start = ((cell - 1) * 2) % n_genes
        g = int(rng.choice([start, (start + 1) % n_genes]))
    else:
        g = int(rng.integers(0, n_genes))
    return f"gene{g}"


def make_sim_dataset(path: Path, *, ds_id: str = None, n_cells: int = 6,
                     per_cell: int = 40, noise: int = 30, seed: int = 0,
                     tier: str = "quick", baysor_overrides: dict | None = None,
                     oracle: float = 0.99, interior_frac: float = 1.0) -> Path:
    rng = np.random.default_rng(seed)
    rows = []
    truth = {}
    for cell in range(1, n_cells + 1):
        cx = 4.0 + (cell % 3) * 9.0
        cy = 4.0 + (cell // 3) * 9.0
        truth[cell] = (cx, cy)
        n = rng.poisson(per_cell) + 5
        xs = rng.normal(cx, 1.6, n)
        ys = rng.normal(cy, 1.6, n)
        for x, y in zip(xs, ys):
            is_noise = rng.random() < 0.05
            rows.append((x, y, _gene_for(0 if is_noise else cell, rng),
                         0 if is_noise else cell))
    for _ in range(noise):
        rows.append((rng.uniform(0, 30), rng.uniform(0, 30),
                     _gene_for(0, rng), 0))
    df = pd.DataFrame(rows, columns=["x", "y", "gene", "cell"])
    df["interior"] = rng.random(len(df)) >= (1.0 - interior_frac)
    df = df.sort_values(["y", "x"]).reset_index(drop=True)
    baysor = dict(BAYSOR_DEFAULTS)
    baysor.update(baysor_overrides or {})
    meta = {
        "id": ds_id or path.name,
        "kind": "sim",
        "tier": tier,
        "platform": "fixture",
        "source": {"url": "fixture", "retrieved": "2026-01-01"},
        "crop": {"bbox_um": [0, 0, 30, 30], "z_range_um": None, "note": "test"},
        "stats": {"n_molecules": int(len(df)),
                  "n_genes": int(df["gene"].nunique()),
                  "area_um2": 900.0},
        "difficulty": {"cell_density": "medium", "gene_panel": "tiny", "notes": ""},
        "baysor": baysor,
        "images": [],
        "truth": {"generator": "fixtures.make_sim_dataset",
                  "params": {"seed": seed}, "seed": seed,
                  "oracle_accuracy": oracle},
    }
    return _write(path, df, meta)


def make_real_dataset(path: Path, *, ds_id: str = None, n_cells: int = 6,
                      per_cell: int = 40, noise: int = 30, seed: int = 1,
                      tier: str = "quick", baysor_overrides: dict | None = None,
                      with_vendor: bool = True) -> Path:
    rng = np.random.default_rng(seed)
    rows = []
    for cell in range(1, n_cells + 1):
        cx = 4.0 + (cell % 3) * 9.0
        cy = 4.0 + (cell // 3) * 9.0
        n = rng.poisson(per_cell) + 5
        for x, y in zip(rng.normal(cx, 1.6, n), rng.normal(cy, 1.6, n)):
            is_noise = rng.random() < 0.05
            vendor = "" if is_noise else f"V{cell}"
            rows.append((x, y, _gene_for(0 if is_noise else cell, rng),
                         0 if is_noise else cell, vendor))
    for _ in range(noise):
        rows.append((rng.uniform(0, 30), rng.uniform(0, 30),
                     _gene_for(0, rng), 0, ""))
    cols = ["x", "y", "gene", "cell"]
    if with_vendor:
        cols.append("cell_vendor")
    df = pd.DataFrame(rows, columns=cols)
    df = df.sort_values(["y", "x"]).reset_index(drop=True)
    baysor = dict(BAYSOR_DEFAULTS)
    baysor.update(baysor_overrides or {})
    meta = {
        "id": ds_id or path.name,
        "kind": "real",
        "tier": tier,
        "platform": "fixture",
        "source": {"url": "fixture", "retrieved": "2026-01-01"},
        "crop": {"bbox_um": [0, 0, 30, 30], "z_range_um": None, "note": "test"},
        "stats": {"n_molecules": int(len(df)),
                  "n_genes": int(df["gene"].nunique()),
                  "area_um2": 900.0, "n_vendor_cells": n_cells},
        "difficulty": {"cell_density": "medium", "gene_panel": "tiny", "notes": ""},
        "baysor": baysor,
        "images": [],
        "truth": None,
    }
    return _write(path, df, meta)


def make_audit_fixture(path: Path, *, seed: int = 0) -> dict:
    """Small synthetic dataset the real cellAdmix audit can score (~3 s).

    48 cells in four 2x2 type blocks (48 genes, type-specific expression,
    350 molecules/cell + 800 background). Leakage is injected with an
    exposure gradient: cells near a block seam receive source-type marker
    molecules proportional to the number of source-type cells among their
    16 nearest neighbours, so exposed bins carry significantly more source
    markers than the zero-exposure baseline and pairs get *detected*.

    Writes ``molecules.parquet``, ``assignment.parquet`` (mol_index, cell,
    confidence; 0 = unassigned) and ``celltypes.parquet``; returns their paths.
    """
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(seed)
    n_types, cells_per_type = 4, 12
    per_cell, genes_per_type = 350, 12
    n_genes = n_types * genes_per_type
    centers, cell_type = {}, {}
    cid = 0
    for t in range(n_types):
        bx, by = (t % 2) * 40.0, (t // 2) * 40.0
        for i in range(4):
            for j in range(3):
                cid += 1
                centers[cid] = (bx + 5 + i * 9 + rng.uniform(-1, 1),
                                by + 5 + j * 11 + rng.uniform(-1, 1))
                cell_type[cid] = t
    prefer = {t: np.arange(t * genes_per_type, (t + 1) * genes_per_type)
              for t in range(n_types)}
    rows, labels = [], []
    for c, (cx, cy) in centers.items():
        w = np.full(n_genes, 1.0)
        w[prefer[cell_type[c]]] = 8.0
        p = w / w.sum()
        g = rng.choice(n_genes, size=per_cell, p=p)
        rows.extend(zip(cx + rng.normal(0, 2.5, per_cell),
                        cy + rng.normal(0, 2.5, per_cell),
                        [f"g{i}" for i in g]))
        labels.extend([c] * per_cell)
    df = pd.DataFrame(rows, columns=["x", "y", "gene"])
    df["cell"] = labels

    # exposure-correlated leakage across block seams (16-NN, self excluded)
    cids = sorted(centers)
    tree = cKDTree(np.array([centers[c] for c in cids]))
    _, idx = tree.query(np.array([centers[c] for c in cids]), k=16)
    leak_rows, leak_labels = [], []
    for ci, c in enumerate(cids):
        neigh = [cids[j] for j in idx[ci] if cids[j] != c][:15]
        expo: dict[int, int] = {}
        for nn in neigh:
            s = cell_type[nn]
            if s != cell_type[c]:
                expo[s] = expo.get(s, 0) + 1
        cx, cy = centers[c]
        for s, e in expo.items():
            n_leak = int(round(8 * e))
            if n_leak <= 0:
                continue
            w = np.zeros(n_genes)
            w[prefer[s]] = 1.0
            g = rng.choice(n_genes, size=n_leak, p=w / w.sum())
            leak_rows.extend(zip(cx + rng.normal(0, 1.5, n_leak),
                                 cy + rng.normal(0, 1.5, n_leak),
                                 [f"g{i}" for i in g]))
            leak_labels.extend([c] * n_leak)
    leak = pd.DataFrame(leak_rows, columns=["x", "y", "gene"])
    leak["cell"] = leak_labels
    df = pd.concat([df, leak], ignore_index=True)

    nb = 800
    bg = pd.DataFrame({"x": rng.uniform(0, 80, nb), "y": rng.uniform(0, 80, nb),
                       "gene": [f"g{i}" for i in rng.integers(0, n_genes, nb)],
                       "cell": [0] * nb})
    df = pd.concat([df, bg], ignore_index=True)
    df = df.sort_values(["y", "x"]).reset_index(drop=True)

    path.mkdir(parents=True, exist_ok=True)
    df[["x", "y", "gene"]].to_parquet(path / "molecules.parquet", index=False)
    pd.DataFrame({
        "mol_index": np.arange(len(df), dtype=np.int64),
        "cell": df["cell"].astype(np.int64),
        "confidence": np.full(len(df), 0.9),
    }).to_parquet(path / "assignment.parquet", index=False)
    pd.DataFrame({
        "cell": [str(c) for c in sorted(centers)],
        "celltype": [f"type{cell_type[c]}" for c in sorted(centers)],
    }).to_parquet(path / "celltypes.parquet", index=False)
    return {"dir": path,
            "molecules": path / "molecules.parquet",
            "assignment": path / "assignment.parquet",
            "celltypes": path / "celltypes.parquet"}


def corrupt(cells: np.ndarray, fraction: float, seed: int = 7) -> np.ndarray:
    """Randomly reassign ``fraction`` of molecules to wrong cells."""
    rng = np.random.default_rng(seed)
    out = cells.copy()
    n = len(out)
    k = int(round(fraction * n))
    idx = rng.choice(n, size=k, replace=False)
    others = np.unique(out[out > 0])
    if len(others) == 0:
        others = np.array([1])
    for i in idx:
        choices = others[others != out[i]]
        if len(choices) == 0:
            choices = others
        out[i] = int(rng.choice(choices))
    return out


def write_assignment(rep_dir: Path, cells: np.ndarray,
                     confidence: np.ndarray | None = None) -> Path:
    rep_dir.mkdir(parents=True, exist_ok=True)
    if confidence is None:
        confidence = np.full(len(cells), np.nan)
    df = pd.DataFrame({
        "mol_index": np.arange(len(cells), dtype=np.int64),
        "cell": cells.astype(np.int64),
        "confidence": confidence.astype(np.float64),
    })
    out = rep_dir / "assignment.parquet"
    df.to_parquet(out, index=False)
    return out


def make_rep_records(rep_dirs_cells: dict[Path, np.ndarray],
                     celladmix_status: str = "absent") -> tuple[list, dict]:
    """Build (rep_records, rep_dir_for) for compute_dataset_metrics tests."""
    import common
    records, dirs = [], {}
    for i, (rep_dir, cells) in enumerate(
            sorted(rep_dirs_cells.items(),
                   key=lambda kv: int(kv[0].name.replace("rep", "") or 0))):
        assign = rep_dir / "assignment.parquet"
        records.append({
            "rep": i, "status": "ok", "exit_code": 0,
            "wall_s": 1.0 + 0.1 * i, "peak_rss_kb": 100_000,
            "command": ["baysor", "run"], "command_str": "baysor run",
            "assignment": f"{rep_dir.name}/assignment.parquet",
            "assignment_sha256": common.sha256_file(assign),
            "n_cells": int(len(np.unique(cells[cells > 0]))),
            "n_assigned": int((cells > 0).sum()),
            "seg_source": "molecules.parquet",
            "celladmix": {"status": celladmix_status},
        })
        dirs[i] = rep_dir

    def rep_dir_for(rec):
        return dirs[rec["rep"]]

    return records, rep_dir_for


FAKE_BINARY = {
    "path": "/nonexistent/baysor",
    "sha256": "0" * 64,
    "version_info": "Baysor (fixture)",
    "label": "fixturelabel",
}


def make_run(root: Path, run_id: str, ds_dir: Path,
             assignments: list[np.ndarray], *, threads: int = 6,
             celladmix_rates: list | None = None,
             wall_s: float = 1.0, binary: dict | None = None,
             statuses: list[str] | None = None,
             inputs: dict | str | None = "auto") -> dict:
    """Write a synthetic run (assignments + run.json + metrics.json) for tests.

    ``assignments[k]`` is the cell vector for replicate k; it must have one
    entry per input molecule. ``inputs`` controls the dataset content hashes
    recorded under ``metrics["inputs"]``: ``"auto"`` hashes the actual
    ``molecules.parquet``/``meta.json``, a dict is used verbatim and ``None``
    omits the block. Returns the metrics dict.
    """
    import common
    import run as runner
    ds = common.load_dataset(ds_dir)
    assert ds is not None, ds_dir
    records, dirs = {}, {}
    for k, cells in enumerate(assignments):
        rep_dir = root / "runs" / run_id / ds.id / f"rep{k}"
        write_assignment(rep_dir, cells)
        status = (statuses[k] if statuses else "ok")
        rec = {
            "rep": k, "status": status,
            "exit_code": 0 if status == "ok" else 1,
            "wall_s": wall_s + 0.1 * k, "peak_rss_kb": 100_000,
            "command": ["baysor", "run"], "command_str": "baysor run",
            "assignment": f"rep{k}/assignment.parquet",
            "threads": threads,
            "binary_sha256": dict(FAKE_BINARY, **(binary or {}))["sha256"],
            "scale_factor": 1.0,
        }
        if status == "ok":
            rec["assignment_sha256"] = common.sha256_file(
                rep_dir / "assignment.parquet")
            rec["n_cells"] = int(len(np.unique(cells[cells > 0])))
            rec["n_assigned"] = int((cells > 0).sum())
            rec["seg_source"] = "molecules.parquet"
            if celladmix_rates is not None:
                rec["celladmix"] = {"status": "ok",
                                    "total_admixture_rate": celladmix_rates[k]}
            else:
                rec["celladmix"] = {"status": "absent"}
        records[k] = rec
        dirs[k] = rep_dir

    def rep_dir_for(rec):
        return dirs[rec["rep"]]

    mjson = runner.compute_dataset_metrics(
        ds, run_id, rep_dir_for, [records[k] for k in sorted(records)],
        dict(FAKE_BINARY, **(binary or {})), threads)
    if isinstance(inputs, dict):
        mjson["inputs"] = dict(inputs)
    elif inputs == "auto":
        mjson["inputs"] = {
            "molecules_sha256": common.sha256_file(ds_dir / "molecules.parquet"),
            "meta_sha256": common.sha256_file(ds_dir / "meta.json"),
        }
    elif inputs is None:
        # Simulate a pre-provenance run: run.py now records the hashes itself.
        mjson.pop("inputs", None)
    common.write_json(root / "runs" / run_id / ds.id / "metrics.json", mjson)
    return mjson
