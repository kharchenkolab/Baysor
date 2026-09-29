#!/usr/bin/env python
"""cellAdmix admixture audit as a benchmark metric for Baysor segmentations.

Runs the cellAdmix-core admixture audit (``fit.audit_admixture()``) on a
contract-format molecule table plus a segmentation, and writes a JSON result
with the dataset-wide admixture rate, the top per-pair rates, counts, applied
filters, the pinned cellAdmix commit and every parameter.

The audit is an exposure-gradient estimate of leaked molecules per ordered
cell-type pair and needs no reference: for each ordered pair (S -> T) target
cells are stratified by the number of source-type cells among their nearest
neighbours, the excess of source markers in exposed over unexposed target
cells is extrapolated by the markers' transcriptome share, and the total rate
is sum_S->T A_{S->T} / M. Lower admixture means a cleaner segmentation.

Cell typing (needed by the audit) comes either from ``--celltypes`` or from
cellAdmix's quick clustering run with a fixed seed
(``--cluster-resolution``). Typing depends on the segmentation; see README.md
for the comparability discussion (transfer baseline labels vs re-cluster).

Unassigned molecules are excluded, matching ``keep_unassigned = FALSE``.

Example:
    python audit.py --molecules molecules.parquet --cell-column cell_vendor \
        --out audit_vendor.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# Pinned cellAdmix-core commit this benchmark runs against (see install.sh).
CELLADMIX_COMMIT = "7d3fe7ae70c61d2b9e57469d38d9a88fcf6ac14d"
CELLADMIX_REPO = "https://github.com/kharchenkolab/cellAdmix-core"

# Values the tabular loader treats as "no cell" (src/input_store.cpp).
UNASSIGNED_TOKENS = {"", "0", "NA", "NaN", "null", "UNASSIGNED", "unassigned", "cell_0"}

# Audit parameters. n_pool deviates from cellAdmix's default of 20: on small
# panels (e.g. 376-gene Xenium) many genes crowd at contrast ~ 2, so the
# 20-gene cutoff reshuffles under small perturbations and the pool's
# transcriptome coverage (the extrapolation factor) swings wildly between
# segmentations (observed 0.018 vs 0.500 for one pair), destroying
# comparability. 60 stabilizes it; see README.md.
AUDIT_PARAM_DEFAULTS = {
    "neighbor_k": 15,
    "n_pool": 60,
    "q_thresh": 0.01,
    "min_excess": 200,
    "min_target_cells": 200,
    "min_reference_cells": 100,
}


def sha256_file(path: Path, chunk: int = 1 << 22) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def load_molecules(molecules: Path) -> pd.DataFrame:
    """Read x/y/gene plus any requested columns from the molecule table."""
    return pd.read_parquet(molecules)


def resolve_assignment(args, molecules: pd.DataFrame) -> tuple[np.ndarray, dict]:
    """Per-molecule cell labels (str, '' = unassigned) + provenance info."""
    info: dict = {}
    if args.cell_column:
        labels = molecules[args.cell_column].fillna("").astype(str).to_numpy()
        info["assignment"] = {"kind": "cell_column", "column": args.cell_column}
    else:
        assign = pd.read_parquet(args.assignment)
        n_mol = len(molecules)
        if len(assign) != n_mol:
            raise SystemExit(
                f"--assignment has {len(assign)} rows but --molecules has {n_mol}; "
                "the assignment must be in the same row order")
        col = args.assignment_col
        if col not in assign.columns:
            raise SystemExit(f"--assignment-col {col!r} not in {sorted(assign.columns)}")
        codes = assign[col].to_numpy()
        if "cell_label" in assign.columns:
            labels = assign["cell_label"].fillna("").astype(str).to_numpy()
            label_col = "cell_label"
        else:
            labels = np.where(pd.notna(codes) & (codes != 0), codes.astype(str), "")
            label_col = None
        info["assignment"] = {
            "kind": "assignment_file",
            "path": str(args.assignment),
            "column": col,
            "label_column": label_col,
            "sha256": sha256_file(args.assignment),
        }
    labels = labels.astype(object)
    labels[np.asarray([str(v) in UNASSIGNED_TOKENS for v in labels])] = ""
    return labels, info


def filter_rows(molecules: pd.DataFrame, labels: np.ndarray) -> tuple[pd.DataFrame, np.ndarray, dict]:
    """Drop rows violating the contract (finite x/y, non-null gene) in sync."""
    keep = (
        np.isfinite(molecules["x"].to_numpy(dtype=float))
        & np.isfinite(molecules["y"].to_numpy(dtype=float))
        & molecules["gene"].notna().to_numpy()
        & (molecules["gene"].astype(str).to_numpy() != "")
    )
    dropped = int((~keep).sum())
    if dropped:
        molecules = molecules[keep].reset_index(drop=True)
        labels = labels[keep]
    return molecules, labels, {"malformed_rows_dropped": dropped}


def write_audit_molecules(src: pd.DataFrame, labels: np.ndarray, out_path: Path) -> int:
    """Write the x/y/gene/cell table consumed by the tabular store builder."""
    assigned = labels != ""
    table = pd.DataFrame({
        "x": src["x"].to_numpy(dtype="float64")[assigned],
        "y": src["y"].to_numpy(dtype="float64")[assigned],
        "gene": src["gene"].astype(str).to_numpy()[assigned],
        "cell": labels[assigned].astype(str),
    })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(out_path, index=False)
    return len(table)


def top_pairs(pairs: pd.DataFrame, n: int) -> list[dict]:
    """Top-n detected pairs by estimated admixed molecules (ties by rate)."""
    if pairs.empty or "detected" not in pairs.columns:
        return []
    det = pairs[pairs["detected"]].sort_values(
        ["admixed_molecules", "rate"], ascending=False, kind="mergesort").head(n)
    records = []
    for _, row in det.iterrows():
        records.append({
            "source": str(row["source"]),
            "target": str(row["target"]),
            "rate": float(row["rate"]),
            "admixed_molecules": float(row["admixed_molecules"]),
            "excess": float(row["excess"]),
            "coverage": float(row["coverage"]),
            "q_value": float(row["q_value"]),
            "n_exposed": int(row["n_exposed"]),
            "n_reference": int(row["n_reference"]),
            "n_markers": int(row["n_markers"]),
            "n_strict": int(row["n_strict"]),
        })
    return records


def compute_metrics(pairs: pd.DataFrame, total_molecules: float) -> dict:
    """Dataset-wide totals from the per-pair audit table."""
    if pairs.empty or "detected" not in pairs.columns:
        return {
            "total_admixture_rate": 0.0,
            "total_admixture_molecules": 0.0,
            "n_pairs_evaluated": 0,
            "n_pairs_detected": 0,
            "status": "no_pairs_evaluated",
        }
    det = pairs[pairs["detected"]]
    total = float(det["admixed_molecules"].sum()) if len(det) else 0.0
    denom = max(float(total_molecules), 1.0)
    return {
        "total_admixture_rate": total / denom,
        "total_admixture_molecules": total,
        "n_pairs_evaluated": int(len(pairs)),
        "n_pairs_detected": int(len(det)),
        "status": "ok" if len(det) else "no_detected_pairs",
    }


def run_audit(args) -> dict:
    """Execute the full audit pipeline and return the result dict."""
    import store as store_mod
    from transfer import read_celltypes

    timings: dict[str, float] = {}
    t_total = time.perf_counter()

    molecules_path = Path(args.molecules)
    required = {"x", "y", "gene"}
    import pyarrow.parquet as pq

    have = set(pq.read_schema(molecules_path).names)
    if not required <= have:
        raise SystemExit(f"molecules table lacks required columns: {sorted(required - have)}")
    if bool(args.cell_column) == bool(args.assignment):
        raise SystemExit("exactly one of --cell-column / --assignment is required")

    t0 = time.perf_counter()
    molecules = pd.read_parquet(
        molecules_path,
        columns=["x", "y", "gene"] + ([args.cell_column] if args.cell_column else []),
    )
    n_input = len(molecules)
    labels, assignment_info = resolve_assignment(args, molecules)
    molecules, labels, row_filters = filter_rows(molecules, labels)
    timings["load"] = time.perf_counter() - t0

    n_unassigned = int(np.sum(labels == ""))
    filters = {
        "unassigned_molecules_dropped": n_unassigned,  # keep_unassigned = FALSE
        **row_filters,
        "untyped_cell_molecules_dropped": 0,
    }

    # Optional pre-given typing: restrict the universe to typed cells now so
    # the store (and later the run) never contains an untyped cell.
    annotation = None
    typing_mode = "quick_cluster"
    if args.celltypes:
        typing_mode = "celltypes"
        types = read_celltypes(args.celltypes,
                               cell_col=args.celltype_cell_col,
                               type_col=args.celltype_type_col)
        typed_set = set(types.index)
        typed = np.array([str(v) in typed_set for v in labels])
        keep = (labels == "") | typed  # unassigned rows go away in the next step anyway
        filters["untyped_cell_molecules_dropped"] = int((~keep).sum())
        molecules = molecules[keep].reset_index(drop=True)
        labels = labels[keep]
        annotation = types

    work_dir = Path(args.work_dir) if args.work_dir else Path(args.out).with_suffix("").parent / (Path(args.out).stem + "_work")
    work_dir.mkdir(parents=True, exist_ok=True)
    temp_path = work_dir / "molecules_with_cell.parquet"
    n_written = write_audit_molecules(molecules, labels, temp_path)

    # --- input store -----------------------------------------------------
    t0 = time.perf_counter()
    ds = store_mod.build_tabular_dataset(
        temp_path, work_dir, annotation=annotation, num_threads=args.threads)
    ds.ensure_store(force=args.overwrite)
    timings["store"] = time.perf_counter() - t0
    if annotation is not None:
        annotation = ds.annotation  # normalized by read_annotation

    # --- typing ----------------------------------------------------------
    t0 = time.perf_counter()
    if annotation is None:
        for _ in range(3):
            types = store_mod.quick_cluster(
                ds.store_dir,
                resolution=args.cluster_resolution,
                seed=args.seed,
                min_molecules=args.cluster_min_molecules,
            )
            store_cells = set(store_mod.read_store_cells(ds.store_dir))
            missing = store_cells - set(types.index)
            if not missing:
                break
            # Cells the clustering gate rejected would become a pseudo-type
            # in the audit; drop their molecules and rebuild once.
            typed_set = set(types.index)
            keep = np.array([str(v) in typed_set for v in labels])
            filters["untyped_cell_molecules_dropped"] += int((~keep).sum())
            molecules = molecules[keep].reset_index(drop=True)
            labels = labels[keep]
            n_written = write_audit_molecules(molecules, labels, temp_path)
            ds.ensure_store(force=True)
        else:
            raise SystemExit("clustering left cells unlabeled after repeated filtering")
        annotation = types
        ds.set_annotation(types)
        annotation = ds.annotation
    timings["typing"] = time.perf_counter() - t0

    run_cells = set(store_mod.read_store_cells(ds.store_dir))
    untyped_run_cells = run_cells - set(annotation.index)
    if untyped_run_cells:
        raise SystemExit(
            f"{len(untyped_run_cells)} store cells lack a cell type; "
            "the audit would treat them as a pseudo-type")

    # --- fit -------------------------------------------------------------
    t0 = time.perf_counter()
    fit = ds.fit(
        nmf_n_runs=args.nmf_n_runs,
        seed=args.seed,
        num_threads=args.threads,
        run_id=f"fit_seed{seed_id(args.seed)}",
        verbose=False,
    )
    timings["fit"] = time.perf_counter() - t0
    fit_cells = set(str(c) for c in fit.cell_factors()["cell_id"])
    missing_fit = fit_cells - set(annotation.index)
    if missing_fit:
        raise SystemExit(f"{len(missing_fit)} fitted cells lack a cell type")

    # --- audit -----------------------------------------------------------
    t0 = time.perf_counter()
    audit = fit.audit_admixture(
        neighbor_k=args.neighbor_k,
        n_pool=args.n_pool,
        q_thresh=args.q_thresh,
        min_excess=args.min_excess,
        min_target_cells=args.min_target_cells,
        min_reference_cells=args.min_reference_cells,
    )
    pairs = audit.pairs()
    matrix, genes, cells = fit.counts()
    timings["audit"] = time.perf_counter() - t0

    total_molecules = float(matrix.sum())
    metrics = compute_metrics(pairs, total_molecules)
    top = top_pairs(pairs, args.top_pairs)

    used_types = annotation.reindex([str(c) for c in cells])
    n_types = int(used_types.dropna().astype(str).nunique())

    if args.save_celltypes:
        out = Path(args.save_celltypes)
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({
            "cell": annotation.index.astype(str),
            "celltype": annotation.astype(str).to_numpy(),
        }).to_parquet(out, index=False)

    timings["total"] = time.perf_counter() - t_total

    import celladmix as ca
    from celladmix import _core

    commit_verified = None
    src_clone = Path(args.celladmix_source) if args.celladmix_source else None
    if src_clone and (src_clone / ".git").exists():
        import subprocess

        head = subprocess.run(
            ["git", "-C", str(src_clone), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip()
        commit_verified = head == CELLADMIX_COMMIT

    fit_params = fit.manifest.get("pipeline_options", {}) or {}
    result = {
        "task": "celladmix-admixture-audit",
        "schema_version": 1,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "inputs": {
            "molecules": str(molecules_path),
            "molecules_sha256": sha256_file(molecules_path),
            "molecules_bytes": molecules_path.stat().st_size,
            "celltypes": str(args.celltypes) if args.celltypes else None,
            "celltypes_sha256": sha256_file(args.celltypes) if args.celltypes else None,
            "image": str(args.image) if args.image else None,
            "image_note": "recorded for membrane scoring; not used by the audit itself",
            **assignment_info,
        },
        "celladmix": {
            "repo": CELLADMIX_REPO,
            "commit": CELLADMIX_COMMIT,
            "commit_verified_against_source_clone": commit_verified,
            "package_version": ca.__version__,
            "core_version": _core.core_version(),
        },
        "parameters": {
            "seed": args.seed,
            "threads": args.threads,
            "typing": typing_mode,
            "cluster_resolution": args.cluster_resolution if typing_mode == "quick_cluster" else None,
            "cluster_min_molecules": args.cluster_min_molecules if typing_mode == "quick_cluster" else None,
            "nmf_n_runs": args.nmf_n_runs,
            "fit": {k: fit_params.get(k) for k in
                    ("rank", "nmf_variant", "nmf_init", "seed", "ncv_k", "graph_k") if k in fit_params},
            "audit": {
                "neighbor_k": args.neighbor_k,
                "n_pool": args.n_pool,
                "q_thresh": args.q_thresh,
                "min_excess": args.min_excess,
                "min_target_cells": args.min_target_cells,
                "min_reference_cells": args.min_reference_cells,
            },
            "top_pairs": args.top_pairs,
        },
        "filters": filters,
        "counts": {
            "n_molecules_input": n_input,
            "n_molecules_written": n_written,
            "n_molecules_used": int(total_molecules),
            "n_cells": int(len(cells)),
            "n_genes": int(len(genes)),
            "n_cell_types": n_types,
        },
        "metrics": metrics,
        "pairs_top": top,
        "runtime_seconds": {k: round(v, 3) for k, v in timings.items()},
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
    }
    return result


def seed_id(seed) -> str:
    """Filesystem/run-id safe seed token."""
    return str(seed).replace("/", "_")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--molecules", type=Path, required=True,
                        help="contract-format molecules parquet (x, y, gene, ...)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--assignment", type=Path,
                       help="per-molecule assignment parquet ('cell' int, 0 = unassigned, "
                            "same row order as --molecules; optional 'cell_label' str column)")
    group.add_argument("--cell-column", help="segmentation column inside --molecules "
                                             "(e.g. cell_vendor or cell)")
    parser.add_argument("--assignment-col", default="cell",
                        help="int code column inside --assignment (default: cell)")
    parser.add_argument("--celltypes", type=Path,
                        help="cell -> type table for fixed typing (else quick clustering)")
    parser.add_argument("--celltype-cell-col")
    parser.add_argument("--celltype-type-col")
    parser.add_argument("--image", type=Path,
                        help="optional membrane-stain image (recorded, not used by the audit)")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument("--cluster-resolution", type=float, default=1.0,
                        help="quick-clustering resolution; higher = more cell types (default 1.0)")
    parser.add_argument("--cluster-min-molecules", type=int, default=1,
                        help="min molecules for a cell to be clustered (default 1: full coverage)")
    parser.add_argument("--nmf-n-runs", type=int, default=1,
                        help="NMF restarts; the audit is factorization-independent, 1 is enough")
    parser.add_argument("--neighbor-k", type=int, default=AUDIT_PARAM_DEFAULTS["neighbor_k"])
    parser.add_argument("--n-pool", type=int, default=AUDIT_PARAM_DEFAULTS["n_pool"])
    parser.add_argument("--q-thresh", type=float, default=AUDIT_PARAM_DEFAULTS["q_thresh"])
    parser.add_argument("--min-excess", type=int, default=AUDIT_PARAM_DEFAULTS["min_excess"])
    parser.add_argument("--min-target-cells", type=int,
                        default=AUDIT_PARAM_DEFAULTS["min_target_cells"])
    parser.add_argument("--min-reference-cells", type=int,
                        default=AUDIT_PARAM_DEFAULTS["min_reference_cells"])
    parser.add_argument("--top-pairs", type=int, default=20)
    parser.add_argument("--out", type=Path, required=True, help="output JSON path")
    parser.add_argument("--work-dir", type=Path,
                        help="scratch dir for the temp table/store/run "
                             "(default: <out stem>_work next to --out)")
    parser.add_argument("--save-celltypes", type=Path,
                        help="write the typing actually used (cell, celltype) parquet")
    parser.add_argument("--overwrite", action="store_true",
                        help="force rebuilding the input store")
    parser.add_argument("--celladmix-source", type=Path,
                        help="cellAdmix-core checkout to verify the pinned commit against")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.image and not Path(args.image).exists():
        raise SystemExit(f"--image not found: {args.image}")
    result = run_audit(args)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
        fh.write("\n")
    summary = {
        "out": str(out),
        "status": result["metrics"]["status"],
        "total_admixture_rate": result["metrics"]["total_admixture_rate"],
        "total_admixture_molecules": result["metrics"]["total_admixture_molecules"],
        "n_pairs_detected": result["metrics"]["n_pairs_detected"],
        "n_cells": result["counts"]["n_cells"],
        "n_molecules_used": result["counts"]["n_molecules_used"],
        "n_cell_types": result["counts"]["n_cell_types"],
        "runtime_seconds": result["runtime_seconds"]["total"],
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
