#!/usr/bin/env python3
"""Audit the **vendor** segmentation of a baseline's real datasets.

For every selected dataset of an existing baseline this script

1. builds the vendor assignment table (``cell_vendor`` factorized to int,
   0 = unassigned) next to the baseline's data files,
2. transfers the baseline's saved cell types onto the vendor segmentation
   (``celladmix/transfer.py`` — the same typing every run scores with), and
3. runs ``celladmix/audit.py`` with the baseline's ``fixed_pairs.json``.

The result lands in
``$BAYSOR_BENCH_DATA/baselines/<name>/<id>/vendor_audit.json`` and is read
by ``baseline_summary.py`` so Baysor's and the vendor's
``total_admixture_rate`` are directly comparable (same typing, same pairs).

Usage::

    vendor_audit.py --baseline bugfixes-35e8a7e \
        --datasets 'xenium_*_admix' [--threads 6] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common                      # noqa: E402
import celladmix as camix           # noqa: E402


def vendor_assignment(molecules_path: Path, out_path: Path) -> dict:
    """Factorize ``cell_vendor`` into an assignment table (0 = unassigned)."""
    df = pd.read_parquet(molecules_path, columns=["cell_vendor"])
    raw = df["cell_vendor"].fillna("").astype(str).str.strip()
    uniq = sorted(v for v in raw.unique() if v != "")
    mapping = {v: i + 1 for i, v in enumerate(uniq)}
    cell = np.array([mapping.get(v, 0) for v in raw], dtype=np.int64)
    out = pd.DataFrame({"mol_index": np.arange(len(cell), dtype=np.int64),
                        "cell": cell,
                        "confidence": np.full(len(cell), np.nan)})
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_path, index=False)
    return {"n_cells": len(uniq), "n_assigned": int((cell > 0).sum())}


def audit_dataset(ds_id: str, root: Path, name: str, threads: int,
                  repo: Path, force: bool = False) -> dict:
    base_dir = root / "baselines" / name / ds_id
    out_json = base_dir / "vendor_audit.json"
    if out_json.is_file() and not force:
        return {"dataset": ds_id, "status": "skipped (exists)"}
    for need in ("celltypes.parquet", "fixed_pairs.json"):
        if not (base_dir / need).is_file():
            return {"dataset": ds_id, "status": "failed",
                    "reason": f"baseline lacks {need}"}
    # locate the dataset directory
    ds_dir = next((root / kind / ds_id for kind in ("real", "sim")
                   if (root / kind / ds_id / "molecules.parquet").is_file()),
                  None)
    if ds_dir is None:
        return {"dataset": ds_id, "status": "failed",
                "reason": "dataset not found under the data root"}
    assign = base_dir / "vendor_assignment.parquet"
    stats = vendor_assignment(ds_dir / "molecules.parquet", assign)
    celltypes = base_dir / "vendor_celltypes.parquet"
    tr = camix.run_transfer(ds_dir / "molecules.parquet",
                            base_dir / "rep0" / "assignment.parquet",
                            base_dir / "celltypes.parquet",
                            assign, celltypes,
                            report=base_dir / "vendor_transfer.json",
                            repo=repo)
    if tr.get("status") != camix.STATUS_OK:
        return {"dataset": ds_id, "status": "failed",
                "reason": f"typing transfer failed: {tr.get('reason') or tr.get('stderr_tail', '')[-400:]}"}
    audit = camix.run_audit(ds_dir / "molecules.parquet", assign, out_json,
                            repo=repo, threads=threads,
                            celltypes=celltypes,
                            fixed_pairs=base_dir / "fixed_pairs.json")
    return {"dataset": ds_id, "status": audit.get("status"),
            "total_admixture_rate": audit.get("total_admixture_rate"),
            "n_pairs_evaluated": audit.get("n_pairs_evaluated"),
            "n_cells_vendor": stats["n_cells"],
            "typed_cells": tr.get("n_transferred") or tr.get("n_cells"),
            "audit_json": str(out_json)}


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--datasets", default=None,
                    help="comma-separated ids/globs (default: every real "
                         "dataset of the baseline that has audit typing)")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--repo", default=None)
    ap.add_argument("--force", action="store_true", help="rerun existing audits")
    args = ap.parse_args(argv)

    repo = Path(args.repo).resolve() if args.repo else common.repo_root()
    root = common.data_root(args.data_root)
    base_dir = root / "baselines" / args.baseline
    if not base_dir.is_dir():
        print(f"error: baseline data not found at {base_dir}", file=sys.stderr)
        return 2

    baseline_ids = sorted(d.name for d in base_dir.iterdir() if d.is_dir())
    if args.datasets:
        try:
            selected = [d.id for d in common.select_datasets(root, args.datasets,
                                                             kind="real")]
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        selected = [i for i in selected if i in baseline_ids]
    else:
        selected = [i for i in baseline_ids
                    if (base_dir / i / "celltypes.parquet").is_file()
                    and (base_dir / i / "fixed_pairs.json").is_file()]
    if not selected:
        print("error: nothing selected", file=sys.stderr)
        return 2

    rc = 0
    for ds_id in selected:
        res = audit_dataset(ds_id, root, args.baseline, args.threads, repo,
                            force=args.force)
        print(res, flush=True)
        if res["status"] == "failed":
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
