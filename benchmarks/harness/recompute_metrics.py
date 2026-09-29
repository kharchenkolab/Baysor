#!/usr/bin/env python3
"""Recompute metrics.json from stored assignment tables (no Baysor rerun).

Re-reads each replicate's ``assignment.parquet`` for every dataset of the
given run(s), recomputes the ``sim`` / ``real`` metric blocks with the
current ``metrics.py`` definitions and rewrites ``metrics.json`` in place.
Runtime, provenance (binary, threads, replicates), per-rep records and the
cellAdmix audit block are preserved untouched.

Dataset content hashes (``inputs.molecules_sha256`` / ``inputs.meta_sha256``)
are recorded when the run predates the runner's own hashing; existing
``inputs`` blocks are kept, because only the runner knows what the binary
actually read. Hashes are computed from the *current* dataset files, so only
use this on runs whose datasets have not been regenerated since.

Typical calibration workflow::

    recompute_metrics.py --run harness-val1 --run rev-same-sim ...
    baseline.py create --run-id harness-val1 --name harness-dev --force
    compare.py --run-id rev-same-sim --baseline harness-dev --expect same
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
import metrics as m                # noqa: E402
import run as runner               # noqa: E402


def recompute_dataset(run_root: Path, metrics_path: Path,
                      dry_run: bool = False) -> dict:
    """Recompute one dataset's metrics.json; returns the new dict."""
    out = common.read_json(metrics_path)
    ds_id = out["dataset"]["id"]
    ds_dir = Path(out["dataset"].get("path") or (run_root.parent.parent / ds_id))
    ds = common.load_dataset(ds_dir)
    if ds is None:
        raise FileNotFoundError(f"dataset not found (needs molecules.parquet "
                                f"+ meta.json): {ds_dir}")

    records = [r for r in out.get("reps", []) if r.get("status") == "ok"]
    if not records:
        print(f"  {ds_id}: no successful replicates, left unchanged")
        return out

    cells = []
    for rec in records:
        rel = rec.get("assignment") or f"rep{rec['rep']}/assignment.parquet"
        path = run_root / ds_id / rel
        cells.append(common.assignment_cells(path))

    if ds.kind == "sim":
        molecules = pd.read_parquet(ds.molecules_path)
        truth = molecules["cell"].to_numpy(np.int64)
        interior = (molecules["interior"].to_numpy(bool)
                    if "interior" in molecules.columns else None)
        oracle = (ds.meta.get("truth") or {}).get("oracle_accuracy")
        per_rep = [m.sim_metrics(c, truth, interior, oracle) for c in cells]
        mean, sd = runner._aggregate_dicts(per_rep)
        out["sim"] = {
            "oracle_accuracy": oracle,
            "per_rep": per_rep,
            "mean": mean,
            "sd": sd,
            "n_metric_reps": len(per_rep),
        }
    else:
        pairs, per_pair = [], []
        ok_reps = [(r["rep"], c) for r, c in zip(records, cells)]
        for i in range(len(ok_reps)):
            for j in range(i + 1, len(ok_reps)):
                pairs.append([ok_reps[i][0], ok_reps[j][0]])
                per_pair.append(m.real_pair_metrics(ok_reps[j][1], ok_reps[i][1]))
        real = out.get("real") or {}
        real["n_metric_reps"] = len(ok_reps)
        if per_pair:
            mean, sd = runner._aggregate_dicts(per_pair)
            real["rep_agreement"] = {"pairs": pairs, "per_pair": per_pair,
                                     "mean": mean, "sd": sd}
        if "cell_vendor" in _parquet_columns(ds.molecules_path):
            vend = runner.vendor_labels(pd.read_parquet(
                ds.molecules_path, columns=["cell_vendor"]))
        else:
            vend = None
        if vend is not None:
            per_rep_v = [m.real_pair_metrics(c, vend) for c in cells]
            vm, vsd = runner._aggregate_dicts(per_rep_v)
            real["vs_vendor"] = {"per_rep": per_rep_v, "mean": vm, "sd": vsd,
                                 "information_only": True}
        out["real"] = real

    if not out.get("inputs"):
        out["inputs"] = {
            "molecules_sha256": common.sha256_file(ds.molecules_path),
            "meta_sha256": common.sha256_file(ds_dir / "meta.json"),
        }
    out["metrics_recomputed"] = {
        "script": "recompute_metrics.py",
        "at": common.utc_now(),
        "metrics_module": "metrics.py (current definitions)",
    }
    if not dry_run:
        common.write_json(metrics_path, out)
    return out


def _parquet_columns(path: Path) -> list[str]:
    import pyarrow.parquet as pq
    return pq.read_schema(path).names


def recompute_run(run_id: str, root: Path, dry_run: bool = False) -> int:
    run_root = root / "runs" / run_id
    if not run_root.is_dir():
        print(f"error: run '{run_id}' not found under {run_root}", file=sys.stderr)
        return 2
    files = sorted(run_root.glob("*/metrics.json"))
    if not files:
        print(f"error: no metrics.json under {run_root}", file=sys.stderr)
        return 2
    print(f"run '{run_id}': recomputing {len(files)} dataset(s)"
          + (" (dry run)" if dry_run else ""))
    for mf in files:
        out = recompute_dataset(run_root, mf, dry_run=dry_run)
        kind = out["dataset"]["kind"]
        if kind == "sim":
            mean = (out.get("sim") or {}).get("mean") or {}
            print(f"  {out['dataset']['id']}: accuracy_1to1="
                  f"{mean.get('accuracy_1to1')} ari_assigned="
                  f"{mean.get('ari_assigned')}")
        else:
            ra = (((out.get("real") or {}).get("rep_agreement")) or {}).get("mean") or {}
            print(f"  {out['dataset']['id']}: ari_assigned="
                  f"{ra.get('ari_assigned')} frac_matched="
                  f"{ra.get('frac_cells_matched')}")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True,
                    help="run id (repeatable)")
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--dry-run", action="store_true",
                    help="recompute but do not write")
    args = ap.parse_args(argv)
    root = common.data_root(args.data_root)
    rc = 0
    for run_id in args.run:
        rc = max(rc, recompute_run(run_id, root, dry_run=args.dry_run))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
