"""Sanity check: run the release Baysor binary on simulated datasets.

Default: one trivial, one sparse and one dense st-recoverability dataset,
runs ``baysor run`` with the ``meta.baysor`` parameters and the ``prior``
column as ``:prior``, and records:

  * wall time and exit code (a run must finish and produce output);
  * a quick molecule-assignment accuracy: predicted cells are matched to true
    cells by majority (per predicted label, its most frequent true cell;
    label 0 / noise never matches), evaluated on interior molecules — a
    plausibility check next to the oracle accuracy from ``meta.truth``;
  * a one-to-one (Hungarian) matched accuracy: the maximum-weight injective
    matching of predicted to true labels over the same molecules.  Unlike the
    majority match, one predicted cell cannot get credit for two true cells
    (a pure split is penalised).  This helper serves the sanity check only;
    the benchmark metric lives in ``harness/metrics.py``.

``--ids ..._noprior --report ...`` runs the no-prior variants (the ``:prior``
argument and prior confidence are dropped, everything else unchanged).

Full metrics are BENCH-HARNESS's job; this only answers "does it run and is
the answer in the plausible range?".

Outputs go to ``$BAYSOR_BENCH_DATA/runs/bench-sim-sanity/<id>/``; a small
report is written to ``$BAYSOR_BENCH_DATA/results/simulate/sanity_check.json``
(local, not committed).

Usage::

    python benchmarks/simulate/sanity.py            # the three default datasets
    python benchmarks/simulate/sanity.py --ids sim_tiled_same_g100
    python benchmarks/simulate/sanity.py \\
        --ids strec_sparse_s1_disjoint_noprior sim_circles_gaps_g100_noprior \\
              strec_dense_s2_merfish_noprior \\
        --report $BAYSOR_BENCH_DATA/results/simulate/sanity_check_noprior.json
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import common  # noqa: E402

DEFAULT_BINARY = "/home/vpetukhov/Projects/Baysor/.bench-data/binaries/baysor-bugfixes-35e8a7e"
DEFAULT_IDS = [
    "sim_circles_gaps_g100",           # trivial
    "strec_sparse_s2_disjoint",        # sparse st-recoverability
    "strec_dense_s2_merfish",          # dense + sigma=2 + realistic
]
REPORT = common.data_root() / "results" / "simulate" / "sanity_check.json"
MAX_THREADS = 6


def majority_match_accuracy(pred: np.ndarray, true: np.ndarray,
                            mask: np.ndarray) -> dict:
    """Match predicted to true cells by majority and score correctness.

    For every predicted label, its most frequent true cell (ties: smallest id)
    becomes the match; predicted label 0 / noise never matches.  Evaluated on
    ``mask`` molecules; returns accuracy plus match statistics.
    """
    p = np.asarray(pred)[mask]
    t = np.asarray(true)[mask]
    if len(p) == 0:
        raise ValueError("empty mask")
    pred_labels = np.unique(p)
    match: dict[int, int] = {}
    correct = 0
    n_multi = 0
    for lab in pred_labels:
        sel = p == lab
        vals, counts = np.unique(t[sel], return_counts=True)
        best = int(vals[np.argmax(counts)])
        match[int(lab)] = best
        if vals.size > 1 and lab > 0:
            n_multi += 1
        if lab > 0 and best > 0:
            correct += int(np.sum(sel & (t == best)))
    return {
        "accuracy": float(correct / len(p)),
        "n_evaluated": int(len(p)),
        "n_predicted_labels": int(pred_labels.size),
        "n_labels_with_mixed_truth": int(n_multi),
        "n_predicted_zero": int(np.sum(p == 0)),
    }


def hungarian_match_accuracy(pred: np.ndarray, true: np.ndarray,
                             mask: np.ndarray) -> dict:
    """One-to-one (Hungarian) matched accuracy of a predicted segmentation.

    Builds the ``(pred label, true label)`` co-occurrence table over ``mask``
    molecules and solves the maximum-weight injective assignment with
    ``scipy.optimize.linear_sum_assignment``.  Molecules inside a matched
    (predicted, true) label pair count as correct; unmapped labels contribute
    0.  A noise prediction (label 0) is only ever matched to true background,
    and a cell prediction never receives credit on background molecules.  A
    pure split — one true cell covered by two predicted labels — can score at
    most the mass of the single best label.  Sanity helper only; the
    benchmark's primary one-to-one metric lives in ``harness/metrics.py``.

    Returns accuracy plus assignment statistics.
    """
    from scipy.optimize import linear_sum_assignment

    p = np.asarray(pred)[mask]
    t = np.asarray(true)[mask]
    if len(p) == 0:
        raise ValueError("empty mask")
    pred_labels = np.unique(p)
    true_labels = np.unique(t)
    pi = {int(v): i for i, v in enumerate(pred_labels)}
    ti = {int(v): i for i, v in enumerate(true_labels)}
    table = np.zeros((len(pred_labels), len(true_labels)), dtype=np.int64)
    np.add.at(table, ([pi[int(a)] for a in p], [ti[int(b)] for b in t]), 1)
    # noise predictions are only ever correct on true background, and a cell
    # prediction never receives credit on background molecules
    if 0 in pi and 0 in ti:
        corner = table[pi[0], ti[0]]
        table[pi[0], :] = 0
        table[:, ti[0]] = 0
        table[pi[0], ti[0]] = corner
    else:
        if 0 in pi:
            table[pi[0], :] = 0
        if 0 in ti:
            table[:, ti[0]] = 0
    row, col = linear_sum_assignment(-table)
    correct = int(sum(int(table[r, c]) for r, c in zip(row, col)))
    return {
        "accuracy": float(correct / len(p)),
        "n_evaluated": int(len(p)),
        "n_matched_labels": int(len(row)),
        "n_predicted_labels": int(len(pred_labels)),
        "n_true_labels": int(len(true_labels)),
        "n_predicted_zero": int(np.sum(p == 0)),
    }


def _read_predictions(out_dir: Path, df: pd.DataFrame) -> np.ndarray:
    """Per-molecule predicted cell from ``segmentation.csv``, aligned to ``df``.

    Labels ``cell_<n>`` become integers; noise / empty becomes 0.  Alignment
    prefers row order (verified against coordinates), falling back to a
    coordinate merge.
    """
    seg = pd.read_csv(out_dir / "segmentation.csv")
    if len(seg) != len(df):
        raise ValueError(f"segmentation has {len(seg)} rows, input has {len(df)}")
    lab = seg["cell"].astype(str)
    noise = seg.get("is_noise", pd.Series(0, index=seg.index)).astype(bool)
    num = lab.str.extract(r"(\d+)", expand=False)
    pred = np.where(noise | num.isna(), 0, pd.to_numeric(num, errors="coerce")
                    .fillna(0)).astype(np.int64)
    x, y = seg["x"].to_numpy(float), seg["y"].to_numpy(float)
    dx, dy = df["x"].to_numpy(float), df["y"].to_numpy(float)
    if np.allclose(x, dx, rtol=1e-5, atol=1e-3) and np.allclose(y, dy, rtol=1e-5, atol=1e-3):
        return pred  # same row order as the input
    # fallback: merge on rounded coordinates (unique by construction here)
    key = lambda a, b: np.round(a, 3) + 1j * np.round(b, 3)  # noqa: E731
    order = np.argsort(key(x, y), kind="stable")
    want = np.argsort(key(dx, dy), kind="stable")
    aligned = np.empty_like(pred)
    aligned[want] = pred[order]
    if not (np.allclose(x[order], dx[want], rtol=1e-5, atol=1e-3)
            and np.allclose(y[order], dy[want], rtol=1e-5, atol=1e-3)):
        raise ValueError("could not align segmentation rows to the input")
    return aligned


def run_dataset(dataset_id: str, binary: str, run_root: Path) -> dict:
    data_dir = common.data_root() / "sim" / dataset_id
    meta = json.loads((data_dir / "meta.json").read_text())
    df = pd.read_parquet(data_dir / "molecules.parquet")
    b = meta["baysor"]

    out_dir = run_root / dataset_id
    out_dir.mkdir(parents=True, exist_ok=True)
    prior = b["prior"]
    if prior not in ("column", "none"):
        raise ValueError(f"{dataset_id}: unsupported baysor.prior {prior!r}")
    cmd = [binary, "run", str(data_dir / "molecules.parquet")]
    if prior == "column":
        cmd.append(":prior")
    cmd += ["-x", "x", "-y", "y", "-g", "gene",
            "-s", f"{b['scale_um']:.4f}",
            "-m", str(b["min_molecules_per_cell"])]
    if prior == "column":
        cmd += ["--prior-segmentation-confidence", str(b["prior_confidence"])]
    cmd += ["-o", str(out_dir) + os.sep]
    if "z" in df.columns:
        cmd += ["-z", "z"]
    cmd += list(b.get("extra_args") or [])

    env = dict(os.environ,
               OMP_NUM_THREADS=str(MAX_THREADS),
               OPENBLAS_NUM_THREADS=str(MAX_THREADS),
               MKL_NUM_THREADS=str(MAX_THREADS))
    t0 = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    wall = time.monotonic() - t0

    record = {
        "id": dataset_id,
        "command": cmd,
        "prior_mode": prior,
        "threads_env": {k: env[k] for k in ("OMP_NUM_THREADS",)},
        "exit_code": proc.returncode,
        "wall_time_s": round(wall, 1),
        "outputs": sorted(p.name for p in out_dir.glob("segmentation*")),
    }
    if proc.returncode != 0:
        record["error"] = (proc.stderr or proc.stdout)[-2000:]
        return record
    if not (out_dir / "segmentation.csv").exists():
        record["error"] = "no segmentation.csv produced"
        return record

    pred = _read_predictions(out_dir, df)
    interior = df["interior"].to_numpy(bool)
    true = df["cell"].to_numpy(np.int64)
    record["n_output_molecules"] = int(len(pred))
    record["assignment_accuracy_interior"] = majority_match_accuracy(
        pred, true, interior)
    record["one_to_one_accuracy_interior"] = hungarian_match_accuracy(
        pred, true, interior)
    record["assignment_accuracy_all"] = majority_match_accuracy(
        pred, true, np.ones(len(true), bool))
    record["one_to_one_accuracy_all"] = hungarian_match_accuracy(
        pred, true, np.ones(len(true), bool))
    truth = meta.get("truth") or {}
    if truth.get("oracle_accuracy") is not None:
        record["oracle_accuracy"] = truth["oracle_accuracy"]
        record["naive_accuracy"] = truth.get("naive_accuracy")
    return record


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--binary", default=DEFAULT_BINARY)
    p.add_argument("--ids", nargs="+", default=DEFAULT_IDS)
    p.add_argument("--run-root", default=None,
                   help="default: $BAYSOR_BENCH_DATA/runs/bench-sim-sanity")
    p.add_argument("--report", default=str(REPORT))
    args = p.parse_args(argv)

    run_root = (Path(args.run_root) if args.run_root
                else common.data_root() / "runs" / "bench-sim-sanity")
    records = []
    for dataset_id in args.ids:
        print(f"== {dataset_id}")
        rec = run_dataset(dataset_id, args.binary, run_root)
        records.append(rec)
        acc = rec.get("assignment_accuracy_interior", {}).get("accuracy")
        o2o = rec.get("one_to_one_accuracy_interior", {}).get("accuracy")
        print(f"   exit={rec['exit_code']} wall={rec['wall_time_s']}s "
              f"maj-acc={acc} one-to-one={o2o} "
              f"oracle={rec.get('oracle_accuracy')}")

    report = {
        "binary": args.binary,
        "max_threads": MAX_THREADS,
        "run_root": str(run_root),
        "datasets": records,
    }
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {args.report}")
    ok = all(r["exit_code"] == 0 and "error" not in r for r in records)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
