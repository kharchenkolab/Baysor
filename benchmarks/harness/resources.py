#!/usr/bin/env python3
"""Extract per-dataset resource usage from finished benchmark runs.

Reads the ``/usr/bin/time -v`` block stored in every replicate's
``baysor.log`` (User/System time, Percent of CPU, Elapsed, Maximum resident
set size), the ``run.json`` status and, for real datasets, the cellAdmix
audit wall time from ``celladmix.json`` (``runtime_seconds.total``). No
Baysor rerun is involved: everything comes from runs that already exist
under ``$BAYSOR_BENCH_DATA/runs/``.

The result is ``$BAYSOR_BENCH_DATA/baselines/bugfixes-35e8a7e/resources.csv``
with one row per dataset:

============================  ====================================================
column                        meaning
============================  ====================================================
``dataset``                   dataset id
``molecules``, ``genes``      ``meta.json`` ``stats`` (for CPU per molecule)
``cpu6_mean_s``               6-thread CPU time (user + sys), mean over replicates
``cpu6_sd_s``                 sample SD of the same (0.0 for one replicate)
``wall6_mean_s``              6-thread wall time, mean over replicates
``peak_rss6_kb``              6-thread peak RSS, **max** over replicates
``wall1_s``                   1-thread wall time (1 replicate)
``peak_rss1_kb``              1-thread peak RSS
``audit_wall_s``              cellAdmix audit wall time, mean over replicates
                              (real datasets only; the audit RSS was never
                              measured and is not reported anywhere)
``cpu_s_per_1k_mol``          ``cpu6_mean_s`` per 1000 molecules
============================  ====================================================

Values that were never measured stay **empty** in the CSV and are rendered
as ``TODO`` by ``inventory.py``; they are never guessed or back-filled.

Usage::

    resources.py                    # rewrite the default CSV from the runs
    resources.py --out /tmp/r.csv   # write elsewhere
    resources.py --check            # exit 1 when the CSV is stale
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common                      # noqa: E402
from run import parse_time_v       # noqa: E402

DEFAULT_RUN6 = "benchbase-b"       # 6 threads, 3 replicates, all datasets
DEFAULT_RUN1 = "benchbase-t1"      # 1 thread, 1 replicate, quick tier
# output path relative to the data root ($BAYSOR_BENCH_DATA)
DEFAULT_OUT = "baselines/bugfixes-35e8a7e/resources.csv"

CSV_COLUMNS = ("dataset", "molecules", "genes",
               "cpu6_mean_s", "cpu6_sd_s", "wall6_mean_s", "peak_rss6_kb",
               "wall1_s", "peak_rss1_kb", "audit_wall_s", "cpu_s_per_1k_mol")


# ---------------------------------------------------------------------------
# aggregation helpers
# ---------------------------------------------------------------------------

def _mean_sd(values: list[float]) -> tuple[Optional[float], Optional[float]]:
    """(mean, sample SD); SD is 0.0 for a single value, None when empty."""
    vals = [v for v in values if v is not None]
    if not vals:
        return None, None
    mean = sum(vals) / len(vals)
    if len(vals) == 1:
        return mean, 0.0
    var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
    return mean, var ** 0.5


def _max(values: list[Optional[float]]) -> Optional[float]:
    vals = [v for v in values if v is not None]
    return max(vals) if vals else None


# ---------------------------------------------------------------------------
# run parsing
# ---------------------------------------------------------------------------

def read_replicate(rep_dir: Path) -> Optional[dict]:
    """Resource record of one replicate, or None when it cannot be used.

    The timing comes from the ``/usr/bin/time -v`` block in ``baysor.log``;
    ``run.json`` supplies the status (only ``ok`` replicates count) and
    serves as a fallback for wall/RSS when the log is gone.
    """
    run_json = rep_dir / "run.json"
    if not run_json.is_file():
        return None
    try:
        rec = common.read_json(run_json)
    except (OSError, ValueError):
        return None
    if rec.get("status") != "ok":
        return None
    out = {"user_s": None, "sys_s": None, "cpu_percent": None,
           "wall_s": None, "peak_rss_kb": None, "audit_wall_s": None}
    log = rep_dir / "baysor.log"
    if log.is_file():
        tv = parse_time_v(log.read_text(errors="replace"))
        out["user_s"] = tv["user_s"]
        out["sys_s"] = tv["sys_s"]
        out["cpu_percent"] = tv["cpu_percent"]
        out["wall_s"] = tv["wall_s"]
        out["peak_rss_kb"] = tv["peak_rss_kb"]
    if out["wall_s"] is None:
        out["wall_s"] = rec.get("wall_s_time_v", rec.get("wall_s"))
    if out["peak_rss_kb"] is None:
        out["peak_rss_kb"] = rec.get("peak_rss_kb")
    cam = rep_dir / "celladmix.json"
    if cam.is_file():
        try:
            total = (common.read_json(cam).get("runtime_seconds") or {}).get("total")
            out["audit_wall_s"] = float(total) if total is not None else None
        except (OSError, ValueError, TypeError):
            out["audit_wall_s"] = None
    return out


def _rep_dirs(run_dir: Optional[Path], dataset: str) -> list[Path]:
    if run_dir is None:
        return []
    base = run_dir / dataset
    if not base.is_dir():
        return []
    return sorted(p for p in base.glob("rep*") if p.is_dir())


def dataset_resources(ds_dir: Path, run6: Optional[Path],
                      run1: Optional[Path]) -> dict:
    """Aggregate resource usage of one dataset across the two runs."""
    meta = common.read_json(ds_dir / "meta.json")
    stats = meta.get("stats") or {}
    n_mol = stats.get("n_molecules")
    n_genes = stats.get("n_genes")

    reps6 = [r for p in _rep_dirs(run6, ds_dir.name)
             if (r := read_replicate(p)) is not None]
    cpu6 = [r["user_s"] + r["sys_s"] for r in reps6
            if r["user_s"] is not None and r["sys_s"] is not None]
    wall6 = [r["wall_s"] for r in reps6 if r["wall_s"] is not None]
    rss6 = [r["peak_rss_kb"] for r in reps6 if r["peak_rss_kb"] is not None]
    audits = [r["audit_wall_s"] for r in reps6
              if r["audit_wall_s"] is not None]

    reps1 = [r for p in _rep_dirs(run1, ds_dir.name)
             if (r := read_replicate(p)) is not None]
    wall1 = [r["wall_s"] for r in reps1 if r["wall_s"] is not None]
    rss1 = [r["peak_rss_kb"] for r in reps1 if r["peak_rss_kb"] is not None]

    cpu_mean, cpu_sd = _mean_sd(cpu6)
    per_1k = None
    if cpu_mean is not None and n_mol:
        per_1k = cpu_mean / (float(n_mol) / 1000.0)
    return {
        "dataset": ds_dir.name,
        "molecules": int(n_mol) if n_mol is not None else None,
        "genes": int(n_genes) if n_genes is not None else None,
        "cpu6_mean_s": cpu_mean,
        "cpu6_sd_s": cpu_sd,
        "wall6_mean_s": _mean_sd(wall6)[0],
        "peak_rss6_kb": _max(rss6),
        "wall1_s": _mean_sd(wall1)[0],
        "peak_rss1_kb": _max(rss1),
        "audit_wall_s": _mean_sd(audits)[0],
        "cpu_s_per_1k_mol": per_1k,
    }


def collect_rows(root: Path, run6: Optional[str], run1: Optional[str]) -> list[dict]:
    runs = root / "runs"
    run6_dir = runs / run6 if run6 else None
    run1_dir = runs / run1 if run1 else None
    if run6_dir is not None and not run6_dir.is_dir():
        raise FileNotFoundError(f"run not found: {run6_dir}")
    if run1_dir is not None and not run1_dir.is_dir():
        raise FileNotFoundError(f"run not found: {run1_dir}")
    ds_dirs = [d for k in ("sim", "real") for d in sorted((root / k).glob("*"))
               if d.is_dir() and (d / "meta.json").is_file()]
    rows = []
    for d in ds_dirs:
        # a dataset with no resources at all still gets a (mostly empty) row
        row = dataset_resources(d, run6_dir, run1_dir)
        if any(row[c] is not None for c in CSV_COLUMNS if c != "dataset"):
            rows.append(row)
    return sorted(rows, key=lambda r: r["dataset"])


# ---------------------------------------------------------------------------
# CSV I/O
# ---------------------------------------------------------------------------

_NUM_COLS = {c: "%.3f" for c in CSV_COLUMNS if c not in ("dataset", "molecules",
                                                        "genes", "peak_rss6_kb",
                                                        "peak_rss1_kb")}
_NUM_COLS["dataset"] = "%s"
_NUM_COLS["molecules"] = "%d"
_NUM_COLS["genes"] = "%d"
_NUM_COLS["peak_rss6_kb"] = "%d"
_NUM_COLS["peak_rss1_kb"] = "%d"


def render_csv(rows: list[dict]) -> str:
    lines = [",".join(CSV_COLUMNS)]
    for r in rows:
        cells = []
        for col in CSV_COLUMNS:
            v = r.get(col)
            if v is None:
                cells.append("")
            else:
                cells.append(_NUM_COLS[col] % v)
        lines.append(",".join(cells))
    return "\n".join(lines) + "\n"


def load_csv(path: Path) -> dict[str, dict]:
    """Parse a resources CSV into ``{dataset: {column: value|None}}``."""
    if not path.is_file():
        return {}
    out: dict[str, dict] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split(",")
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if not parts or not parts[0]:
                continue
            row = dict(zip(header, parts))
            rec: dict = {}
            for col in header:
                raw = row.get(col, "")
                if col == "dataset" or raw == "":
                    rec[col] = raw if col == "dataset" else None
                elif col in ("molecules", "genes", "peak_rss6_kb",
                             "peak_rss1_kb"):
                    rec[col] = int(raw)
                else:
                    rec[col] = float(raw)
            out[rec["dataset"]] = rec
    return out


# ---------------------------------------------------------------------------
# human-readable formatting (shared with inventory.py)
# ---------------------------------------------------------------------------

def fmt_duration(seconds: Optional[float]) -> str:
    """``36.9 s`` / ``16.4 min`` / ``TODO`` for missing values."""
    if seconds is None:
        return "TODO"
    s = float(seconds)
    if abs(s) < 60:
        return f"{s:.1f} s"
    return f"{s / 60:.1f} min"


def fmt_mean_sd(mean: Optional[float], sd: Optional[float]) -> str:
    """``77.4 ± 0.4 s`` — always seconds (SD stays readable), ``TODO``
    when the mean is missing."""
    if mean is None:
        return "TODO"
    if sd is None:
        return fmt_duration(mean)
    return f"{mean:.1f} ± {sd:.1f} s"


def fmt_bytes_kb(kb: Optional[float]) -> str:
    """``274.7 MB`` / ``3.2 GB`` / ``TODO`` for missing values."""
    if kb is None:
        return "TODO"
    mb = float(kb) / 1024.0
    if mb < 1024:
        return f"{mb:.1f} MB"
    return f"{mb / 1024:.2f} GB"


def fmt_wall_ram(wall: Optional[float], rss: Optional[float]) -> str:
    """Combined ``51.3 s / 268.1 MB`` cell for the 1-thread columns."""
    if wall is None and rss is None:
        return "TODO"
    return f"{fmt_duration(wall)} / {fmt_bytes_kb(rss)}"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=None,
                    help="data root (default $BAYSOR_BENCH_DATA or <repo>/.bench-data)")
    ap.add_argument("--run6", default=DEFAULT_RUN6,
                    help=f"6-thread run id (default {DEFAULT_RUN6})")
    ap.add_argument("--run1", default=DEFAULT_RUN1,
                    help=f"1-thread run id (default {DEFAULT_RUN1}; '' disables)")
    ap.add_argument("--out", default=None,
                    help=f"output CSV (default <data-root>/{DEFAULT_OUT})")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if the CSV is out of date")
    args = ap.parse_args(argv)

    root = common.data_root(args.data_root)
    out = Path(args.out) if args.out else root / DEFAULT_OUT

    try:
        rows = collect_rows(root, args.run6, args.run1 or None)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if not rows:
        print(f"error: no datasets with runs under {root}", file=sys.stderr)
        return 2
    text = render_csv(rows)

    if args.check:
        current = out.read_text() if out.is_file() else ""
        if current != text:
            print(f"{out} is out of date (rerun resources.py without --check)",
                  file=sys.stderr)
            return 1
        print(f"{out} is up to date ({len(rows)} datasets)")
        return 0

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"wrote {out} ({len(rows)} datasets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
