#!/usr/bin/env python3
"""Create and list committed benchmark baselines.

``baseline.py create --run-id R --name NAME``

* copies each dataset's ``metrics.json`` from ``runs/<R>/`` into
  ``benchmarks/baselines/NAME/<dataset>.json`` (committed, small), and
* keeps the per-replicate assignment tables under
  ``$BAYSOR_BENCH_DATA/baselines/NAME/<dataset>/rep<k>/assignment.parquet``
  (never committed), recording their sha256 in the JSON.

Because Baysor is stochastic above one thread, a baseline normally requires a
run with >= 3 replicates so per-metric mean and SD (the noise floor) can be
stored; use ``--allow-incomplete`` to override (deterministic 1-thread runs,
fixtures, quick smoke baselines).
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

NAME_RE = re.compile(r"[A-Za-z0-9._-]+")
MIN_REPLICATES = 3


def create(run_id: str, name: str, root: Path, baselines_dir: Path,
           allow_incomplete: bool = False, force: bool = False) -> int:
    if not NAME_RE.fullmatch(name):
        print(f"error: invalid baseline name {name!r}", file=sys.stderr)
        return 2
    run_root = root / "runs" / run_id
    if not run_root.is_dir():
        print(f"error: run '{run_id}' not found under {run_root}", file=sys.stderr)
        return 2

    metrics_files = sorted(run_root.glob("*/metrics.json"))
    if not metrics_files:
        print(f"error: no */metrics.json under {run_root}", file=sys.stderr)
        return 2

    out_dir = baselines_dir / name
    if out_dir.exists() and any(out_dir.glob("*.json")) and not force:
        print(f"error: baseline '{name}' already exists at {out_dir}; use --force",
              file=sys.stderr)
        return 2

    assignments_root = root / "baselines" / name
    written, problems = [], []
    for mf in metrics_files:
        m = common.read_json(mf)
        ds_id = m["dataset"]["id"]
        reps = m.get("replicates", 0)
        ok_reps = [r for r in m.get("reps", []) if r.get("status") == "ok"]
        if reps < MIN_REPLICATES and not allow_incomplete:
            problems.append(
                f"{ds_id}: only {reps} replicate(s); a noise floor needs >= "
                f"{MIN_REPLICATES} (use --allow-incomplete for deterministic/"
                f"fixture runs)")
            continue

        # copy assignment tables, record their sha256
        ds_assign_root = assignments_root / ds_id
        if ds_assign_root.exists() and force:
            shutil.rmtree(ds_assign_root)
        sha_by_rep: dict[str, str] = {}
        for rec in ok_reps:
            src = run_root / ds_id / f"rep{rec['rep']}" / "assignment.parquet"
            if not src.is_file():
                problems.append(f"{ds_id}: missing assignment {src}")
                continue
            dst = ds_assign_root / f"rep{rec['rep']}" / "assignment.parquet"
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            sha = common.sha256_file(dst)
            recorded = rec.get("assignment_sha256")
            if recorded and recorded != sha:
                problems.append(f"{ds_id} rep{rec['rep']}: sha256 mismatch "
                                f"({recorded} != {sha})")
                continue
            sha_by_rep[str(rec["rep"])] = sha

        m["baseline"] = {
            "name": name,
            "run_id": run_id,
            "created": common.utc_now(),
            "noise_floor_replicates": len(ok_reps),
            "noise_floor_valid": len(ok_reps) >= MIN_REPLICATES,
            "assignments_dir": f"$BAYSOR_BENCH_DATA/baselines/{name}/{ds_id}/",
            "assignment_sha256": sha_by_rep,
        }
        out_path = out_dir / f"{ds_id}.json"
        common.write_json(out_path, m)
        written.append(out_path)

    if problems:
        for p in problems:
            print(f"error: {p}", file=sys.stderr)
        for p in written:
            p.unlink(missing_ok=True)
        return 2
    if not written:
        print("error: nothing written", file=sys.stderr)
        return 2

    print(f"baseline '{name}': {len(written)} dataset(s) -> {out_dir}")
    for p in written:
        m = common.read_json(p)
        n_ok = sum(1 for r in m.get("reps", []) if r.get("status") == "ok")
        print(f"  {m['dataset']['id']}: {n_ok} rep(s), "
              f"{m['dataset'].get('n_molecules')} molecules, "
              f"assign sha256: {list(m['baseline']['assignment_sha256'].values())[:1]}"
              f"{'...' if len(m['baseline']['assignment_sha256']) > 1 else ''}")
    print(f"assignment tables (not committed): {assignments_root}")
    return 0


def list_baselines(baselines_dir: Path) -> int:
    if not baselines_dir.is_dir():
        print(f"no baselines directory at {baselines_dir}")
        return 0
    names = sorted(d.name for d in baselines_dir.iterdir() if d.is_dir())
    if not names:
        print("no baselines")
        return 0
    for name in names:
        files = sorted((baselines_dir / name).glob("*.json"))
        print(f"{name}: {len(files)} dataset(s)")
        for f in files:
            m = common.read_json(f)
            bl = m.get("baseline", {})
            print(f"  {m['dataset']['id']}: run={bl.get('run_id')} "
                  f"reps={bl.get('noise_floor_replicates')} "
                  f"valid_floor={bl.get('noise_floor_valid')}")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("create", help="create a baseline from a run")
    c.add_argument("--run-id", required=True)
    c.add_argument("--name", required=True)
    c.add_argument("--data-root", default=None)
    c.add_argument("--baselines-dir", default=None,
                   help="default <repo>/benchmarks/baselines")
    c.add_argument("--allow-incomplete", action="store_true",
                   help=f"allow baselines with fewer than {MIN_REPLICATES} replicates")
    c.add_argument("--force", action="store_true", help="overwrite an existing baseline")

    l = sub.add_parser("list", help="list baselines")
    l.add_argument("--baselines-dir", default=None)

    args = ap.parse_args(argv)
    repo = common.repo_root()
    baselines_dir = Path(args.baselines_dir) if getattr(args, "baselines_dir", None) \
        else repo / "benchmarks" / "baselines"

    if args.cmd == "create":
        root = common.data_root(args.data_root)
        return create(args.run_id, args.name, root, baselines_dir,
                      allow_incomplete=args.allow_incomplete, force=args.force)
    return list_baselines(baselines_dir)


if __name__ == "__main__":
    raise SystemExit(main())
