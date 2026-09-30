#!/usr/bin/env python3
"""Create and list benchmark baselines.

``baseline.py create --run-id R --name NAME``

* copies each dataset's ``metrics.json`` from ``runs/<R>/`` into
  ``$BAYSOR_BENCH_DATA/baselines/NAME/<dataset>.json`` (local, never
  committed), and
* keeps the per-replicate assignment tables, the saved cell types
  (``celltypes.parquet``) and the fixed audit pair set
  (``fixed_pairs.json``) under ``$BAYSOR_BENCH_DATA/baselines/NAME/<dataset>/``,
  recording their sha256 in the metric JSON.

Only *successful* replicates count towards the noise floor. Because Baysor is
stochastic above one thread, a normal baseline requires a run with >= 3
successful replicates so per-metric mean and SD (the noise floor) can be
stored; use ``--allow-incomplete`` to override (fixtures, smoke baselines).

``--identical`` creates the exact baseline flavour for ``--expect
identical``: a 1-thread run (bitwise-deterministic), 1+ replicates, recorded
as ``flavour: identical``.

Everything is staged in temp directories and swapped in atomically at the
end: an error (or ``--force`` over an existing baseline) never leaves stale
files behind and never deletes the old baseline unless the new one is fully
in place.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

NAME_RE = re.compile(r"[A-Za-z0-9._-]+")
MIN_REPLICATES = 3


def _collect_dataset(run_root: Path, ds_root: Path, mf: Path, data_out: Path,
                     name: str, run_id: str, identical: bool,
                     allow_incomplete: bool, problems: list[str]) -> Optional[dict]:
    """Stage one dataset's baseline files under ``data_out`` and return the
    amended metrics dict (or None, appending to ``problems``)."""
    m = common.read_json(mf)
    ds_id = m["dataset"]["id"]
    ok_reps = [r for r in m.get("reps", []) if r.get("status") == "ok"]
    n_ok = len(ok_reps)

    if identical:
        threads = m.get("threads")
        if threads != 1:
            problems.append(f"{ds_id}: --identical requires a 1-thread run "
                            f"(this run used threads={threads})")
            return None
        if n_ok < 1:
            problems.append(f"{ds_id}: no successful replicates")
            return None
    elif n_ok < MIN_REPLICATES and not allow_incomplete:
        problems.append(
            f"{ds_id}: only {n_ok} successful replicate(s); a noise floor needs >= "
            f"{MIN_REPLICATES} (use --allow-incomplete for deterministic/"
            f"fixture runs)")
        return None

    # copy assignment tables of successful replicates, record their sha256
    sha_by_rep: dict[str, str] = {}
    for rec in ok_reps:
        src = ds_root / f"rep{rec['rep']}" / "assignment.parquet"
        if not src.is_file():
            problems.append(f"{ds_id}: missing assignment {src}")
            continue
        dst = data_out / f"rep{rec['rep']}" / "assignment.parquet"
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        sha = common.sha256_file(dst)
        recorded = rec.get("assignment_sha256")
        if recorded and recorded != sha:
            problems.append(f"{ds_id} rep{rec['rep']}: sha256 mismatch "
                            f"({recorded} != {sha})")
            continue
        sha_by_rep[str(rec["rep"])] = sha

    # saved cell types (rep0's audit typing: quick-cluster anchor or transfer)
    celltypes_sha = None
    src_ct = ds_root / "rep0" / "celltypes.parquet"
    if src_ct.is_file():
        dst_ct = data_out / "celltypes.parquet"
        dst_ct.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_ct, dst_ct)
        celltypes_sha = common.sha256_file(dst_ct)

    # fixed pair set taken from rep0's audit (pairs_top), used by later runs
    # so total_admixture_rate is comparable across runs
    fixed_sha, n_fixed = None, None
    src_audit = ds_root / "rep0" / "celladmix.json"
    if src_audit.is_file():
        try:
            audit = common.read_json(src_audit)
        except (OSError, ValueError):
            audit = None
        if isinstance(audit, dict):
            pairs = [{"source": str(p["source"]), "target": str(p["target"])}
                     for p in (audit.get("pairs_top") or [])
                     if isinstance(p, dict) and "source" in p and "target" in p]
            if pairs:
                dst_fp = data_out / "fixed_pairs.json"
                dst_fp.parent.mkdir(parents=True, exist_ok=True)
                with open(dst_fp, "w") as fh:
                    json.dump(pairs, fh, indent=2)
                    fh.write("\n")
                fixed_sha = common.sha256_file(dst_fp)
                n_fixed = len(pairs)

    m["baseline"] = {
        "name": name,
        "run_id": run_id,
        "created": common.utc_now(),
        "flavour": "identical" if identical else "noise_floor",
        "identical": bool(identical),
        "noise_floor_replicates": n_ok,
        "noise_floor_valid": n_ok >= MIN_REPLICATES,
        "assignments_dir": f"$BAYSOR_BENCH_DATA/baselines/{name}/{ds_id}/",
        "assignment_sha256": sha_by_rep,
        "celltypes": (f"$BAYSOR_BENCH_DATA/baselines/{name}/{ds_id}/celltypes.parquet"
                      if celltypes_sha else None),
        "celltypes_sha256": celltypes_sha,
        "fixed_pairs": (f"$BAYSOR_BENCH_DATA/baselines/{name}/{ds_id}/fixed_pairs.json"
                        if fixed_sha else None),
        "fixed_pairs_sha256": fixed_sha,
        "n_fixed_pairs": n_fixed,
    }
    return m


def _commit(swaps: list[tuple[Path, Path]], force: bool) -> None:
    """Rename each ``tmp`` into ``final``, rolling back on any error.

    Existing ``final`` directories are moved aside first and only deleted
    once every swap succeeded, so a failure never loses the old baseline.
    """
    backups: list[tuple[Path, Path]] = []
    placed: list[Path] = []
    try:
        for tmp, final in swaps:
            final.parent.mkdir(parents=True, exist_ok=True)
            if final.exists():
                if not force:
                    raise RuntimeError(f"{final} already exists; use --force")
                old = final.with_name(f".{final.name}.old-{os.getpid()}")
                shutil.rmtree(old, ignore_errors=True)
                final.rename(old)
                backups.append((old, final))
            tmp.rename(final)
            placed.append(final)
    except BaseException:
        for final in placed:
            shutil.rmtree(final, ignore_errors=True)
        for old, final in backups:
            if old.exists():
                old.rename(final)
        raise
    for old, _ in backups:
        shutil.rmtree(old, ignore_errors=True)


def create(run_id: str, name: str, root: Path, baselines_dir: Path,
           allow_incomplete: bool = False, force: bool = False,
           identical: bool = False) -> int:
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
    # a directory holding only hand-written docs (README.md) is not a
    # baseline yet: creating one over it must work without --force
    docs_only = out_dir.is_dir() and not any(out_dir.glob("*.json"))
    if out_dir.exists() and not docs_only and not force:
        print(f"error: baseline '{name}' already exists at {out_dir}; use --force",
              file=sys.stderr)
        return 2

    assignments_root = root / "baselines" / name
    # the default layout stores metric JSONs and assignment data in the same
    # directory (<data-root>/baselines/NAME); keep them separable for
    # --baselines-dir overrides (tests, fixtures)
    try:
        merged = baselines_dir.resolve() == (root / "baselines").resolve()
    except OSError:
        merged = baselines_dir == root / "baselines"
    pid = os.getpid()
    if merged:
        json_tmp = data_tmp = baselines_dir / f".{name}.tmp-{pid}"
        shutil.rmtree(json_tmp, ignore_errors=True)
        json_tmp.mkdir(parents=True, exist_ok=True)
    else:
        json_tmp = baselines_dir / f".{name}.tmp-{pid}"
        data_tmp = root / "baselines" / f".{name}.tmp-{pid}"
        for d in (json_tmp, data_tmp):
            shutil.rmtree(d, ignore_errors=True)
            d.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    problems: list[str] = []
    # carry over hand-written docs (README.md, SUMMARY.md, resources.csv,
    # ...) across recreations: the atomic swap replaces the whole directory
    if out_dir.is_dir():
        for extra in out_dir.iterdir():
            if extra.is_file() and extra.suffix != ".json":
                shutil.copy2(extra, json_tmp / extra.name)
    try:
        for mf in metrics_files:
            ds_id = common.read_json(mf)["dataset"]["id"]
            data_out = data_tmp / ds_id
            m = _collect_dataset(run_root, run_root / ds_id, mf, data_out,
                                 name, run_id, identical, allow_incomplete,
                                 problems)
            if m is None:
                continue
            out_path = json_tmp / f"{ds_id}.json"
            common.write_json(out_path, m)
            written.append(out_path)

        if problems:
            for p in problems:
                print(f"error: {p}", file=sys.stderr)
            print("error: baseline not modified (old baseline kept)", file=sys.stderr)
            return 2
        if not written:
            print("error: nothing written", file=sys.stderr)
            return 2

        swaps = ([(json_tmp, out_dir)] if merged else
                 [(json_tmp, out_dir), (data_tmp, assignments_root)])
        _commit(swaps, force or docs_only)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        print("error: baseline not modified (old baseline kept)", file=sys.stderr)
        return 2
    finally:
        shutil.rmtree(json_tmp, ignore_errors=True)
        shutil.rmtree(data_tmp, ignore_errors=True)

    print(f"baseline '{name}': {len(written)} dataset(s) -> {out_dir}")
    for p in sorted(out_dir.glob("*.json")):
        m = common.read_json(p)
        bl = m["baseline"]
        print(f"  {m['dataset']['id']}: {bl['noise_floor_replicates']} ok rep(s), "
              f"{m['dataset'].get('n_molecules')} molecules, "
              f"flavour={bl['flavour']}, "
              f"assign sha256: {list(bl['assignment_sha256'].values())[:1]}"
              f"{'...' if len(bl['assignment_sha256']) > 1 else ''}"
              f"{', celltypes+fixed-pairs stored' if bl.get('celltypes_sha256') else ''}")
    print(f"baseline data: {assignments_root}")
    return 0


def create_suite(suite_name: str, name: str, run_id: str, root: Path,
                 manifest: Optional[str] = None) -> int:
    """Freeze every run-id group of a suite run as a baseline
    (``bench.sh --suite ... --create-baseline NAME``).

    One baseline per group: the multi-thread group as ``NAME``, the
    1-thread bitwise group as ``NAME-t1`` (identical flavour). Each group
    is created with ``--force`` (an explicit (re)create, swapped in
    atomically), ``--allow-incomplete`` when its steps run fewer than 3
    replicates, and ``--identical`` when its steps expect ``identical``.
    Returns the first non-zero ``create`` exit code (2 on setup errors).
    """
    import suites
    try:
        suite = suites.resolve(suite_name, manifest)
        run_ids = suites.group_run_ids(suite, run_id)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    baselines_dir = root / "baselines"
    seen: dict[str, str] = {}
    rc = 0
    for group, rid in run_ids.items():
        steps = [s for s in suite.steps if s.group == group]
        threads = {s.threads for s in steps}
        if len(threads) != 1:
            print(f"error: group '{group}' mixes thread counts "
                  f"{sorted(threads)}; cannot derive a baseline name -- "
                  "create the baselines with baseline.py create yourself",
                  file=sys.stderr)
            return 2
        t = threads.pop()
        bname = name if t != 1 else f"{name}-t1"
        if bname in seen:
            print(f"error: groups '{seen[bname]}' and '{group}' both map to "
                  f"baseline '{bname}' -- create them with baseline.py "
                  "create yourself", file=sys.stderr)
            return 2
        seen[bname] = group
        identical = any(s.expect == "identical" for s in steps)
        reps = min(s.replicates for s in steps)
        print(f"== baseline {bname}: group {group} (run {rid}, {t} "
              f"thread(s) x {reps} rep(s)"
              f"{', identical' if identical else ''}) ==")
        r = create(rid, bname, root, baselines_dir,
                   allow_incomplete=reps < MIN_REPLICATES, force=True,
                   identical=identical)
        if r != 0:
            rc = r
    return rc


def list_baselines(baselines_dir: Path) -> int:
    if not baselines_dir.is_dir():
        print(f"no baselines directory at {baselines_dir}")
        return 0
    names = sorted(d.name for d in baselines_dir.iterdir()
                   if d.is_dir() and not d.name.startswith("."))
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
                  f"valid_floor={bl.get('noise_floor_valid')} "
                  f"flavour={bl.get('flavour', 'noise_floor')}")
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
                   help="default <data-root>/baselines")
    c.add_argument("--allow-incomplete", action="store_true",
                   help=f"allow baselines with fewer than {MIN_REPLICATES} "
                        f"successful replicates")
    c.add_argument("--force", action="store_true", help="overwrite an existing baseline")
    c.add_argument("--identical", action="store_true",
                   help="exact baseline flavour for --expect identical: requires "
                        "a 1-thread run (bitwise-deterministic), 1+ replicates")

    l = sub.add_parser("list", help="list baselines")
    l.add_argument("--baselines-dir", default=None)
    l.add_argument("--data-root", default=None)

    s = sub.add_parser("create-suite",
                       help="freeze every group of a suite run as baselines "
                            "(multi-thread group as NAME, 1-thread bitwise "
                            "group as NAME-t1)")
    s.add_argument("--suite", required=True)
    s.add_argument("--name", required=True)
    s.add_argument("--run-id", required=True)
    s.add_argument("--data-root", default=None)
    s.add_argument("--manifest", default=None,
                   help="suite manifest path (default "
                        "<repo>/benchmarks/datasets/suites.yaml)")

    args = ap.parse_args(argv)
    baselines_dir = Path(args.baselines_dir) \
        if getattr(args, "baselines_dir", None) \
        else common.baselines_root(getattr(args, "data_root", None))

    if args.cmd == "create-suite":
        root = common.data_root(args.data_root)
        return create_suite(args.suite, args.name, args.run_id, root,
                            manifest=args.manifest)
    if args.cmd == "create":
        root = common.data_root(args.data_root)
        return create(args.run_id, args.name, root, baselines_dir,
                      allow_incomplete=args.allow_incomplete, force=args.force,
                      identical=args.identical)
    return list_baselines(baselines_dir)


if __name__ == "__main__":
    raise SystemExit(main())
