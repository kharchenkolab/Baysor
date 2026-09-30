#!/usr/bin/env python3
"""Resolve benchmark suites from the manifest ``benchmarks/datasets/suites.yaml``.

A *suite* is an ordered list of *steps*; every step is one ``run.py``
invocation (dataset spec, threads, replicates, timeout, cellAdmix audit,
metric options) with an expected comparison mode (``expect``) against a
``baseline``. Steps that share a ``group`` execute in the same run folder,
so a suite produces one run folder per group (see ``group_run_ids``); the
folder of the group holding the suite's ``identical`` step keeps the base
run-id verbatim, because at 1 thread Baysor's output depends on the output
path length (run-ids must stay <= 17 characters, see the harness README).

Manifest schema (all keys except ``name``/``datasets``/``expect``/
``baseline`` have defaults)::

    suites:
      regular:
        description: ...
        budget_minutes: 20
        resources: baselines/bugfixes-35e8a7e/resources.csv   # estimates
        steps:
          - name: exact
            group: exact                # default: the step name
            datasets: quick|full|all|<ids/globs>
            threads: 1                  # default 6
            replicates: 1               # default 1
            timeout: 1800               # seconds per replicate, default 1800
            celladmix: false            # default true
            celltypes_from: baseline    # optional, default none
            no_ami: true                # skip informational AMI, default false
            expect: identical           # identical | same | improved
            baseline: bugfixes-35e8a7e-t1

Usage::

    suites.py --list                       # available suites
    suites.py --suite regular              # resolved plan + time estimate
    suites.py --suite release --json       # machine-readable plan
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common                      # noqa: E402
import resources as resmod         # noqa: E402

DEFAULT_MANIFEST = "benchmarks/datasets/suites.yaml"
EXPECTS = ("identical", "same", "improved")
_RUN_ID_RE = re.compile(r"[A-Za-z0-9._-]+")


# ---------------------------------------------------------------------------
# manifest loading
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Step:
    name: str
    group: str
    datasets: str
    threads: int
    replicates: int
    timeout: float
    celladmix: bool
    celltypes_from: Optional[str]
    no_ami: bool
    expect: str
    baseline: str


@dataclass(frozen=True)
class Suite:
    name: str
    description: str
    budget_minutes: Optional[float]
    resources: Optional[str]      # CSV path relative to benchmarks/
    steps: tuple


def default_manifest(repo: Optional[Path] = None) -> Path:
    repo = repo or common.repo_root()
    return repo / DEFAULT_MANIFEST


def load(path: Optional[Path] = None) -> dict:
    path = Path(path) if path else default_manifest()
    if not path.is_file():
        raise ValueError(f"suite manifest not found: {path}")
    with open(path) as fh:
        doc = yaml.safe_load(fh) or {}
    if not isinstance(doc, dict) or not isinstance(doc.get("suites"), dict):
        raise ValueError(f"{path}: expected a top-level 'suites' mapping")
    return doc


def _step_from(entry: dict, index: int) -> Step:
    if not isinstance(entry, dict):
        raise ValueError(f"step #{index + 1} is not a mapping")
    name = str(entry.get("name") or "")
    if not name:
        raise ValueError(f"step #{index + 1} lacks a 'name'")
    if "datasets" not in entry:
        raise ValueError(f"step '{name}' lacks 'datasets'")
    ds = entry["datasets"]
    spec = ds if isinstance(ds, str) else ",".join(str(t) for t in ds)
    expect = str(entry.get("expect") or "")
    if expect not in EXPECTS:
        raise ValueError(f"step '{name}': 'expect' must be one of "
                         f"{', '.join(EXPECTS)} (got {expect!r})")
    if not entry.get("baseline"):
        raise ValueError(f"step '{name}' lacks 'baseline'")
    group = str(entry.get("group") or name)
    if not _RUN_ID_RE.fullmatch(group):
        raise ValueError(f"step '{name}': group {group!r} is not run-id safe")
    return Step(
        name=name,
        group=group,
        datasets=spec,
        threads=int(entry.get("threads", 6)),
        replicates=int(entry.get("replicates", 1)),
        timeout=float(entry.get("timeout", 1800)),
        celladmix=bool(entry.get("celladmix", True)),
        celltypes_from=(str(entry["celltypes_from"])
                        if entry.get("celltypes_from") else None),
        no_ami=bool(entry.get("no_ami", False)),
        expect=expect,
        baseline=str(entry["baseline"]),
    )


def resolve(name: str, manifest_path: Optional[str] = None) -> Suite:
    """Load the manifest and validate one suite; ValueError on any problem."""
    path = Path(manifest_path) if manifest_path else default_manifest()
    doc = load(path)
    suites = doc["suites"]
    if name not in suites:
        raise ValueError(f"unknown suite {name!r} in {path} "
                         f"(available: {', '.join(sorted(suites))})")
    entry = suites[name] or {}
    if not isinstance(entry, dict):
        raise ValueError(f"{path}: suite '{name}' is not a mapping")
    raw_steps = entry.get("steps")
    if not isinstance(raw_steps, list) or not raw_steps:
        raise ValueError(f"{path}: suite '{name}' has no steps")
    steps = tuple(_step_from(s, i) for i, s in enumerate(raw_steps))
    names = [s.name for s in steps]
    if len(set(names)) != len(names):
        raise ValueError(f"{path}: suite '{name}' has duplicate step names")
    # steps sharing a group must be comparable in one run folder
    for group in dict.fromkeys(s.group for s in steps):
        g = [s for s in steps if s.group == group]
        if len({(s.expect, s.baseline) for s in g}) != 1:
            raise ValueError(
                f"{path}: suite '{name}', group '{group}': steps must share "
                f"one expect/baseline pair (got "
                f"{sorted({(s.expect, s.baseline) for s in g})})")
    return Suite(
        name=name,
        description=str(entry.get("description") or ""),
        budget_minutes=(float(entry["budget_minutes"])
                        if entry.get("budget_minutes") is not None else None),
        resources=(str(entry["resources"]) if entry.get("resources") else None),
        steps=steps,
    )


def list_suites(manifest_path: Optional[str] = None) -> list[str]:
    doc = load(Path(manifest_path) if manifest_path else None)
    return sorted(doc["suites"])


def step_spec(step: Step) -> str:
    """The ``--datasets`` spec of a step (tier name, ids and/or globs)."""
    return step.datasets


# ---------------------------------------------------------------------------
# run-id groups
# ---------------------------------------------------------------------------

def group_run_ids(suite: Suite, base_run_id: str) -> dict[str, str]:
    """Map each step group of the suite to its run folder under ``runs/``.

    One group -> the base run-id verbatim. Multiple groups -> the group
    holding the suite's (first) ``identical`` step keeps the base run-id
    (1-thread determinism depends on the output-path length), every other
    group gets ``<base>-<group>``.
    """
    if not _RUN_ID_RE.fullmatch(base_run_id):
        raise ValueError(f"run-id {base_run_id!r} may only contain letters, "
                         f"digits, '.', '_', '-'")
    groups = list(dict.fromkeys(s.group for s in suite.steps))
    if len(groups) == 1:
        return {groups[0]: base_run_id}
    primary = next((g for g in groups
                    if any(s.expect == "identical" for s in suite.steps
                           if s.group == g)), groups[0])
    out = {}
    for g in groups:
        rid = base_run_id if g == primary else f"{base_run_id}-{g}"
        if not _RUN_ID_RE.fullmatch(rid):
            raise ValueError(f"group '{g}' produces an invalid run-id {rid!r}")
        out[g] = rid
    return out


# ---------------------------------------------------------------------------
# dataset expansion / estimation
# ---------------------------------------------------------------------------

def expand_tolerant(root: Path, spec: str) -> list[str]:
    """Dataset ids of a spec; unmatched tokens are skipped (used by the
    inventory, which must not fail on a partial data root)."""
    all_ds = common.discover_datasets(root)
    by_id = {d.id: d for d in all_ds}
    out: list[str] = []
    for tok in (t.strip() for t in (spec or "").split(",")):
        if not tok:
            continue
        if tok in ("quick", "full", "all"):
            matched = [d for d in all_ds if tok == "all" or d.tier == tok]
        elif tok in by_id:
            matched = [by_id[tok]]
        else:
            matched = [d for d in all_ds if fnmatch.fnmatch(d.id, tok)]
        out.extend(d.id for d in matched)
    return list(dict.fromkeys(out))


def expand_strict(root: Path, spec: str) -> list[str]:
    """Expand a spec with ``common.select_datasets`` semantics (raises on
    unknown ids — this is what dry-run validation uses)."""
    return [d.id for d in common.select_datasets(root, spec)]


def resources_path(suite: Suite, repo: Optional[Path] = None) -> Optional[Path]:
    if not suite.resources:
        return None
    repo = repo or common.repo_root()
    return repo / "benchmarks" / suite.resources


def estimate(suite: Suite, root: Path,
             resources_csv: Optional[Path]) -> dict:
    """Estimated wall/CPU seconds per step from the committed resources CSV.

    Wall = measured Baysor wall per replicate x replicates + cellAdmix audit
    time (steps with ``celladmix``); CPU = measured CPU time (user+sys) the
    same way, and wall as a proxy for CPU on 1-thread steps (~100 % usage).
    Metrics/typing bookkeeping (seconds per dataset) is not included.
    """
    table = resmod.load_csv(resources_csv) if resources_csv else {}
    steps = []
    totals = {"wall_s": 0.0, "cpu_s": 0.0, "resources_found": True}
    for step in suite.steps:
        wall = cpu = 0.0
        missing = []
        ids = expand_strict(root, step.datasets)
        for ds_id in ids:
            row = table.get(ds_id) or {}
            one_wall = row.get("wall1_s" if step.threads == 1
                               else "wall6_mean_s")
            one_cpu = row.get("wall1_s" if step.threads == 1
                              else "cpu6_mean_s")
            if one_wall is None:
                missing.append(ds_id)
                continue
            wall += one_wall * step.replicates
            cpu += (one_cpu if one_cpu is not None else one_wall) \
                * step.replicates
            audit = row.get("audit_wall_s")
            if step.celladmix and audit:
                wall += audit * step.replicates
                cpu += audit * step.replicates
        if missing:
            totals["resources_found"] = False
        steps.append({"step": step.name, "group": step.group,
                      "datasets": len(ids), "threads": step.threads,
                      "replicates": step.replicates,
                      "wall_s": wall, "cpu_s": cpu,
                      "missing_resources": missing})
        totals["wall_s"] += wall
        totals["cpu_s"] += cpu
    return {"steps": steps, "total": totals,
            "resources_csv": str(resources_csv) if resources_csv else None}


def plan_text(suite: Suite, base_run_id: str, selections: dict,
              resources_csv: Optional[Path] = None,
              run_ids: Optional[dict] = None,
              root: Optional[Path] = None) -> str:
    """Human-readable resolution of a suite: steps, run-ids, dataset lists
    and the estimated wall/CPU time from the resources CSV."""
    run_ids = run_ids or group_run_ids(suite, base_run_id)
    root = root or common.data_root()
    est = estimate(suite, root, resources_csv) if resources_csv else None
    lines = [f"suite {suite.name}: {suite.description}"]
    if suite.budget_minutes:
        lines.append(f"  design budget: <= {suite.budget_minutes:g} min wall")
    if est and not est["total"]["resources_found"]:
        lines.append("  WARNING: no resource data for some datasets "
                     "(estimates incomplete; regenerate resources.csv)")
    for i, step in enumerate(suite.steps):
        ids = selections.get(step.name)
        lines.append(
            f"  step {step.name}: group={step.group} "
            f"run-id={run_ids[step.group]} threads={step.threads} "
            f"replicates={step.replicates} timeout={step.timeout:g}s "
            f"celladmix={'on' if step.celladmix else 'off'} "
            f"celltypes-from={step.celltypes_from or '-'} "
            f"ami={'skipped' if step.no_ami else 'computed'} "
            f"expect={step.expect} baseline={step.baseline}")
        if ids is not None:
            lines.append(f"    datasets ({len(ids)}): " + ", ".join(ids))
        if est:
            e = est["steps"][i]
            lines.append(
                f"    estimated: wall ~{resmod.fmt_duration(e['wall_s'])}, "
                f"cpu ~{resmod.fmt_duration(e['cpu_s'])}"
                + (f" ({len(e['missing_resources'])} datasets without "
                   f"resource data)" if e["missing_resources"] else ""))
    if est:
        t = est["total"]
        lines.append(f"  estimated total: wall ~{resmod.fmt_duration(t['wall_s'])}, "
                     f"cpu ~{resmod.fmt_duration(t['cpu_s'])}"
                     f" (from {est['resources_csv']}; excludes metrics/"
                     f"typing bookkeeping)")
    return "\n".join(lines)


def membership(root: Path, repo: Optional[Path] = None,
               manifest_path: Optional[str] = None) -> dict[str, set[str]]:
    """Per-suite dataset-id membership, tolerant of a partial data root."""
    path = Path(manifest_path) if manifest_path else default_manifest(repo)
    if not path.is_file():
        return {}
    doc = load(path)
    out: dict[str, set[str]] = {}
    for name in doc["suites"]:
        suite = resolve(name, str(path))
        ids: set[str] = set()
        for step in suite.steps:
            ids.update(expand_tolerant(root, step.datasets))
        out[name] = ids
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=None,
                    help=f"manifest path (default <repo>/{DEFAULT_MANIFEST})")
    ap.add_argument("--list", action="store_true", help="list available suites")
    ap.add_argument("--suite", default=None, help="suite name to resolve")
    ap.add_argument("--run-id", default="RUNID",
                    help="base run id used in the plan (default RUNID)")
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--json", action="store_true",
                    help="emit the plan as JSON")
    args = ap.parse_args(argv)

    try:
        if args.list or not args.suite:
            for name in list_suites(args.manifest):
                print(name)
            return 0
        suite = resolve(args.suite, args.manifest)
        root = common.data_root(args.data_root)
        selections = {s.name: expand_strict(root, s.datasets)
                      for s in suite.steps}
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.json:
        run_ids = group_run_ids(suite, args.run_id)
        est = estimate(suite, root, resources_path(suite))
        print(json.dumps({
            "suite": suite.name,
            "description": suite.description,
            "budget_minutes": suite.budget_minutes,
            "run_ids": run_ids,
            "steps": [dict(vars(s), datasets=selections[s.name],
                           run_id=run_ids[s.group])
                      for s in suite.steps],
            "estimate": est,
        }, indent=2))
        return 0
    print(plan_text(suite, args.run_id, selections,
                    resources_csv=resources_path(suite), root=root))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
