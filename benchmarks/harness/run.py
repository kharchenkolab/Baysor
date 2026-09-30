#!/usr/bin/env python3
"""Run a Baysor binary over benchmark datasets and collect metrics.

For every selected dataset and replicate this runner

* builds the ``baysor run`` command from the dataset's ``meta.json``,
* executes it under ``/usr/bin/time -v`` (wall time, CPU user/system time,
  percent of CPU, peak RSS) with ``OMP_NUM_THREADS`` pinned to ``--threads``,
* normalizes the segmentation output into ``assignment.parquet``
  (molecule index in input order, cell id with 0 = unassigned/noise,
  assignment confidence),
* records command, exit code, binary sha256, version info and git label, and
* computes per-dataset metrics into ``metrics.json`` (sim vs truth, real
  replicate-vs-replicate agreement, optional cellAdmix audit).

With ``--suite NAME`` the dataset selection and run configuration come from
the suite manifest (``benchmarks/datasets/suites.yaml``): every step of the
suite runs in order, each group of steps sharing a run root (see
``suites.py``). ``--dry-run`` resolves and prints the plan without touching
the binary.

Outputs live under ``$BAYSOR_BENCH_DATA/runs/<run_id>/<dataset>/rep<k>/``.
"""
from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common                      # noqa: E402
import metrics as m                # noqa: E402
import celladmix as camix          # noqa: E402

DEFAULT_THREADS = 6

# pristine guard value of metrics.AMI_MAX_LABELS (run_body switches it off
# and back on per invocation for --no-ami)
_AMI_DEFAULT = m.AMI_MAX_LABELS


# ---------------------------------------------------------------------------
# binary probing / command construction
# ---------------------------------------------------------------------------

def probe_binary(baysor: Path) -> dict:
    """Capture ``baysor run --help`` / ``baysor --help`` and detect flags."""
    def _help(args):
        try:
            p = subprocess.run([str(baysor)] + args, capture_output=True, text=True,
                               timeout=60)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise SystemExit(f"failed to probe {baysor} {' '.join(args)}: {exc}")
        return (p.stdout or "") + (p.stderr or "")

    help_run = _help(["run", "--help"])
    help_main = _help(["--help"])
    flags = {name: (f"--{name}" in help_run) for name in (
        "output-style", "prior-segmentation-confidence", "scale-std",
        "skip-ncv-color", "force-2d")}
    version_lines = [ln for ln in help_main.splitlines() if ln.strip()][:6]
    return {
        "path": str(baysor),
        "sha256": common.sha256_file(baysor),
        "version_info": "\n".join(version_lines) or None,
        "help_run": help_run,
        "flags": flags,
        "probed_at": common.utc_now(),
    }


def build_command(baysor: Path, ds: common.Dataset, seg_dir: Path,
                  probe: dict, repo: Path, scale_factor: float = 1.0) -> list[str]:
    """Assemble the ``baysor run`` command line from ``meta.json``.

    ``scale_factor`` multiplies ``baysor.scale_um`` (used to build degraded
    runs for validation; must stay 1.0 for normal benchmark runs).

    ``baysor.extra_args`` is appended verbatim after ``-c <config>`` (so the
    dataset's explicit flags override the config); the builder skips any of
    its own flags that ``extra_args`` already provides, keeping each option
    exactly once as CLI11 requires.
    """
    cfg = ds.baysor_cfg
    extra = list(cfg.get("extra_args") or [])
    extra_flags = {tok.split("=", 1)[0] for tok in extra if tok.startswith("-")}
    cmd: list[str] = [str(baysor), "run", str(ds.molecules_path)]

    def provided(flag_group) -> bool:
        return bool(set(flag_group) & extra_flags)

    def add(flag_group, *values) -> None:
        if provided(flag_group):
            return   # extra_args carries it (verbatim, after -c)
        cmd.extend([next(iter(flag_group)), *map(str, values)])

    prior = cfg.get("prior", "none")
    prior = "none" if prior in (None, "") else str(prior)
    if prior == "column":
        cmd.append(":prior")
    elif prior.startswith("image:"):
        img = ds.path / prior[len("image:"):]
        if not img.is_file():
            raise FileNotFoundError(f"prior image for {ds.id} not found: {img}")
        cmd.append(str(img))
    elif prior != "none":
        raise ValueError(f"{ds.id}: unsupported baysor.prior value {prior!r}")

    add(("-x", "--x-column"), "x")
    add(("-y", "--y-column"), "y")
    if ds.has_z:
        add(("-z", "--z-column"), "z")
    add(("-g", "--gene-column"), "gene")

    if cfg.get("scale_um") is not None:
        scale = float(cfg["scale_um"]) * scale_factor
        add(("-s", "--scale"), scale)
    if cfg.get("scale_std") is not None:
        add(("--scale-std",), cfg["scale_std"])
    if prior != "none" and cfg.get("prior_confidence") is not None:
        if probe["flags"].get("prior-segmentation-confidence", True):
            add(("--prior-segmentation-confidence",), cfg["prior_confidence"])
    if cfg.get("min_molecules_per_cell") is not None:
        add(("-m", "--min-molecules-per-cell"), cfg["min_molecules_per_cell"])
    if cfg.get("config"):
        cfg_path = repo / str(cfg["config"])
        if not cfg_path.is_file():
            raise FileNotFoundError(f"{ds.id}: config {cfg_path} not found")
        add(("-c", "--config"), cfg_path)
    cmd += extra

    if probe["flags"].get("output-style", False):
        add(("--output-style",), "parquet")
    add(("-o", "--output"), seg_dir)
    return cmd


# ---------------------------------------------------------------------------
# process execution with timeout + /usr/bin/time -v
# ---------------------------------------------------------------------------

def _parse_elapsed(text: str) -> Optional[float]:
    """Parse GNU time 'h:mm:ss' / 'm:ss' / 'ss.cc' elapsed strings."""
    parts = text.split(":")
    try:
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        if len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
        return float(parts[0])
    except ValueError:
        return None


def parse_time_v(stderr: str) -> dict:
    """Extract wall time, CPU usage and peak RSS from ``/usr/bin/time -v``
    output (all fields are None when the block is absent)."""
    out: dict = {"wall_s": None, "peak_rss_kb": None, "exit_status": None,
                 "user_s": None, "sys_s": None, "cpu_percent": None}
    mo = re.search(r"Maximum resident set size \(kbytes\): (\d+)", stderr)
    if mo:
        out["peak_rss_kb"] = int(mo.group(1))
    mo = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", stderr)
    if mo:
        out["wall_s"] = _parse_elapsed(mo.group(1))
    mo = re.search(r"Exit status: (\d+)", stderr)
    if mo:
        out["exit_status"] = int(mo.group(1))
    mo = re.search(r"User time \(seconds\): (\S+)", stderr)
    if mo:
        try:
            out["user_s"] = float(mo.group(1))
        except ValueError:
            pass
    mo = re.search(r"System time \(seconds\): (\S+)", stderr)
    if mo:
        try:
            out["sys_s"] = float(mo.group(1))
        except ValueError:
            pass
    mo = re.search(r"Percent of CPU this job got: (\d+)%", stderr)
    if mo:
        out["cpu_percent"] = int(mo.group(1))
    return out


def execute(cmd: list[str], env: dict, timeout: Optional[float]) -> dict:
    """Run ``cmd`` in its own process group, killing it on timeout."""
    start = time.monotonic()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, start_new_session=True, env=env)
    timed_out = False
    try:
        out, err = proc.communicate(timeout=timeout if timeout else None)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        out, err = proc.communicate()
        timed_out = True
    wall = time.monotonic() - start
    return {"returncode": proc.returncode, "stdout": out or "", "stderr": err or "",
            "wall_s": wall, "timed_out": timed_out}


# ---------------------------------------------------------------------------
# one replicate
# ---------------------------------------------------------------------------

def prepare_typing(ds: common.Dataset, rep: int, rep_dir: Path, assign_path: Path,
                   args, repo: Path, root: Path) -> dict:
    """Resolve the cell typing passed to the audit for this replicate.

    Returns ``{"celltypes": Path|None, "fixed_pairs": Path|None,
    "save_celltypes": Path|None, "source": str, "error": str|None}``.

    * ``--celltypes-from BASELINE``: transfer the baseline's saved types
      (clustered once when the baseline was created) onto this replicate's
      segmentation with ``transfer.py`` and score the baseline's fixed pair
      set, so rates are comparable across runs.
    * otherwise: replicate 0 is quick-clustered once (``--save-celltypes``)
      and replicates >= 1 reuse that typing via transfer; the run is marked
      ``quick_cluster``-anchored in metrics.json.
    """
    out = {"celltypes": None, "fixed_pairs": None, "save_celltypes": None,
           "source": "quick_cluster", "error": None}
    if args.celltypes_from:
        out["source"] = "baseline_transferred"
        base_dir = root / "baselines" / args.celltypes_from / ds.id
        tr = camix.run_transfer(
            ds.molecules_path, base_dir / "rep0" / "assignment.parquet",
            base_dir / "celltypes.parquet", assign_path,
            rep_dir / "celltypes.parquet",
            report=rep_dir / "celltypes_transfer.json", repo=repo)
        if tr.get("status") == camix.STATUS_OK:
            out["celltypes"] = rep_dir / "celltypes.parquet"
        else:
            out["error"] = "baseline typing transfer failed: " + str(
                tr.get("reason") or (tr.get("stderr_tail") or "")[-400:])
            return out
        fixed = base_dir / "fixed_pairs.json"
        if fixed.is_file():
            out["fixed_pairs"] = fixed
        return out
    if rep == 0:
        out["save_celltypes"] = rep_dir / "celltypes.parquet"
        return out
    # replicates >= 1: anchor on replicate 0's typing instead of re-clustering
    run_root = root / "runs" / args.run_id / ds.id
    anchor_types = run_root / "rep0" / "celltypes.parquet"
    anchor_assign = run_root / "rep0" / "assignment.parquet"
    if anchor_types.is_file() and anchor_assign.is_file():
        tr = camix.run_transfer(ds.molecules_path, anchor_assign, anchor_types,
                                assign_path, rep_dir / "celltypes.parquet",
                                report=rep_dir / "celltypes_transfer.json",
                                repo=repo)
        if tr.get("status") == camix.STATUS_OK:
            out["celltypes"] = rep_dir / "celltypes.parquet"
            out["source"] = "run_rep0"
        else:
            # fall back to quick clustering for this replicate; the audit's
            # parameters.typing records which path was taken
            pass
    return out


def run_replicate(ds: common.Dataset, rep: int, rep_dir: Path, baysor: Path,
                  probe: dict, args, repo: Path, root: Path) -> dict:
    """Run one replicate; returns the rep record for metrics.json."""
    seg_dir = rep_dir / "seg"
    seg_dir.mkdir(parents=True, exist_ok=True)
    cmd = build_command(baysor, ds, seg_dir, probe, repo,
                        scale_factor=args.scale_factor)
    full_cmd = ["/usr/bin/time", "-v"] + cmd
    env = os.environ.copy()
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        env[var] = str(args.threads)

    print(f"[{ds.id}] rep{rep}: {shlex.join(cmd)}", flush=True)
    res = execute(full_cmd, env, args.timeout)
    tv = parse_time_v(res["stderr"])
    log_path = rep_dir / "baysor.log"
    with open(log_path, "w") as fh:
        fh.write(f"# command: {shlex.join(full_cmd)}\n"
                 f"# exit: {res['returncode']} timed_out: {res['timed_out']}\n"
                 f"--- stdout ---\n{res['stdout']}\n"
                 f"--- stderr ---\n{res['stderr']}\n")
    if res["timed_out"]:
        status = "timeout"
    elif res["returncode"] == 0:
        status = "ok"
    else:
        status = "failed"

    record = {
        "rep": rep,
        "status": status,
        "exit_code": res["returncode"],
        "command": cmd,
        "command_str": shlex.join(cmd),
        "wall_s": round(res["wall_s"], 3),
        "wall_s_time_v": tv["wall_s"],
        "peak_rss_kb": tv["peak_rss_kb"],
        "cpu_user_s": tv["user_s"],
        "cpu_sys_s": tv["sys_s"],
        "cpu_percent": tv["cpu_percent"],
        "threads": args.threads,
        "binary_sha256": probe["sha256"],
        "scale_factor": args.scale_factor,
        "log": "baysor.log",
        "stderr_tail": res["stderr"][-8000:],
    }

    if status == "ok":
        assign_df = common.normalize_assignment(ds.path, seg_dir)
        assign_path = rep_dir / "assignment.parquet"
        assign_df.to_parquet(assign_path, index=False)
        cells = assign_df["cell"].to_numpy(np.int64)
        record["assignment"] = str(assign_path.relative_to(rep_dir.parent))
        record["assignment_sha256"] = common.sha256_file(assign_path)
        record["n_cells"] = int(len(np.unique(cells[cells > 0])))
        record["n_assigned"] = int((cells > 0).sum())
        # molecules without a segmentation row were dropped by Baysor's loader
        record["n_loader_filtered"] = int(assign_df["confidence"].isna().sum())
        seg_file = common._find_segmentation_file(seg_dir)
        record["seg_source"] = seg_file.name
        # the cellAdmix admixture audit applies to real data only
        if ds.kind == "real" and camix.audit_available(repo) \
                and not args.no_celladmix:
            image = None
            images = ds.meta.get("images") or []
            if images:
                cand = ds.path / images[0].get("file", "")
                if cand.is_file():
                    image = cand
            typing = prepare_typing(ds, rep, rep_dir, assign_path, args, repo, root)
            if typing["error"]:
                record["celladmix"] = {"status": camix.STATUS_UNAVAILABLE,
                                       "reason": typing["error"],
                                       "typing_source": typing["source"]}
            else:
                audit = camix.run_audit(
                    ds.molecules_path, assign_path, rep_dir / "celladmix.json",
                    image=image, repo=repo, threads=args.threads,
                    celltypes=typing["celltypes"],
                    fixed_pairs=typing["fixed_pairs"],
                    save_celltypes=typing["save_celltypes"])
                audit["typing_source"] = typing["source"]
                audit["fixed_pairs"] = str(typing["fixed_pairs"]) \
                    if typing["fixed_pairs"] else None
                record["celladmix"] = audit
        elif ds.kind == "real" and args.no_celladmix:
            record["celladmix"] = {"status": camix.STATUS_DISABLED}
        elif ds.kind == "real":
            record["celladmix"] = {"status": camix.STATUS_ABSENT}

    common.write_json(rep_dir / "run.json", record)
    print(f"[{ds.id}] rep{rep}: {status} "
          f"exit={record['exit_code']} wall={record['wall_s']}s "
          f"rss={record['peak_rss_kb']}kB", flush=True)
    return record


def reusable_run(prev: dict, probe: dict, args) -> tuple[bool, str]:
    """Whether an existing run.json may be reused by ``--skip-existing``.

    The check is per replicate: status, binary sha256, thread count and
    scale factor must all match the current invocation, otherwise the
    replicate is rerun.
    """
    if prev.get("status") != "ok":
        return False, f"status={prev.get('status')}"
    if prev.get("binary_sha256") != probe["sha256"]:
        return False, "binary sha256 differs"
    if prev.get("threads") != args.threads:
        return False, f"threads {prev.get('threads')} != {args.threads}"
    if prev.get("scale_factor") != args.scale_factor:
        return False, f"scale factor {prev.get('scale_factor')} != {args.scale_factor}"
    return True, ""


# ---------------------------------------------------------------------------
# metrics aggregation
# ---------------------------------------------------------------------------

def _mean_sd(values: list[float]) -> tuple[float, float]:
    arr = np.asarray([v for v in values if v is not None], dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return float("nan"), float("nan")
    mean = float(arr.mean())
    sd = float(arr.std(ddof=1)) if len(arr) > 1 else 0.0
    return mean, sd


def _aggregate_dicts(dicts: list[dict]) -> dict:
    keys = sorted({k for d in dicts for k in d})
    mean, sd = {}, {}
    for k in keys:
        vals = [d.get(k) for d in dicts]
        mean[k], sd[k] = _mean_sd([v for v in vals])
    return mean, sd


def vendor_labels(molecules: pd.DataFrame) -> Optional[np.ndarray]:
    """Factorize the ``cell_vendor`` column to int64 labels (0 = unassigned)."""
    if "cell_vendor" not in molecules.columns:
        return None
    raw = molecules["cell_vendor"].fillna("").astype(str).str.strip()
    uniq = sorted(v for v in raw.unique() if v != "")
    mapping = {v: i + 1 for i, v in enumerate(uniq)}
    return np.array([mapping.get(v, 0) for v in raw], dtype=np.int64)


def aggregate_dataset(ds: common.Dataset, run_id: str, rep_records: list[dict],
                      binary: dict, threads: int) -> dict:
    """Runtime / provenance part of metrics.json (no quality metrics yet)."""
    molecules = pd.read_parquet(ds.molecules_path)
    ok = [r for r in rep_records if r["status"] == "ok"]
    out: dict = {
        "schema": 1,
        "run_id": run_id,
        "dataset": {
            "id": ds.id, "kind": ds.kind, "tier": ds.tier, "path": str(ds.path),
            "n_molecules": int(len(molecules)),
            "n_genes": int(molecules["gene"].nunique()) if "gene" in molecules else None,
        },
        "inputs": {
            "molecules": str(ds.molecules_path),
            "molecules_sha256": common.sha256_file(ds.molecules_path),
            "meta": str(ds.path / "meta.json"),
            "meta_sha256": common.sha256_file(ds.path / "meta.json"),
        },
        "binary": {k: binary.get(k) for k in ("path", "sha256", "version_info")},
        "label": binary.get("label"),
        "threads": threads,
        "replicates": len(rep_records),
        "created": common.utc_now(),
        "reps": rep_records,
        "runtime": {
            "wall_s": [r.get("wall_s") for r in rep_records],
            "peak_rss_kb": [r.get("peak_rss_kb") for r in rep_records],
        },
        "failures": [r for r in rep_records if r["status"] != "ok"],
    }
    wall_ok = [r.get("wall_s") for r in ok if r.get("wall_s") is not None]
    rss_ok = [r.get("peak_rss_kb") for r in ok if r.get("peak_rss_kb") is not None]
    out["runtime"]["wall_s_mean"], out["runtime"]["wall_s_sd"] = _mean_sd(wall_ok)
    out["runtime"]["peak_rss_kb_mean"], out["runtime"]["peak_rss_kb_sd"] = _mean_sd(rss_ok)
    return out


def compute_dataset_metrics(ds: common.Dataset, run_id: str, rep_dir_for,
                            rep_records: list[dict], binary: dict,
                            threads: int, typing: Optional[dict] = None) -> dict:
    """Full metrics.json: runtime info + sim or real metric blocks."""
    out = aggregate_dataset(ds, run_id, rep_records, binary, threads)
    ok = [r for r in rep_records if r["status"] == "ok"]
    if not ok:
        return out
    molecules = pd.read_parquet(ds.molecules_path)
    cells = [common.assignment_cells(rep_dir_for(r) / "assignment.parquet")
             for r in ok]

    if ds.kind == "sim":
        truth = molecules["cell"].to_numpy(np.int64)
        interior = (molecules["interior"].to_numpy(bool)
                    if "interior" in molecules.columns else None)
        truth_meta = ds.meta.get("truth") or {}
        oracle = truth_meta.get("oracle_accuracy")
        per_rep = [m.sim_metrics(c, truth, interior, oracle) for c in cells]
        mean, sd = _aggregate_dicts(per_rep)
        out["sim"] = {
            "oracle_accuracy": oracle,
            "per_rep": per_rep,
            "mean": mean,
            "sd": sd,
            "n_metric_reps": len(per_rep),
        }
    else:
        pairs, per_pair = [], []
        ok_reps = [(r["rep"], c) for r, c in zip(ok, cells)]
        for i in range(len(ok_reps)):
            for j in range(i + 1, len(ok_reps)):
                pairs.append([ok_reps[i][0], ok_reps[j][0]])
                per_pair.append(m.real_pair_metrics(ok_reps[j][1], ok_reps[i][1]))
        real: dict = {
            "rep_agreement": {},
            "n_metric_reps": len(ok_reps),
        }
        if per_pair:
            mean, sd = _aggregate_dicts(per_pair)
            real["rep_agreement"] = {"pairs": pairs, "per_pair": per_pair,
                                     "mean": mean, "sd": sd}
        vend = vendor_labels(molecules)
        if vend is not None:
            per_rep_v = [m.real_pair_metrics(c, vend) for c in cells]
            vm, vsd = _aggregate_dicts(per_rep_v)
            real["vs_vendor"] = {"per_rep": per_rep_v, "mean": vm, "sd": vsd,
                                 "information_only": True}
        cam = []
        totals: list = []
        capable: list = []
        for r in ok:
            audit = r.get("celladmix") or {}
            entry = {"status": audit.get("status", camix.STATUS_DISABLED)}
            if audit.get("typing_source"):
                entry["typing_source"] = audit["typing_source"]
            if audit.get("status") == camix.STATUS_OK:
                entry.update({k: v for k, v in audit.items()
                              if k not in ("status", "typing_source")})
                totals.append(audit.get("total_admixture_rate"))
                capable.append(audit.get("admixture_capable"))
            else:
                entry["reason"] = audit.get("reason")
                totals.append(None)   # aligned with `ok` run replicates; never 0.0
            cam.append(entry)
        statuses = [c["status"] for c in cam]
        if statuses and all(s == camix.STATUS_OK for s in statuses):
            ad_status = camix.STATUS_OK
        elif any(s == camix.STATUS_FAILED for s in statuses):
            ad_status = camix.STATUS_FAILED
        elif statuses and all(s == statuses[0] for s in statuses):
            ad_status = statuses[0]
        elif any(s == camix.STATUS_OK for s in statuses):
            ad_status = "partial"
        else:
            ad_status = statuses[0] if statuses else camix.STATUS_DISABLED
        totals_f = [float(t) for t in totals if t is not None
                    and not np.isnan(float(t))]
        cap_f = [bool(c) for c in capable if c is not None]
        real["celladmix"] = {
            "per_rep": cam,
            "status": ad_status,
            "typing": typing,
            # True only when every scored replicate has >= 2000 cells;
            # comparisons should gate admixture on this flag
            "admixture_capable": (all(cap_f) if cap_f else None),
            "mean_total": (_mean_sd(totals_f)[0] if totals_f else None),
            "sd_total": (_mean_sd(totals_f)[1] if totals_f else None),
            "per_rep_total": totals,
        }
        out["real"] = real
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _select_datasets(root: Path, args, ap) -> list:
    """Resolve ``--datasets`` (or a suite step's spec) to dataset objects."""
    try:
        selected = common.select_datasets(root, args.datasets, kind=args.kind)
    except ValueError as exc:
        ap.error(str(exc))
    if not selected:
        ap.error(f"no datasets matched {args.datasets!r} under {root}")
    return selected


def _validate_celltypes(selected, args, root, ap) -> None:
    """The ``--celltypes-from`` baseline must carry typing for every selected
    real dataset (setup error, exit 2)."""
    if not (args.celltypes_from and not args.no_celladmix):
        return
    base_root = root / "baselines" / args.celltypes_from
    if not base_root.is_dir():
        ap.error(f"baseline '{args.celltypes_from}' not found at {base_root} "
                 f"(create it with baseline.py create)")
    missing = []
    for ds in selected:
        if ds.kind != "real":
            continue
        if not (base_root / ds.id / "celltypes.parquet").is_file():
            missing.append(f"{ds.id}: {base_root / ds.id / 'celltypes.parquet'}")
        if not (base_root / ds.id / "rep0" / "assignment.parquet").is_file():
            missing.append(f"{ds.id}: rep0 assignment table")
    if missing:
        ap.error("baseline '" + args.celltypes_from + "' lacks typing for:\n  "
                 + "\n  ".join(missing)
                 + "\nrecreate it from a run whose real datasets have audit "
                   "celltypes, or drop --celltypes-from")


def _resolve_baysor(args, ap) -> Path:
    if not args.baysor:
        ap.error("--baysor PATH is required (only --dry-run may omit it)")
    baysor = Path(args.baysor).expanduser().resolve()
    if not baysor.is_file():
        ap.error(f"baysor binary not found: {baysor}")
    if not os.access(baysor, os.X_OK):
        ap.error(f"baysor binary not executable: {baysor}")
    return baysor


def _step_args(args, step, run_id: str):
    """Copy of the CLI namespace configured for one suite step."""
    import suites
    ns = argparse.Namespace(**vars(args))
    ns.datasets = suites.step_spec(step)
    ns.threads = step.threads
    ns.replicates = step.replicates
    ns.timeout = step.timeout
    ns.run_id = run_id
    ns.no_celladmix = (not step.celladmix) or args.no_celladmix
    ns.celltypes_from = args.celltypes_from or step.celltypes_from
    ns.no_ami = args.no_ami or step.no_ami
    return ns


def _print_run_plan(args, selected) -> None:
    print(f"plan: run-id={args.run_id} threads={args.threads} "
          f"replicates={args.replicates} timeout={args.timeout or 'none'}s "
          f"celladmix={'off' if args.no_celladmix else 'on'} "
          f"celltypes-from={args.celltypes_from or '-'} "
          f"ami={'skipped' if args.no_ami else 'computed'}")
    print(f"datasets ({len(selected)}): "
          + ", ".join(d.id for d in selected))


def run_suite(args, ap, repo: Path, root: Path) -> int:
    """Run (or, with ``--dry-run``, resolve and print) a suite's steps."""
    import suites
    for flag, val in (("--threads", args.threads),
                      ("--replicates", args.replicates),
                      ("--timeout", args.timeout),
                      ("--datasets", args.datasets)):
        if val is not None:
            ap.error(f"{flag} cannot be combined with --suite; the suite "
                     f"manifest defines it per step")
    suite = None
    try:
        suite = suites.resolve(args.suite, args.manifest)
    except ValueError as exc:
        ap.error(str(exc))
    steps = suite.steps
    if args.step:
        steps = [s for s in steps if s.name == args.step]
        if not steps:
            ap.error(f"no step named {args.step!r} in suite '{suite.name}' "
                     f"(steps: {', '.join(s.name for s in suite.steps)})")
    group_ids = suites.group_run_ids(suite, args.run_id)
    prepared = []
    for step in steps:
        step_args = _step_args(args, step, group_ids[step.group])
        selected = _select_datasets(root, step_args, ap)
        _validate_celltypes(selected, step_args, root, ap)
        prepared.append((step, step_args, selected))
    if args.dry_run:
        selections = {st.name: [d.id for d in sel]
                      for st, _, sel in prepared}
        print(suites.plan_text(suite, args.run_id, selections,
                               resources_csv=suites.resources_path(suite, root),
                               run_ids=group_ids, root=root))
        return 0
    baysor = _resolve_baysor(args, ap)
    label = args.label or _default_label(repo)
    probe = probe_binary(baysor)
    probe["label"] = label
    any_failure = False
    for step, step_args, selected in prepared:
        run_root = root / "runs" / step_args.run_id
        run_root.mkdir(parents=True, exist_ok=True)
        common.write_json(run_root / "_suite.json", {
            "suite": suite.name, "base_run_id": args.run_id,
            "run_id": step_args.run_id, "group": step.group,
            "steps": [s.name for s in suite.steps if s.group == step.group]})
        print(f"== suite {suite.name} / step {step.name} "
              f"(run-id {step_args.run_id}, threads={step.threads}, "
              f"replicates={step.replicates}, "
              f"celladmix={'off' if step_args.no_celladmix else 'on'}) ==",
              flush=True)
        rc = run_body(step_args, selected, probe, repo, root)
        any_failure = any_failure or rc != 0
    return 1 if any_failure else 0


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baysor", default=None,
                    help="path to the Baysor binary "
                         "(required unless --dry-run)")
    ap.add_argument("--datasets", default=None,
                    help="dataset ids/globs (comma-separated) or tier: quick|full|all "
                         "(with --suite: from the manifest)")
    ap.add_argument("--kind", choices=["sim", "real"],
                    help="restrict to sim or real datasets")
    ap.add_argument("--run-id", required=True, help="name of the run (output folder; "
                    "with --suite the base id, one folder per step group)")
    ap.add_argument("--threads", type=int, default=None,
                    help=f"threads for Baysor (OMP_NUM_THREADS; default "
                         f"{DEFAULT_THREADS}; with --suite: per step)")
    ap.add_argument("--replicates", type=int, default=None,
                    help="number of repeated runs per dataset (default 1; "
                         "with --suite: per step)")
    ap.add_argument("--timeout", type=float, default=None,
                    help="per-replicate timeout in seconds (0 = none; with "
                         "--suite: per step)")
    ap.add_argument("--suite", default=None, metavar="NAME",
                    help="run a suite from the manifest (benchmarks/datasets/"
                         "suites.yaml, e.g. regular|release): dataset ids, "
                         "threads, replicates, timeouts, audit and metrics "
                         "options come from the manifest")
    ap.add_argument("--step", default=None, metavar="STEP",
                    help="with --suite: run only this named step")
    ap.add_argument("--manifest", default=None,
                    help="suite manifest path (default <repo>/benchmarks/"
                         "datasets/suites.yaml)")
    ap.add_argument("--data-root", default=None,
                    help="data root (default $BAYSOR_BENCH_DATA or <repo>/.bench-data)")
    ap.add_argument("--label", default=None,
                    help="git SHA recorded with the run (default: current HEAD)")
    ap.add_argument("--no-celladmix", action="store_true",
                    help="skip the cellAdmix admixture audit")
    ap.add_argument("--no-ami", action="store_true",
                    help="skip the informational AMI metrics (no gate uses "
                         "AMI; saves ~30-60 s per small-cell-count sim "
                         "replicate of metrics time)")
    ap.add_argument("--celltypes-from", default=None, metavar="BASELINE",
                    help="transfer this baseline's saved cell types onto every "
                         "replicate (transfer.py) and audit the baseline's fixed "
                         "pair set; without it the run is quick-cluster-anchored "
                         "and marked as such in metrics.json")
    ap.add_argument("--scale-factor", type=float, default=1.0,
                    help="multiply baysor.scale_um from meta.json (degraded-"
                         "run experiments only; default 1.0)")
    ap.add_argument("--skip-existing", action="store_true",
                    help="skip replicates that already have a successful run.json")
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve and print the run plan (datasets, threads, "
                         "replicates, estimated time) without executing Baysor")
    args = ap.parse_args(argv)

    if args.step and not args.suite:
        ap.error("--step requires --suite")
    if args.replicates is not None and args.replicates < 1:
        ap.error("--replicates must be >= 1")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.run_id):
        ap.error("--run-id may only contain letters, digits, '.', '_', '-'")

    repo = common.repo_root()
    root = common.data_root(args.data_root)
    if args.suite:
        return run_suite(args, ap, repo, root)

    # single-run mode
    if args.threads is None:
        args.threads = DEFAULT_THREADS
    if args.replicates is None:
        args.replicates = 1
    if args.timeout is None:
        args.timeout = 0.0

    selected = _select_datasets(root, args, ap)
    _validate_celltypes(selected, args, root, ap)

    if args.dry_run:
        _print_run_plan(args, selected)
        return 0
    baysor = _resolve_baysor(args, ap)
    label = args.label or _default_label(repo)
    probe = probe_binary(baysor)
    probe["label"] = label
    return run_body(args, selected, probe, repo, root)


def run_body(args, selected: list, probe: dict, repo: Path, root: Path) -> int:
    """Execute one run root: every selected dataset x replicate, then metrics."""
    m.AMI_MAX_LABELS = 0 if args.no_ami else _AMI_DEFAULT
    baysor = Path(probe["path"])
    run_root = root / "runs" / args.run_id
    run_root.mkdir(parents=True, exist_ok=True)
    common.write_json(run_root / "_binary.json", probe)
    record_selection(run_root, selected, args)

    typing = None
    if args.celltypes_from:
        typing = {"mode": "baseline_transferred", "baseline": args.celltypes_from}
    elif any(ds.kind == "real" for ds in selected) and not args.no_celladmix:
        typing = {"mode": "quick_cluster", "baseline": None}

    any_failure = False
    for ds in selected:
        ds_root = run_root / ds.id
        rep_records = []
        for k in range(args.replicates):
            rep_dir = ds_root / f"rep{k}"
            existing = rep_dir / "run.json"
            if args.skip_existing and existing.is_file():
                prev = common.read_json(existing)
                ok_reuse, why = reusable_run(prev, probe, args)
                if ok_reuse:
                    print(f"[{ds.id}] rep{k}: reusing existing run", flush=True)
                    rep_records.append(prev)
                    continue
                print(f"[{ds.id}] rep{k}: not reusing existing run ({why}); "
                      f"rerunning", flush=True)
            if rep_dir.exists():
                shutil.rmtree(rep_dir)   # never mix artifacts of two binaries
            rep_dir.mkdir(parents=True, exist_ok=True)
            rec = run_replicate(ds, k, rep_dir, baysor, probe, args, repo, root)
            rep_records.append(rec)
        if any(r["status"] != "ok" for r in rep_records):
            any_failure = True

        def rep_dir_for(rec, _ds_root=ds_root):
            return _ds_root / f"rep{rec['rep']}"

        mjson = compute_dataset_metrics(ds, args.run_id, rep_dir_for, rep_records,
                                        probe, args.threads, typing=typing)
        if args.no_ami:
            mjson["metric_options"] = {"ami": "skipped (--no-ami; "
                                                 "informational only)"}
        common.write_json(ds_root / "metrics.json", mjson)
        kind = ds.kind
        if kind == "sim" and mjson.get("sim"):
            acc = mjson["sim"]["mean"].get("matched_accuracy")
            print(f"[{ds.id}] matched_accuracy mean={acc:.4f} "
                  f"sd={mjson['sim']['sd'].get('matched_accuracy'):.4f} "
                  f"({mjson['sim']['n_metric_reps']} reps)", flush=True)
        elif kind == "real" and mjson.get("real", {}).get("rep_agreement"):
            ra = mjson["real"]["rep_agreement"]["mean"]
            print(f"[{ds.id}] rep agreement: ari={ra.get('molecule_ari'):.4f} "
                  f"assigned={ra.get('assigned_agreement'):.4f}", flush=True)

    print(f"run '{args.run_id}': "
          f"{'OK' if not any_failure else 'FINISHED WITH FAILURES'} -> {run_root}")
    return 0 if not any_failure else 1


def record_selection(run_root: Path, selected: list, args) -> dict:
    """Merge this invocation's dataset selection into ``_selection.json``.

    ``compare.py`` reads it so that a deliberate subset run does not fail on
    baseline datasets outside its selection (they are reported as skipped),
    while a dataset the run *intended* to execute but that has no
    ``metrics.json`` still fails the comparison.
    """
    sel_path = run_root / "_selection.json"
    try:
        sel = common.read_json(sel_path) if sel_path.is_file() else {}
    except (OSError, ValueError):
        sel = {}
    if not isinstance(sel, dict):
        sel = {}
    sel["datasets"] = sorted(set(sel.get("datasets") or [])
                             | {d.id for d in selected})
    inv = sel.get("invocations") or []
    if not isinstance(inv, list):
        inv = []
    inv.append({"spec": args.datasets, "kind": args.kind,
                "replicates": args.replicates, "threads": args.threads,
                "scale_factor": args.scale_factor, "at": common.utc_now()})
    sel["invocations"] = inv
    common.write_json(sel_path, sel)
    return sel


def _default_label(repo: Path) -> str:
    try:
        sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=30).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"],
                               capture_output=True, text=True, timeout=30).stdout
        return sha + ("+dirty" if dirty.strip() else "")
    except OSError:
        return "unknown"


if __name__ == "__main__":
    raise SystemExit(main())
