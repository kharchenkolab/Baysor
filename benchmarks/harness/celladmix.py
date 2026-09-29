"""Optional integration with the cellAdmix admixture audit (BENCH-CELLADMIX).

The audit wrapper lives in ``benchmarks/celladmix/audit.py`` and is invoked as::

    python benchmarks/celladmix/audit.py --molecules <parquet> \
        --assignment <parquet> --out <json> [--image ...] [--threads N] \
        [--celltypes <parquet>] [--fixed-pairs <json>]

Audit JSON contract (what this adapter relies on):

* ``metrics.total_admixture_rate`` — ``None`` unless ``metrics.status ==
  "ok"`` and ``n_pairs_evaluated > 0``; a rate of 0.0 is never fabricated;
* ``metrics.n_pairs_evaluated`` / ``metrics.n_pairs_detected`` /
  ``metrics.admixture_capable`` (crop has >= 2000 cells);
* ``counts.n_cells``;
* ``pairs_top`` — top detected (source, target) pairs;
* ``parameters.typing`` — ``quick_cluster`` or ``celltypes``.

:func:`normalize_audit` maps that layout onto the harness's flat rep-record
block (top-level ``total_admixture_rate``, ``pairs``, ...). Everything here
degrades gracefully: when the module is missing, disabled or crashing,
callers get ``{"status": ...}`` instead of raising.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

STATUS_OK = "ok"
STATUS_ABSENT = "absent"
STATUS_DISABLED = "disabled"
STATUS_FAILED = "failed"
STATUS_UNAVAILABLE = "unavailable"


def audit_script_path(repo: Optional[Path] = None) -> Path:
    if repo is None:
        repo = Path(__file__).resolve().parents[2]
    return repo / "benchmarks" / "celladmix" / "audit.py"


def transfer_script_path(repo: Optional[Path] = None) -> Path:
    if repo is None:
        repo = Path(__file__).resolve().parents[2]
    return repo / "benchmarks" / "celladmix" / "transfer.py"


def audit_available(repo: Optional[Path] = None) -> bool:
    return audit_script_path(repo).is_file()


def format_command(python: str, script: Path, molecules: Path, assignment: Path,
                   out_json: Path, image: Optional[Path] = None,
                   threads: Optional[int] = None,
                   celltypes: Optional[Path] = None,
                   fixed_pairs: Optional[Path] = None,
                   save_celltypes: Optional[Path] = None,
                   extra_args: Optional[list[str]] = None) -> list[str]:
    cmd = [python, str(script), "--molecules", str(molecules),
           "--assignment", str(assignment), "--out", str(out_json)]
    if image is not None:
        cmd += ["--image", str(image)]
    if threads is not None:
        cmd += ["--threads", str(threads)]
    if celltypes is not None:
        cmd += ["--celltypes", str(celltypes)]
    if fixed_pairs is not None:
        cmd += ["--fixed-pairs", str(fixed_pairs)]
    if save_celltypes is not None:
        cmd += ["--save-celltypes", str(save_celltypes)]
    if extra_args:
        cmd += list(extra_args)
    return cmd


def normalize_audit(result: dict) -> dict:
    """Flatten a raw audit JSON onto the harness contract.

    Raises ``ValueError`` when the JSON does not follow the contract at all
    (callers turn that into ``{"status": "failed"}``). An audit whose status
    is not ``ok``, that evaluated no pairs, or whose rate is null becomes
    ``status: unavailable`` with ``total_admixture_rate: None`` — never 0.0.
    """
    if not isinstance(result, dict) or "metrics" not in result \
            or "counts" not in result:
        raise ValueError("audit JSON lacks the metrics/counts blocks "
                         "(is benchmarks/celladmix/audit.py current?)")
    metrics = result.get("metrics") or {}
    counts = result.get("counts") or {}
    params = result.get("parameters") or {}
    audit_status = metrics.get("status")
    rate = metrics.get("total_admixture_rate")
    evaluated = metrics.get("n_pairs_evaluated") or 0

    out = {
        "audit_status": audit_status,
        "total_admixture_rate": None,
        "total_admixture_molecules": metrics.get("total_admixture_molecules"),
        "pairs": result.get("pairs_top") or [],
        "n_pairs_evaluated": int(evaluated),
        "n_pairs_detected": metrics.get("n_pairs_detected"),
        "n_cells": counts.get("n_cells"),
        "admixture_capable": metrics.get("admixture_capable"),
        "typing": params.get("typing"),
        "seed": params.get("seed"),
        "audit_threads": params.get("threads"),
    }
    if audit_status != "ok" or evaluated <= 0 or rate is None:
        reason = metrics.get("reason")
        if audit_status not in (None, "ok"):
            reason = f"audit status {audit_status!r}" + \
                     (f": {reason}" if reason else "")
        elif evaluated <= 0:
            reason = "no pairs evaluated"
        else:
            reason = "audit returned a null rate"
        out["reason"] = reason
        out["status"] = STATUS_UNAVAILABLE
        return out
    out["total_admixture_rate"] = float(rate)
    out["status"] = STATUS_OK
    return out


def run_audit(molecules: Path, assignment: Path, out_json: Path,
              image: Optional[Path] = None, repo: Optional[Path] = None,
              timeout: float = 1800.0, threads: Optional[int] = None,
              celltypes: Optional[Path] = None,
              fixed_pairs: Optional[Path] = None,
              save_celltypes: Optional[Path] = None,
              extra_args: Optional[list[str]] = None) -> dict:
    """Run the audit wrapper and return its normalized JSON (plus ``status``).

    Never raises for audit problems: returns
    ``{"status": "absent" | "failed" | "unavailable", ...}`` instead.
    """
    script = audit_script_path(repo)
    if not script.is_file():
        return {"status": STATUS_ABSENT,
                "reason": f"{script} not found (BENCH-CELLADMIX not landed yet)"}
    cmd = format_command(sys.executable, script, molecules, assignment, out_json,
                         image=image, threads=threads, celltypes=celltypes,
                         fixed_pairs=fixed_pairs, save_celltypes=save_celltypes,
                         extra_args=extra_args)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": STATUS_FAILED, "reason": f"timeout after {timeout}s"}
    except OSError as exc:
        return {"status": STATUS_FAILED, "reason": str(exc)}
    if proc.returncode != 0 or not out_json.is_file():
        return {"status": STATUS_FAILED, "returncode": proc.returncode,
                "stderr_tail": proc.stderr[-2000:]}
    try:
        with open(out_json) as f:
            result = json.load(f)
    except (OSError, json.JSONDecodeError) as exc:
        return {"status": STATUS_FAILED, "reason": f"bad audit output: {exc}"}
    try:
        normalized = normalize_audit(result)
    except ValueError as exc:
        return {"status": STATUS_FAILED, "reason": str(exc), "raw": result}
    normalized["audit_json"] = str(out_json)
    return normalized


def run_transfer(molecules: Path, baseline_assignment: Path,
                 baseline_celltypes: Path, target_assignment: Path,
                 out_celltypes: Path, report: Optional[Path] = None,
                 repo: Optional[Path] = None, timeout: float = 600.0) -> dict:
    """Transfer baseline cell types onto a target segmentation (transfer.py).

    Returns the transfer report with a ``status`` key; never raises for
    transfer problems.
    """
    script = transfer_script_path(repo)
    if not script.is_file():
        return {"status": STATUS_ABSENT,
                "reason": f"{script} not found (BENCH-CELLADMIX not landed yet)"}
    cmd = [sys.executable, str(script),
           "--molecules", str(molecules),
           "--baseline-assignment", str(baseline_assignment),
           "--baseline-celltypes", str(baseline_celltypes),
           "--target-assignment", str(target_assignment),
           "--out", str(out_celltypes)]
    if report is not None:
        cmd += ["--report", str(report)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"status": STATUS_FAILED, "reason": f"timeout after {timeout}s"}
    except OSError as exc:
        return {"status": STATUS_FAILED, "reason": str(exc)}
    if proc.returncode != 0 or not out_celltypes.is_file():
        return {"status": STATUS_FAILED, "returncode": proc.returncode,
                "stderr_tail": proc.stderr[-2000:]}
    stats: dict = {}
    if report is not None and report.is_file():
        try:
            with open(report) as f:
                stats = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            return {"status": STATUS_FAILED, "reason": f"bad transfer report: {exc}"}
    stats["status"] = STATUS_OK
    return stats
