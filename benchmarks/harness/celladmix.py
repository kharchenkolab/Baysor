"""Optional integration with the cellAdmix admixture audit (BENCH-CELLADMIX).

The audit wrapper lives in ``benchmarks/celladmix/audit.py`` and is expected to
expose::

    python benchmarks/celladmix/audit.py --molecules <parquet> \
        --assignment <parquet> --out <json> [--image ...]

returning a JSON object with ``total_admixture_rate`` and per-pair rates.
Everything here degrades gracefully: when the module is missing, disabled or
crashing, callers get ``{"status": ...}`` instead of raising.
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


def audit_script_path(repo: Optional[Path] = None) -> Path:
    if repo is None:
        repo = Path(__file__).resolve().parents[2]
    return repo / "benchmarks" / "celladmix" / "audit.py"


def audit_available(repo: Optional[Path] = None) -> bool:
    return audit_script_path(repo).is_file()


def format_command(python: str, script: Path, molecules: Path, assignment: Path,
                   out_json: Path, image: Optional[Path] = None) -> list[str]:
    cmd = [python, str(script), "--molecules", str(molecules),
           "--assignment", str(assignment), "--out", str(out_json)]
    if image is not None:
        cmd += ["--image", str(image)]
    return cmd


def run_audit(molecules: Path, assignment: Path, out_json: Path,
              image: Optional[Path] = None, repo: Optional[Path] = None,
              timeout: float = 1800.0) -> dict:
    """Run the audit wrapper and return its JSON (plus a ``status`` key).

    Never raises for audit problems: returns
    ``{"status": "absent" | "failed", ...}`` instead.
    """
    script = audit_script_path(repo)
    if not script.is_file():
        return {"status": STATUS_ABSENT,
                "reason": f"{script} not found (BENCH-CELLADMIX not landed yet)"}
    cmd = format_command(sys.executable, script, molecules, assignment, out_json, image)
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
    if not isinstance(result, dict) or "total_admixture_rate" not in result:
        return {"status": STATUS_FAILED, "reason": "audit JSON lacks total_admixture_rate",
                "raw": result}
    result = dict(result)
    result["status"] = STATUS_OK
    return result
