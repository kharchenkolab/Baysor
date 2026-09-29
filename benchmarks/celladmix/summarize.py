#!/usr/bin/env python
"""Summarize cellAdmix audit validation runs into a small JSON + Markdown report.

Reads the per-variant audit JSONs produced by ``run_validation.sh``, tabulates
the headline metrics, evaluates the monotonicity checks the benchmark metric
must satisfy (worse segmentations must score higher), and reports the
seed-to-seed tolerance measured on the unmodified input.

Usage:
    python summarize.py --variant vendor=vendor_seed1.json \
        --variant border10=border10.json ... \
        --out summary.json --markdown summary.md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


def row(name: str, data: dict) -> dict:
    m = data["metrics"]
    c = data["counts"]
    return {
        "variant": name,
        "total_admixture_rate": m["total_admixture_rate"],
        "total_admixture_molecules": m["total_admixture_molecules"],
        "n_pairs_detected": m["n_pairs_detected"],
        "n_pairs_evaluated": m["n_pairs_evaluated"],
        "n_cells": c["n_cells"],
        "n_molecules_used": c["n_molecules_used"],
        "n_cell_types": c["n_cell_types"],
        "typing": data["parameters"]["typing"],
        "seed": data["parameters"]["seed"],
        "runtime_seconds": data["runtime_seconds"]["total"],
        "status": m["status"],
    }


def monotone(results: dict[str, dict], chain: list[str]) -> dict:
    """Check a strict increasing chain of total_admixture_rate."""
    missing = [n for n in chain if n not in results]
    if missing:
        return {"chain": chain, "checked": False, "missing": missing}
    values = [results[n]["total_admixture_rate"] for n in chain]
    steps = [{"from": a, "to": b, "rate_from": va, "rate_to": vb, "increased": vb > va}
             for a, b, va, vb in zip(chain, chain[1:], values, values[1:])]
    return {
        "chain": chain,
        "checked": True,
        "values": dict(zip(chain, values)),
        "steps": steps,
        "passed": all(s["increased"] for s in steps),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variant", action="append", required=True, metavar="NAME=PATH",
                        help="named audit result JSON (repeatable)")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)

    named: dict[str, Path] = {}
    for spec in args.variant:
        name, _, path = spec.partition("=")
        if not path:
            raise SystemExit(f"--variant must be NAME=PATH, got {spec!r}")
        named[name] = Path(path)

    results = {name: row(name, load(path)) for name, path in named.items()}

    checks: dict[str, dict] = {}
    checks["border_degradation_monotone"] = monotone(results, ["vendor", "border10", "border30"])
    checks["dilation_worse_than_vendor"] = monotone(results, ["vendor", "dilate2"])
    if "vendor" in results and "vendor_seed2" in results:
        a = results["vendor"]["total_admixture_rate"]
        b = results["vendor_seed2"]["total_admixture_rate"]
        checks["seed_stochasticity"] = {
            "rate_seed1": a,
            "rate_seed2": b,
            "abs_difference": abs(a - b),
            "relative_difference": abs(a - b) / max(a, 1e-12),
            "note": "two full runs with different --seed (the seed reaches the NMF "
                    "fit; the audit is factorization-independent) plus a fresh "
                    "quick-clustering draw; this sets the comparison tolerance",
        }
    if "vendor" in results and "vendor_repeat" in results:
        checks["same_seed_deterministic"] = {
            "abs_difference": abs(results["vendor"]["total_admixture_rate"]
                                  - results["vendor_repeat"]["total_admixture_rate"]),
        }
    if "border30" in results and "border30_recluster" in results:
        checks["transfer_vs_recluster_typing"] = {
            "transfer_rate": results["border30"]["total_admixture_rate"],
            "recluster_rate": results["border30_recluster"]["total_admixture_rate"],
            "vendor_rate": results.get("vendor", {}).get("total_admixture_rate"),
        }

    summary = {"results": results, "checks": checks}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    if args.markdown:
        lines = [
            "# cellAdmix admixture audit — validation on the pancreas Xenium crop",
            "",
            "Lower total admixture = cleaner segmentation.",
            "",
            "| variant | typing | seed | total_admixture_rate | admixed molecules | detected/evaluated pairs | cells | molecules used | runtime (s) |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for name, r in results.items():
            lines.append(
                f"| {name} | {r['typing']} | {r['seed']} | {r['total_admixture_rate']:.6f} "
                f"| {r['total_admixture_molecules']:.0f} "
                f"| {r['n_pairs_detected']}/{r['n_pairs_evaluated']} "
                f"| {r['n_cells']} | {r['n_molecules_used']} | {r['runtime_seconds']:.1f} |")
        lines += ["", "## Checks", ""]
        for name, chk in checks.items():
            status = chk.get("passed")
            if status is None:
                status = "n/a"
            lines.append(f"- **{name}**: {'PASS' if status is True else ('FAIL' if status is False else 'info')} — "
                         f"`{json.dumps(chk)}`")
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text("\n".join(lines) + "\n")

    print(json.dumps({"out": str(args.out),
                      "rates": {k: v["total_admixture_rate"] for k, v in results.items()},
                      "checks_passed": {k: v.get("passed") for k, v in checks.items()}},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
