#!/usr/bin/env python3
"""Generate ``SUMMARY.md`` for a committed baseline.

Reads the committed per-dataset metric JSONs of
``benchmarks/baselines/<name>/`` plus each dataset's ``meta.json``
(oracle/naive accuracies) and renders:

* **sim** — one-to-one accuracy (mean +/- SD), ARI over assigned molecules,
  recovery, cell-count ratio, oracle / naive accuracy and the oracle gap;
* **real** — replicate-vs-replicate ARI over assigned molecules, cell count,
  assigned fraction, agreement with the vendor segmentation, and the
  cellAdmix ``total_admixture_rate`` (Baysor vs the vendor segmentation
  when a ``vendor_audit.json`` exists next to the baseline's data files);
* **runtime** — wall time and peak RSS per dataset;
* **highlights** — where Baysor is far below the oracle or the naive
  baseline, prior vs no-prior pairs, and performance outliers.

Usage::

    baseline_summary.py --baseline bugfixes-35e8a7e [--out PATH]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common                      # noqa: E402

ORACLE_GAP_FLAG = 0.05      # accuracy this far below the oracle is "far below"
SLOW_PER_MOL = 0.6          # ms per molecule flagged as a performance outlier
RSS_FLAG_GB = 8.0           # peak RSS flagged as an outlier


def _msd(mean: Optional[float], sd: Optional[float], digits: int = 4) -> str:
    if mean is None or (isinstance(mean, float) and math.isnan(mean)):
        return "—"
    if sd is None or (isinstance(sd, float) and math.isnan(sd)) or sd == 0:
        return f"{mean:.{digits}f}"
    return f"{mean:.{digits}f} ± {sd:.{digits}f}"


def _num(v: Optional[float], digits: int = 4) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "—"
    return f"{v:.{digits}f}"


def load_baseline(baseline_dir: Path) -> dict[str, dict]:
    out = {}
    for p in sorted(baseline_dir.glob("*.json")):
        m = common.read_json(p)
        ds_id = (m.get("dataset") or {}).get("id") or p.stem
        out[ds_id] = m
    return out


def _meta_for(ds_id: str, metrics: dict, data_root: Path) -> dict:
    path = Path((metrics.get("dataset") or {}).get("path")
                or (data_root / "sim" / ds_id / "meta.json"))
    meta_path = path / "meta.json" if path.is_dir() else None
    if meta_path and meta_path.is_file():
        return common.read_json(meta_path)
    for kind in ("sim", "real"):
        cand = data_root / kind / ds_id / "meta.json"
        if cand.is_file():
            return common.read_json(cand)
    return {}


def _vendor_audit_rate(ds_id: str, data_root: Path, name: str) -> Optional[float]:
    p = data_root / "baselines" / name / ds_id / "vendor_audit.json"
    if not p.is_file():
        return None
    try:
        metrics = common.read_json(p).get("metrics") or {}
    except (OSError, ValueError):
        return None
    if metrics.get("status") != "ok":
        return None
    return metrics.get("total_admixture_rate")


def _rep_cell_counts(m: dict) -> list[int]:
    return [r.get("n_cells") for r in m.get("reps", [])
            if r.get("status") == "ok" and r.get("n_cells") is not None]


def _assigned_fraction(m: dict) -> Optional[float]:
    n_as = [r.get("n_assigned") for r in m.get("reps", [])
            if r.get("status") == "ok" and r.get("n_assigned") is not None]
    if not n_as:
        return None
    n_mol = (m.get("dataset") or {}).get("n_molecules")
    if not n_mol:
        return None
    return sum(n_as) / len(n_as) / n_mol


# ---------------------------------------------------------------------------
# sections
# ---------------------------------------------------------------------------

def sim_table(metrics_by_id: dict[str, dict], data_root: Path) -> tuple[str, list[dict]]:
    lines = [
        "| dataset | tier | accuracy_1to1 | ari_assigned | recovery | "
        "cell_count_ratio | oracle | naive | gap (oracle−acc) | wall s | peak RSS GB |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    rows = []
    for ds_id in sorted(metrics_by_id):
        m = metrics_by_id[ds_id]
        sim = m.get("sim") or {}
        mean = sim.get("mean") or {}
        sd = sim.get("sd") or {}
        meta = _meta_for(ds_id, m, data_root)
        truth = meta.get("truth") or {}
        oracle = sim.get("oracle_accuracy")
        naive = truth.get("naive_accuracy")
        acc = mean.get("accuracy_1to1")
        gap = None
        if oracle is not None and acc is not None:
            gap = oracle - acc
        rt = m.get("runtime") or {}
        rss_gb = (rt.get("peak_rss_kb_mean") or 0) / 1024 / 1024
        lines.append("| " + " | ".join([
            f"`{ds_id}`", str((m.get("dataset") or {}).get("tier") or "?"),
            _msd(acc, sd.get("accuracy_1to1")),
            _msd(mean.get("ari_assigned"), sd.get("ari_assigned")),
            _msd(mean.get("recovery_rate"), sd.get("recovery_rate")),
            _msd(mean.get("cell_count_ratio"), sd.get("cell_count_ratio")),
            _num(oracle), _num(naive), _num(gap),
            _msd(rt.get("wall_s_mean"), rt.get("wall_s_sd"), 1),
            _num(rss_gb, 2),
        ]) + " |")
        rows.append({"id": ds_id, "tier": (m.get("dataset") or {}).get("tier"),
                     "acc": acc, "oracle": oracle, "naive": naive, "gap": gap,
                     "wall": rt.get("wall_s_mean"),
                     "rss_gb": rss_gb,
                     "n_molecules": (m.get("dataset") or {}).get("n_molecules"),
                     "prior": ((meta.get("baysor") or {}).get("prior") or "none")})
    return "\n".join(lines), rows


def real_table(metrics_by_id: dict[str, dict], data_root: Path,
               name: str) -> tuple[str, list[dict]]:
    lines = [
        "| dataset | tier | rep ARI (assigned) | cells | assigned frac | "
        "vs vendor ARI | vs vendor cells matched | admixture (Baysor) | "
        "admixture (vendor) | wall s | peak RSS GB |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    rows = []
    for ds_id in sorted(metrics_by_id):
        m = metrics_by_id[ds_id]
        real = m.get("real") or {}
        rep = real.get("rep_agreement") or {}
        rmean, rsd = rep.get("mean") or {}, rep.get("sd") or {}
        vend = (real.get("vs_vendor") or {}).get("mean") or {}
        cam = real.get("celladmix") or {}
        cells = _rep_cell_counts(m)
        cells_s = (f"{sum(cells) / len(cells):,.0f}" if cells else "—")
        v_rate = _vendor_audit_rate(ds_id, data_root, name)
        rt = m.get("runtime") or {}
        rss_gb = (rt.get("peak_rss_kb_mean") or 0) / 1024 / 1024
        lines.append("| " + " | ".join([
            f"`{ds_id}`", str((m.get("dataset") or {}).get("tier") or "?"),
            _msd(rmean.get("ari_assigned"), rsd.get("ari_assigned")),
            cells_s, _num(_assigned_fraction(m)),
            _num(vend.get("ari_assigned")),
            _num(vend.get("frac_cells_matched")),
            _msd(cam.get("mean_total"), cam.get("sd_total"), 5),
            _msd(v_rate, None, 5),
            _msd(rt.get("wall_s_mean"), rt.get("wall_s_sd"), 1),
            _num(rss_gb, 2),
        ]) + " |")
        rows.append({"id": ds_id, "tier": (m.get("dataset") or {}).get("tier"),
                     "ari_rep": rmean.get("ari_assigned"),
                     "cells": sum(cells) / len(cells) if cells else None,
                     "assigned_frac": _assigned_fraction(m),
                     "vs_vendor_ari": vend.get("ari_assigned"),
                     "admx": cam.get("mean_total"), "admx_vendor": v_rate,
                     "wall": rt.get("wall_s_mean"), "rss_gb": rss_gb,
                     "n_molecules": (m.get("dataset") or {}).get("n_molecules")})
    return "\n".join(lines), rows


def runtime_table(metrics_by_id: dict[str, dict]) -> str:
    lines = ["| dataset | kind | tier | molecules | wall s (mean ± SD) | "
             "peak RSS GB |",
             "|---|---|---|---|---|---|"]
    for ds_id in sorted(metrics_by_id):
        m = metrics_by_id[ds_id]
        rt = m.get("runtime") or {}
        lines.append("| " + " | ".join([
            f"`{ds_id}`", str((m.get("dataset") or {}).get("kind") or "?"),
            str((m.get("dataset") or {}).get("tier") or "?"),
            f"{(m.get('dataset') or {}).get('n_molecules') or 0:,}",
            _msd(rt.get("wall_s_mean"), rt.get("wall_s_sd"), 1),
            _num((rt.get("peak_rss_kb_mean") or 0) / 1024 / 1024, 2),
        ]) + " |")
    return "\n".join(lines)


def highlights(sim_rows: list[dict], real_rows: list[dict],
               all_metrics: dict[str, dict]) -> list[str]:
    out: list[str] = []

    # --- far below oracle / naive ----------------------------------------
    below = [r for r in sim_rows
             if r["gap"] is not None and r["gap"] >= ORACLE_GAP_FLAG]
    below.sort(key=lambda r: r["gap"], reverse=True)
    under_naive = [r for r in sim_rows
                   if r["naive"] is not None and r["acc"] is not None
                   and r["acc"] < r["naive"]]
    out.append("### Baysor far below the oracle or the naive baseline")
    out.append("")
    if below:
        out.append(f"* below the oracle by >= {ORACLE_GAP_FLAG:.2f} "
                   "(accuracy_1to1 vs `meta.truth.oracle_accuracy`):")
        for r in below:
            out.append(f"  * `{r['id']}` — acc {_num(r['acc'])}, oracle "
                       f"{_num(r['oracle'])}, gap **{r['gap']:+.3f}**")
    else:
        out.append("* no sim dataset is more than "
                   f"{ORACLE_GAP_FLAG:.2f} below its oracle")
    if under_naive:
        out.append("* below the nearest-nucleus naive baseline:")
        for r in sorted(under_naive, key=lambda r: r["acc"]):
            out.append(f"  * `{r['id']}` — acc {_num(r['acc'])} < naive "
                       f"{_num(r['naive'])}")
    out.append("")

    # --- prior vs no-prior ------------------------------------------------
    out.append("### Prior vs no-prior")
    out.append("")
    by_id = {r["id"]: r for r in sim_rows}
    pairs = []
    for r in sim_rows:
        if r["id"].endswith("_noprior"):
            base = by_id.get(r["id"][: -len("_noprior")])
            if base:
                pairs.append((base, r))
    if pairs:
        out.append("| base (prior) | accuracy | no-prior variant | accuracy "
                   "| Δ (prior − noprior) |")
        out.append("|---|---|---|---|---|")
        for base, nop in sorted(pairs, key=lambda p: p[0]["id"]):
            d = (base["acc"] or 0) - (nop["acc"] or 0)
            out.append(f"| `{base['id']}` | {_num(base['acc'])} | "
                       f"`{nop['id']}` | {_num(nop['acc'])} | {d:+.4f} |")
        deltas = [((b["acc"] or 0) - (n["acc"] or 0)) for b, n in pairs]
        out.append("")
        out.append(f"Prior helps on average by "
                   f"**{sum(deltas) / len(deltas):+.4f}** over {len(pairs)} "
                   "matched pairs (identical molecules, truth and parameters).")
    imprior = [r for r in sim_rows if r["id"].endswith("_imprior")]
    if imprior:
        out.append("")
        out.append("Imperfect-prior variants: "
                   + ", ".join(f"`{r['id']}` acc {_num(r['acc'])}"
                               for r in sorted(imprior, key=lambda r: r["id"])))
    out.append("")

    # --- performance outliers ---------------------------------------------
    out.append("### Performance outliers")
    out.append("")
    timed = [r for r in sim_rows + real_rows if r["wall"]]
    if timed:
        per_mol = [r["wall"] / max(r["n_molecules"] or 1, 1) * 1000
                   for r in timed]          # ms per molecule
        med = sorted(per_mol)[len(per_mol) // 2]
        slow = [r for r, ppm in zip(timed, per_mol)
                if ppm > max(SLOW_PER_MOL, 4 * med)]
        slow.sort(key=lambda r: r["wall"], reverse=True)
        if slow:
            out.append(f"* slow per molecule (> {max(SLOW_PER_MOL, 4 * med):.2f} "
                       "ms/mol, > 4× the median):")
            for r in slow:
                ppm = r["wall"] / max(r["n_molecules"] or 1, 1) * 1000
                out.append(f"  * `{r['id']}` — {r['wall']:.0f} s for "
                           f"{r['n_molecules']:,} molecules "
                           f"({ppm:.2f} ms/mol)")
        big = sorted((r for r in timed if (r["rss_gb"] or 0) >= RSS_FLAG_GB),
                     key=lambda r: -(r["rss_gb"] or 0))
        if big:
            out.append(f"* peak RSS >= {RSS_FLAG_GB:g} GB:")
            for r in big:
                out.append(f"  * `{r['id']}` — {r['rss_gb']:.1f} GB")
        out.append(f"* median costs {med:.2f} ms/molecule across all "
                   f"{len(timed)} datasets.")
    out.append("")

    # --- failures -----------------------------------------------------------
    failures = []
    for ds_id, m in sorted(all_metrics.items()):
        for rec in m.get("failures") or []:
            failures.append((ds_id, rec.get("rep"), rec.get("status"),
                             rec.get("wall_s")))
    out.append("### Failures and timeouts")
    out.append("")
    if failures:
        for ds_id, rep, status, wall in failures:
            out.append(f"* `{ds_id}` rep{rep}: **{status}**"
                       + (f" after {wall:.0f} s" if wall else ""))
    else:
        out.append("* none — every replicate of every dataset succeeded.")
    out.append("")
    return out


# ---------------------------------------------------------------------------
# rendering / CLI
# ---------------------------------------------------------------------------

def render(name: str, metrics_by_id: dict[str, dict], data_root: Path) -> str:
    sim = {k: v for k, v in metrics_by_id.items()
           if (v.get("dataset") or {}).get("kind") == "sim"}
    real = {k: v for k, v in metrics_by_id.items()
            if (v.get("dataset") or {}).get("kind") == "real"}

    first = next(iter(metrics_by_id.values()), {})
    binary = first.get("binary") or {}
    some = list(metrics_by_id.values())
    threads = {m.get("threads") for m in some}
    reps = {m.get("replicates") for m in some}
    bl = (first.get("baseline") or {})
    created = bl.get("created")

    sim_md, sim_rows = sim_table(sim, data_root)
    real_md, real_rows = real_table(real, data_root, name)

    lines = [
        f"# Baseline summary: `{name}`",
        "",
        "Generated by [`harness/baseline_summary.py`](../../harness/"
        "baseline_summary.py) from the committed baseline JSONs:",
        "",
        "```bash",
        f".deps/bench/bin/python benchmarks/harness/baseline_summary.py "
        f"--baseline {name}",
        "```",
        "",
        f"* binary sha256: `{binary.get('sha256')}`",
        f"* label: `{first.get('label')}`",
        f"* threads: {sorted(threads)}; replicates: {sorted(reps)}; "
        f"created: {created}",
        f"* datasets: {len(sim)} sim, {len(real)} real",
        "",
        "Accuracy metrics are mean ± sample SD over replicates "
        "(sim) or over replicate pairs (real). `vs vendor` columns compare "
        "each replicate with `cell_vendor` (informational). `admixture` is "
        "the cellAdmix `total_admixture_rate` (lower = cleaner); the vendor "
        "column is the vendor segmentation scored on the same typing and "
        "fixed pair set (`vendor_audit.json` in the baseline's data "
        "directory).",
        "",
    ]
    if sim:
        lines += ["## Simulated datasets", "", sim_md, ""]
    if real:
        lines += ["## Real datasets", "", real_md, ""]
    lines += ["## Runtime and peak RSS", "", runtime_table(metrics_by_id), ""]
    lines += ["## Highlights", ""]
    lines += highlights(sim_rows, real_rows, metrics_by_id)
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--repo", default=None)
    ap.add_argument("--out", default=None,
                    help="default <repo>/benchmarks/baselines/<name>/SUMMARY.md")
    args = ap.parse_args(argv)

    repo = Path(args.repo).resolve() if args.repo else common.repo_root()
    root = common.data_root(args.data_root)
    baseline_dir = repo / "benchmarks" / "baselines" / args.baseline
    if not baseline_dir.is_dir():
        print(f"error: no such baseline: {baseline_dir}", file=sys.stderr)
        return 2
    metrics_by_id = load_baseline(baseline_dir)
    if not metrics_by_id:
        print(f"error: no baseline JSONs in {baseline_dir}", file=sys.stderr)
        return 2
    text = render(args.baseline, metrics_by_id, root) + "\n"
    out = Path(args.out) if args.out else baseline_dir / "SUMMARY.md"
    out.write_text(text)
    print(f"wrote {out} ({len(metrics_by_id)} datasets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
