#!/usr/bin/env python3
"""Generate ``benchmarks/DATASETS.md`` — the inventory of every dataset.

The inventory is a table of all datasets under ``$BAYSOR_BENCH_DATA``
(``sim/`` and ``real/``), grouped by kind and then by platform (real) or
generator (sim), plus a density x gene-panel coverage matrix per kind.

Columns (see the dataset contract in ``benchmarks/README.md``):

* id, kind, tier, platform/generator, tissue or scenario;
* genes, molecules, area (mm2), cells/mm2;
* density class, gene-panel class;
* 2D/3D, prior type, images;
* ``admixture_capable`` — real datasets with >= 2000 vendor cells (the
  cellAdmix audit's power threshold); the per-run audit flag (Baysor's own
  cell count >= 2000) is recorded in ``metrics.json``;
* source (URL/DOI or generator) and notes (perf-stress, variant tags, ...);
* resource usage (6-thread CPU/wall/peak RAM, 1-thread wall/RAM, cellAdmix
  audit time) from the committed ``baselines/bugfixes-35e8a7e/resources.csv``
  (``TODO`` = never measured) and the ``Suite`` membership column from
  ``datasets/suites.yaml``.

Usage::

    inventory.py                     # (re)write benchmarks/DATASETS.md
    inventory.py --out /tmp/x.md     # write elsewhere
    inventory.py --check             # exit 1 if the file is stale
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

import pyarrow.parquet as pq
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common                      # noqa: E402
import resources as resmod         # noqa: E402
import suites as suites_mod        # noqa: E402
from validate_datasets import (     # noqa: E402
    cell_density_class, gene_panel_class, load_manifest_index)

ADMIXTURE_MIN_CELLS = 2000

# committed resource table + suite manifest feeding the extra columns
DEFAULT_RESOURCES = "benchmarks/baselines/bugfixes-35e8a7e/resources.csv"

PANEL_COLS = ("tiny", "small", "medium", "large", "huge")
DENSITY_ROWS = ("sparse", "medium", "dense")


# ---------------------------------------------------------------------------
# row extraction
# ---------------------------------------------------------------------------

def _scenario(ds: common.Dataset, entry: dict) -> str:
    """Scenario string for a sim dataset (manifest scenario, else a compact
    st-recoverability parameter summary)."""
    if entry and entry.get("scenario"):
        return str(entry["scenario"])
    truth = ds.meta.get("truth") or {}
    params = truth.get("params") or {}
    if "packing_cells_per_mm2" in params:
        model = params.get("model", "?")
        sigma = params.get("sigma_um", "?")
        packing = params.get("packing_cells_per_mm2", "?")
        geom = params.get("geometry", "voronoi")
        out = f"strec {model}, sigma={sigma}, p={packing}"
        if geom and geom != "voronoi":
            out += f", {geom}"
        return out
    src_scenario = (ds.meta.get("source") or {}).get("scenario")
    if src_scenario:
        return str(src_scenario)
    if entry and entry.get("base"):
        return f"derived from {entry['base']}"
    return truth.get("scenario") or ""


def _tissue(ds: common.Dataset, entry: dict) -> str:
    if entry and entry.get("tissue"):
        return str(entry["tissue"])
    src = ds.meta.get("source") or {}
    return str(src.get("original_dataset") or "")[:60]


def _short_source(ds: common.Dataset, entry: dict) -> str:
    src = ds.meta.get("source") or {}
    if ds.kind == "sim":
        gen = Path(str(src.get("generator", "?"))).name
        seed = src.get("seed")
        return f"{gen}, seed {seed}" if seed is not None else gen
    doi = src.get("doi")
    if doi:
        return f"doi:{doi}"
    url = str(src.get("url") or "")
    parsed = urlparse(url)
    head = (parsed.netloc + parsed.path).rstrip("/")
    if len(head) > 44:
        head = head[:41] + "..."
    return head or url


def _notes(ds: common.Dataset) -> str:
    tags: list[str] = []
    notes = str((ds.meta.get("difficulty") or {}).get("notes") or "")
    if "perf-stress" in notes:
        tags.append("perf-stress: default ICA/mrf init hangs on the ~19k-gene "
                    "crop; uses --cluster-method louvain")
    extra = [str(t) for t in (ds.meta.get("baysor") or {}).get("extra_args") or []]
    if "--cluster-method" in extra:
        i = extra.index("--cluster-method")
        if i + 1 < len(extra) and "louvain" not in " ".join(tags):
            tags.append(f"cluster-method {extra[i + 1]} (extra_args)")
    ds_id = ds.id
    if ds_id.endswith("_noprior"):
        tags.append("no-prior variant of "
                    f"{ds_id[: -len('_noprior')]}")
    if ds_id.endswith("_imprior"):
        tags.append("imperfect-prior variant of "
                    f"{ds_id[: -len('_imprior')]}")
    if ds_id.endswith("_admix"):
        tags.append("cellAdmix audit crop (>= 2000 vendor cells)")
    if ds_id.endswith("_3d"):
        tags.append("3D")
    if ds_id.endswith("_full"):
        tags.append("full tier")
    return "; ".join(tags)


def dataset_row(ds: common.Dataset, entry: dict, z_present: bool,
                res_row: Optional[dict] = None,
                suite_of: Optional[dict] = None) -> dict:
    stats = ds.meta["stats"]
    diff = ds.meta.get("difficulty") or {}
    area = float(stats["area_um2"])
    cells_per_mm2 = ((stats.get("n_true_cells") if ds.kind == "sim"
                      else stats.get("n_vendor_cells")) or 0) / area * 1e6
    images = ds.meta.get("images") or []
    prior = str((ds.meta.get("baysor") or {}).get("prior") or "none").lower()
    prior = "column" if prior == "column" else (
        "image" if prior.startswith("image:") else "none")
    if ds.kind == "sim":
        admix = "—"
    else:
        admix = "yes" if stats.get("n_vendor_cells", 0) >= ADMIXTURE_MIN_CELLS \
            else "no"
    platform = str(ds.meta.get("platform") or "")
    generator = Path(str((ds.meta.get("source") or {}).get("generator", ""))).name
    r = res_row or {}
    if ds.kind == "sim":
        audit = "—"          # the cellAdmix audit never applies to sim data
    else:
        audit = resmod.fmt_duration(r.get("audit_wall_s"))
    suites_here = []
    if suite_of:
        for name in ("regular", "release"):
            if ds.id in (suite_of.get(name) or set()):
                suites_here.append(name)
    return {
        "id": ds.id,
        "kind": ds.kind,
        "tier": ds.tier or "?",
        "platform": platform if ds.kind == "real" else generator,
        "tissue": _tissue(ds, entry) if ds.kind == "real"
        else _scenario(ds, entry),
        "genes": int(stats["n_genes"]),
        "molecules": int(stats["n_molecules"]),
        "area_mm2": area / 1e6,
        "cells_per_mm2": cells_per_mm2,
        "density": diff.get("cell_density", "?"),
        "panel": diff.get("gene_panel", "?"),
        "dim": "3D" if z_present else "2D",
        "prior": prior,
        "images": ", ".join(str(im.get("name") or Path(im["file"]).stem)
                            for im in images) if images else "—",
        "admixture_capable": admix,
        "cpu6": resmod.fmt_mean_sd(r.get("cpu6_mean_s"), r.get("cpu6_sd_s")),
        "wall6": resmod.fmt_duration(r.get("wall6_mean_s")),
        "ram6": resmod.fmt_bytes_kb(r.get("peak_rss6_kb")),
        "wallram1": resmod.fmt_wall_ram(r.get("wall1_s"), r.get("peak_rss1_kb")),
        "audit": audit,
        "suite": "+".join(suites_here) if suites_here else "—",
        "source": _short_source(ds, entry),
        "notes": _notes(ds),
        "group": (platform if ds.kind == "real" else generator),
    }


def collect_rows(root: Path, repo: Path,
                 resources_csv: Optional[Path] = None) -> list[dict]:
    manifests = load_manifest_index(repo)
    table = resmod.load_csv(resources_csv) if resources_csv else {}
    try:
        suite_of = suites_mod.membership(root, repo)
    except ValueError as exc:
        print(f"warning: suites not resolved ({exc}); Suite column is '—'",
              file=sys.stderr)
        suite_of = {}
    rows: list[dict] = []
    for ds in common.discover_datasets(root):
        entry = (manifests.get(ds.id) or {}).get("entry") or {}
        z_present = "z" in pq.read_schema(ds.molecules_path).names
        rows.append(dataset_row(ds, entry, z_present,
                                table.get(ds.id), suite_of))
    return rows


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------

_HEADER = ["id", "kind", "tier", "platform/generator", "tissue/scenario",
           "genes", "molecules", "area mm²", "cells/mm²", "density",
           "panel", "dim", "prior", "images", "admixture_capable",
           "CPU time (6 thr)", "wall (6 thr)", "peak RAM (6 thr)",
           "wall/RAM (1 thr)", "audit time",
           "source", "notes", "Suite"]


def _fmt_row(r: dict) -> str:
    return "| " + " | ".join([
        f"`{r['id']}`", r["kind"], r["tier"], r["platform"], r["tissue"],
        f"{r['genes']:,}", f"{r['molecules']:,}", f"{r['area_mm2']:.4g}",
        f"{r['cells_per_mm2']:,.0f}", r["density"], r["panel"], r["dim"],
        r["prior"], r["images"], r["admixture_capable"],
        r["cpu6"], r["wall6"], r["ram6"], r["wallram1"], r["audit"],
        r["source"], r["notes"] or "—", r["suite"],
    ]) + " |"


def _table(rows: list[dict]) -> list[str]:
    out = ["| " + " | ".join(_HEADER) + " |",
           "|" + "|".join(["---"] * len(_HEADER)) + "|"]
    out += [_fmt_row(r) for r in sorted(rows, key=lambda r: r["id"])]
    return out


def _matrix(rows: list[dict]) -> list[str]:
    counts: dict[tuple[str, str], int] = Counter(
        (r["density"], r["panel"]) for r in rows)
    out = ["| density \\ panel | " + " | ".join(PANEL_COLS) + " | total |",
           "|" + "|".join(["---"] * (len(PANEL_COLS) + 1)) + "|"]
    totals = Counter(r["panel"] for r in rows)
    for d in DENSITY_ROWS:
        cells = [str(counts.get((d, p), 0)) for p in PANEL_COLS]
        out.append(f"| **{d}** | " + " | ".join(cells)
                   + f" | {sum(counts.get((d, p), 0) for p in PANEL_COLS)} |")
    out.append("| **total** | "
               + " | ".join(str(totals.get(p, 0)) for p in PANEL_COLS)
               + f" | {len(rows)} |")
    return out


def render(rows: list[dict], root: Path) -> str:
    sim = [r for r in rows if r["kind"] == "sim"]
    real = [r for r in rows if r["kind"] == "real"]
    n_quick = sum(r["tier"] == "quick" for r in rows)
    n_full = sum(r["tier"] == "full" for r in rows)
    admix = [r["id"] for r in real if r["admixture_capable"] == "yes"]

    lines: list[str] = []
    lines += [
        "# Baysor benchmark dataset inventory",
        "",
        "Generated by [`harness/inventory.py`](harness/inventory.py) from the",
        "dataset contract in [`README.md`](README.md); regenerate with",
        "",
        "```bash",
        "BAYSOR_BENCH_DATA=... .deps/bench/bin/python benchmarks/harness/inventory.py",
        "```",
        "",
        f"{len(rows)} datasets ({n_quick} quick, {n_full} full): "
        f"{len(sim)} simulated, {len(real)} real.",
        f"`admixture_capable` = real dataset with >= {ADMIXTURE_MIN_CELLS} vendor",
        "cells (the cellAdmix audit power threshold; the per-run audit flag uses",
        "Baysor's own cell count and lives in `metrics.json`).",
        "",
        "The resource columns (`CPU time (6 thr)`, `wall (6 thr)`, `peak RAM",
        "(6 thr)`, `wall/RAM (1 thr)`, `audit time`) come from the committed",
        "[`baselines/bugfixes-35e8a7e/resources.csv`](baselines/bugfixes-35e8a7e/resources.csv)",
        "generated by [`harness/resources.py`](harness/resources.py) from the",
        "existing `benchbase-b` / `benchbase-t1` runs (no Baysor rerun).",
        "`TODO` marks a value that was never measured and is never guessed;",
        "`—` marks fields that do not apply (cellAdmix audit on simulated data,",
        "suite membership). `Suite` is membership in",
        "[`datasets/suites.yaml`](datasets/suites.yaml) (`regular+release` = in",
        "the fast every-change suite, `release` = full pre-release suite only).",
        "",
    ]

    groups: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        groups["sim" if r["kind"] == "sim" else "real"].append(r)

    # --- simulated, grouped by generator ----------------------------------
    sim_generators = sorted({r["group"] for r in sim})
    lines.append("## Simulated datasets (`$BAYSOR_BENCH_DATA/sim/`)")
    lines.append("")
    for gen in sim_generators:
        subset = [r for r in sim if r["group"] == gen]
        lines.append(f"### generator: `{gen}` ({len(subset)})")
        lines.append("")
        lines += _table(subset)
        lines.append("")

    # --- real, grouped by platform ----------------------------------------
    platforms = sorted({r["group"] for r in real})
    lines.append("## Real datasets (`$BAYSOR_BENCH_DATA/real/`)")
    lines.append("")
    for plat in platforms:
        subset = [r for r in real if r["group"] == plat]
        lines.append(f"### platform: {plat} ({len(subset)})")
        lines.append("")
        lines += _table(subset)
        lines.append("")

    # --- coverage matrices -------------------------------------------------
    lines.append("## Coverage matrix (density class × gene-panel class)")
    lines.append("")
    for title, subset in (("Simulated", sim), ("Real", real)):
        lines.append(f"### {title}")
        lines.append("")
        lines += _matrix(subset)
        lines.append("")

    lines += [
        "## Notes",
        "",
        f"* `admixture_capable` = yes on {len(admix)} real datasets: "
        + ", ".join(f"`{i}`" for i in sorted(admix)),
        "* density / panel classes follow the contract thresholds in "
        "[`README.md`](README.md); `area mm²` is `stats.area_um2` (the crop "
        "bbox area), `cells/mm²` is the vendor count for real datasets and the "
        "true-cell count for simulated ones.",
        "* resource columns: CPU time (user+sys) and wall are means over the "
        "6-thread replicates (CPU with sample SD), peak RAM is the max over "
        "replicates, wall/RAM (1 thr) is the single 1-thread replicate, audit "
        "time is the cellAdmix `runtime_seconds.total` mean over replicates "
        "(6-thread runs); `TODO` = never measured — never guessed or "
        "back-filled.",
        "* `Suite` from [`datasets/suites.yaml`](datasets/suites.yaml): "
        "`regular` = the <= ~20 min suite run on every change/PR, `release` = "
        "the full pre-release suite (every dataset, 3 replicates + 1-thread "
        "quick tier).",
        "* `source` shows the DOI when the manifest records one, otherwise the "
        "downstream URL; simulated datasets show the generator and its seed.",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--repo", default=None)
    ap.add_argument("--resources", default=None,
                    help="resources CSV for the resource columns "
                         f"(default <repo>/{DEFAULT_RESOURCES})")
    ap.add_argument("--out", default=None,
                    help="output path (default <repo>/benchmarks/DATASETS.md)")
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if the file is out of date")
    args = ap.parse_args(argv)

    repo = Path(args.repo).resolve() if args.repo else common.repo_root()
    root = common.data_root(args.data_root)
    out = Path(args.out) if args.out else repo / "benchmarks" / "DATASETS.md"
    resources_csv = Path(args.resources) if args.resources \
        else repo / DEFAULT_RESOURCES

    rows = collect_rows(root, repo, resources_csv)
    if not rows:
        print(f"error: no datasets under {root}", file=sys.stderr)
        return 2
    text = render(rows, root) + "\n"

    if args.check:
        current = out.read_text() if out.is_file() else ""
        if current != text:
            print(f"{out} is out of date (rerun inventory.py without --check)",
                  file=sys.stderr)
            return 1
        print(f"{out} is up to date ({len(rows)} datasets)")
        return 0

    out.write_text(text)
    print(f"wrote {out} ({len(rows)} datasets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
