#!/usr/bin/env python3
"""Fetch and crop public 10x Xenium datasets for the Baysor benchmark suite.

Subcommands (all read ``benchmarks/datasets/real_xenium.yaml`` by default):

``verify``
    HEAD-check every source zip URL and compare its size with the manifest.
``fetch``
    Download the manifest-listed members of each source zip into the cache
    with :mod:`remotezip` (never the whole zip).
``pick``
    Choose crop boxes (``crop.bbox_um``) for entries that do not have one yet,
    from the molecule/cell density histograms of the source bundle.  With
    ``--write`` the chosen boxes are saved into the manifest.
``build``
    Build ``$BAYSOR_BENCH_DATA/real/<id>/`` per the dataset contract:
    ``molecules.parquet``, ``meta.json``, ``reference/`` boundaries, cropped
    focus images where applicable, and a provenance ``README.md``.
``report``
    Print the dataset inventory table (reads ``meta.json`` files and, when
    present, the timing of the sanity Baysor runs).

Example::

    python benchmarks/fetch/xenium.py fetch
    python benchmarks/fetch/xenium.py pick --write
    python benchmarks/fetch/xenium.py build
    python benchmarks/fetch/xenium.py report
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import xenium_common as xc  # noqa: E402


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _datasets(manifest: dict, ids: list[str]) -> list[dict]:
    if not ids:
        return manifest["datasets"]
    known = {d["id"] for d in manifest["datasets"]}
    missing = [i for i in ids if i not in known]
    if missing:
        raise SystemExit(f"unknown dataset ids: {missing}")
    return [d for d in manifest["datasets"] if d["id"] in ids]


def _group_by_source(datasets: list[dict]) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {}
    for ds in datasets:
        groups.setdefault(ds["source"]["url"], []).append(ds)
    return groups


def _source_members(ds: dict, include_images: bool = True) -> list[str]:
    members = list(ds["members"].values())
    if include_images and ds.get("images"):
        members += list(ds["images"]["focus_members"])
    return members


def _source_with_members(ds: dict) -> dict:
    """The dataset's source dict plus its member paths (for the helpers)."""
    return {**ds["source"], "members": ds["members"]}


def _hist_cache(source: dict, cells, min_qv: float, bin_um: float):
    """Load or build the (molecule, cell) density histograms of a source."""
    zip_name = Path(source["url"]).name[:-4]
    x = cells["x_centroid"].to_numpy(dtype=np.float64)
    y = cells["y_centroid"].to_numpy(dtype=np.float64)
    bounds = [
        math.floor((x.min() - 50) / bin_um) * bin_um,
        math.floor((y.min() - 50) / bin_um) * bin_um,
        math.ceil((x.max() + 50) / bin_um) * bin_um,
        math.ceil((y.max() + 50) / bin_um) * bin_um,
    ]
    meta = {"bin_um": bin_um, "min_qv": min_qv, "bounds": bounds}
    hist_path = xc.cache_root() / "hist" / f"{zip_name}.npz"
    if hist_path.exists():
        with np.load(hist_path, allow_pickle=False) as z:
            stored = json.loads(str(z["meta"]))
            if stored == meta:
                return z["mol"], z["cell"], bounds, bin_um
        print("  histogram cache stale, recomputing")
    hist_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"  scanning transcripts for density histogram ({zip_name}) ...", flush=True)
    tp = xc.member_cache_path(source["url"], source["members"]["transcripts"])
    mol = xc.transcript_histogram(tp, min_qv, bounds, bin_um)
    cel = xc.cells_histogram(cells, bounds, bin_um)
    np.savez_compressed(
        hist_path, mol=mol, cell=cel, meta=np.array(json.dumps(meta))
    )
    return mol, cel, bounds, bin_um


def _load_cells(source: dict):
    import pyarrow.parquet as pq

    cp = xc.member_cache_path(source["url"], source["members"]["cells"])
    cols = ["cell_id", "x_centroid", "y_centroid", "cell_area", "nucleus_area"]

    pf = pq.ParquetFile(str(cp))
    use = [c for c in cols if c in pf.schema_arrow.names]
    return pq.read_table(str(cp), columns=use).to_pandas()


# ---------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------


def cmd_verify(args, manifest) -> int:
    ok = True
    for src in xc.unique_sources(manifest):
        status, length = xc.head_info(src["url"])
        expected = src.get("zip_bytes")
        verdict = "OK"
        if status != 200:
            verdict, ok = "FAIL (not reachable)", False
        elif expected and length != expected:
            verdict = f"SIZE MISMATCH (manifest {expected})"
            ok = False
        print(f"[{status}] {length or '?':>13} bytes  {verdict:24}  {src['url']}")
    return 0 if ok else 1


def cmd_fetch(args, manifest) -> int:
    datasets = _datasets(manifest, args.ids)
    total = 0
    for url, group in _group_by_source(datasets).items():
        print(f"source {url}", flush=True)
        members = _source_members(group[0])
        for ds in group[1:]:
            members += _source_members(ds)
        paths = xc.fetch_members(url, members)
        total += sum(p.stat().st_size for p in set(paths.values()))
    print(f"fetched/verified {total / 1e6:.1f} MB (cumulative cache: "
          f"{xc.cached_bytes() / 1e9:.2f} GB)")
    return 0


def cmd_pick(args, manifest) -> int:
    datasets = _datasets(manifest, args.ids)
    bbox_by_id = {
        d["id"]: d["crop"].get("bbox_um") for d in manifest["datasets"]
    }
    bin_um = args.bin_um
    changed = False
    for url, group in _group_by_source(datasets).items():
        needs_hist = [d for d in group if not d["crop"].get("bbox_um")]
        # source members must exist for hist + pick
        any_ds = group[0]
        xc.fetch_members(url, _source_members(any_ds, include_images=False))
        cells = _load_cells(_source_with_members(any_ds))
        if not needs_hist:
            continue
        mol, cel, bounds, bin_um = _hist_cache(
            _source_with_members(any_ds), cells, manifest["defaults"]["min_qv"], bin_um
        )
        for ds in needs_hist:
            crop = ds["crop"]
            within = None
            if crop.get("within"):
                ref = crop["within"]
                within = bbox_by_id.get(ref)
                if within is None:
                    raise SystemExit(
                        f"{ds['id']}: within={ref} has no bbox yet; list it first"
                    )
            res = xc.pick_bbox(
                mol, cel, bin_um, bounds,
                target=crop["target_molecules"],
                min_mols=crop["min_molecules"],
                max_mols=crop["max_molecules"],
                min_cells=crop.get("min_cells", 150),
                min_coverage=crop.get("min_coverage", 0.7),
                density_hint=crop.get("density_hint", "any"),
                within=within,
                seed=args.seed,
            )
            print(
                f"{ds['id']}: bbox={res['bbox_um']} mols={res['mols']} "
                f"cells={res['cells']} cov={res['coverage']:.2f} "
                f"dens={res['cells_per_mm2']:.0f}/mm2 hint={res['density_hint']}"
                + (" RELAXED" if res["relaxed_density"] else "")
            )
            crop["bbox_um"] = res["bbox_um"]
            crop["pick"] = {
                k: res[k]
                for k in ("mols", "cells", "coverage", "cells_per_mm2",
                          "density_hint", "relaxed_density")
            }
            bbox_by_id[ds["id"]] = res["bbox_um"]
            changed = True
    if changed and args.write:
        xc.save_manifest(manifest, args.manifest)
        print(f"manifest updated: {args.manifest}")
    elif changed:
        print("dry run (use --write to save bboxes into the manifest)")
    return 0


def cmd_build(args, manifest) -> int:
    datasets = _datasets(manifest, args.ids)
    defaults = manifest["defaults"]
    min_qv = defaults["min_qv"]
    built = []
    for url, group in _group_by_source(datasets).items():
        any_ds = group[0]
        print(f"source {url}", flush=True)
        xc.fetch_members(url, _source_members(any_ds))
        for ds in group:
            if not ds["crop"].get("bbox_um"):
                raise SystemExit(
                    f"{ds['id']}: crop.bbox_um is null; run `pick --write` first"
                )
        cells = _load_cells(_source_with_members(any_ds))
        tp = xc.member_cache_path(
            url, _source_with_members(any_ds)["members"]["transcripts"]
        )
        bboxes = [ds["crop"]["bbox_um"] for ds in group]
        print(f"  reading {len(bboxes)} crop(s) from transcripts ...", flush=True)
        crops = xc.read_transcript_crops(tp, bboxes, min_qv)
        for ds, crop_tbl in zip(group, crops):
            built.append(_build_one(ds, crop_tbl, cells, defaults))
    print(f"built {len(built)} dataset(s):")
    for line in built:
        print("  " + line)
    return 0


def _build_one(ds: dict, crop_tbl, cells, defaults: dict) -> str:
    import pyarrow.parquet as pq

    src = ds["source"]
    bbox = ds["crop"]["bbox_um"]
    ds_dir = xc.data_root() / "real" / ds["id"]
    ds_dir.mkdir(parents=True, exist_ok=True)

    mol = xc.build_molecule_table(crop_tbl, bbox)
    n_vendor = len(xc.cells_in_bbox(cells, bbox))
    stats = xc.molecule_stats(mol, bbox, n_vendor)
    scale_um, scale_method = xc.estimate_scale_um(cells, bbox)

    gp = xc.member_cache_path(src["url"], ds["members"]["gene_panel"])
    panel_genes = xc.panel_gene_count(gp)

    # reference boundaries for the crop
    ref_dir = ds_dir / "reference"
    ref_dir.mkdir(exist_ok=True)
    for key in ("cell_boundaries", "nucleus_boundaries"):
        if key in ds["members"]:
            bp = xc.member_cache_path(src["url"], ds["members"][key])
            cropped = xc.crop_boundaries(bp, bbox)
            pq.write_table(cropped, ref_dir / f"{key}.parquet", compression="zstd")

    # focus images (if the manifest asks for them)
    images = []
    if ds.get("images"):
        images = _build_images(ds, bbox)

    meta = xc.build_meta(
        dataset_id=ds["id"],
        tier=ds["tier"],
        source=src,
        bbox=bbox,
        note=ds["crop"]["note"],
        stats=stats,
        panel_genes=panel_genes,
        scale_um=scale_um,
        scale_method=scale_method,
        baysor=defaults["baysor"],
        images=images,
        retrieved=src.get("retrieved", manifest_retrieved),
        difficulty_notes=ds.get("difficulty_notes", ds["crop"]["note"]),
    )
    with open(ds_dir / "meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")

    pq.write_table(mol, ds_dir / "molecules.parquet", compression="zstd")

    members_desc = ", ".join(
        f"{k}={v}" for k, v in ds["members"].items()
    )
    extra = [
        "Bundle members used:",
        "",
        f"* {members_desc}",
        f"* source zip: {src['url']}",
        "",
    ]
    (ds_dir / "README.md").write_text(xc.dataset_readme(meta, extra))

    return (
        f"{ds['id']}: {stats['n_molecules']} mol, {stats['n_genes']} genes, "
        f"{stats['n_vendor_cells']} vendor cells, "
        f"{stats['vendor_cells_per_mm2']:.0f}/mm2, scale_um={scale_um}"
    )


def _build_images(ds: dict, bbox: Sequence) -> list[dict]:
    src = ds["source"]
    spec = ds["images"]
    paths = [
        xc.member_cache_path(src["url"], m) for m in spec["focus_members"]
    ]
    channels = spec["channels"]
    # Multi-file OME: one member per channel; any member resolves all channels
    # through its OME XML (sibling files must all be in the cache, they are).
    info = xc.ome_info(paths[0])
    px = info["pixel_size_um"]
    if px is None:
        raise SystemExit(f"{ds['id']}: focus image has no PhysicalSizeX")
    _, height, width = info["shape"]
    x0, y0, x1, y1 = bbox
    c0 = max(0, int(math.floor(x0 / px)))
    r0 = max(0, int(math.floor(y0 / px)))
    c1 = min(width, int(math.ceil(x1 / px)))
    r1 = min(height, int(math.ceil(y1 / px)))
    if (c1 - c0) * px < (x1 - x0) * 0.99 or (r1 - r0) * px < (y1 - y0) * 0.99:
        raise SystemExit(f"{ds['id']}: crop box exceeds the focus image extent")

    out = []
    for name, ch in channels.items():
        arr = xc.read_image_window(paths[0], ch, bbox, px)
        xc.write_tif(xc.data_root() / "real" / ds["id"] / "images" / f"{name}.tif", arr)
        ch_name = info["channels"][ch] if ch < len(info["channels"]) else ch
        out.append({
            "name": name,
            "file": f"images/{name}.tif",
            "pixel_size_um": px,
            "origin_um": [round(c0 * px, 4), round(r0 * px, 4)],
            "source": f"OME channel {ch} ({ch_name}) of "
                      + ", ".join(Path(m).name for m in spec["focus_members"]),
        })
        print(f"  image {name}: channel {ch} ({ch_name}) {arr.shape} "
              f"mean={float(arr.mean()):.1f}")
    return out


manifest_retrieved = ""  # set in main() from the manifest


def cmd_report(args, manifest) -> int:
    rows = []
    for ds in manifest["datasets"]:
        meta_path = xc.data_root() / "real" / ds["id"] / "meta.json"
        if not meta_path.exists():
            rows.append((ds, None, None))
            continue
        with open(meta_path) as fh:
            meta = json.load(fh)
        timing = None
        tpath = xc.data_root() / "runs" / args.run_id / ds["id"] / "timing.json"
        if tpath.exists():
            with open(tpath) as fh:
                timing = json.load(fh)
        rows.append((ds, meta, timing))

    src_bytes: dict[str, int] = {}
    for src in xc.unique_sources(manifest):
        zip_name = Path(src["url"]).name[:-4]
        d = xc.cache_root() / zip_name
        src_bytes[src["url"]] = sum(
            p.stat().st_size for p in d.rglob("*") if p.is_file()
        ) if d.exists() else 0

    header = (
        "| id | tissue | genes | crop size (µm) | molecules | vendor cells | "
        "cells/mm² | density | tier | source URL | bytes | Baysor quick |"
    )
    sep = "|" + "---|" * 12
    print(header)
    print(sep)
    for ds, meta, timing in rows:
        if meta is None:
            print(f"| {ds['id']} | {ds.get('tissue','')} | NOT BUILT | | | | | | "
                  f"{ds['tier']} | {ds['source']['url']} | | |")
            continue
        b = meta["stats"]
        x0, y0, x1, y1 = meta["crop"]["bbox_um"]
        size = f"{x1 - x0:.0f} × {y1 - y0:.0f}"
        t = f"{timing['wall_s']:.1f} s" if timing else "—"
        print(
            f"| {meta['id']} | {ds.get('tissue','')} | {b['n_genes']} | {size} | "
            f"{b['n_molecules']} | {b['n_vendor_cells']} | "
            f"{b['vendor_cells_per_mm2']:.0f} | {meta['difficulty']['cell_density']} | "
            f"{meta['tier']} | {meta['source']['url']} | "
            f"{src_bytes[meta['source']['url']]} | {t} |"
        )
    print()
    print(f"cache total: {xc.cached_bytes() / 1e9:.2f} GB under {xc.cache_root()}")
    return 0


# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    global manifest_retrieved
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-m", "--manifest", default=str(xc.DEFAULT_MANIFEST))
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("verify", help="HEAD-check source URLs")

    p = sub.add_parser("fetch", help="download manifest members into the cache")
    p.add_argument("ids", nargs="*", help="dataset ids (default: all)")

    p = sub.add_parser("pick", help="choose crop bboxes")
    p.add_argument("ids", nargs="*", help="dataset ids (default: all)")
    p.add_argument("--write", action="store_true", help="save bboxes into the manifest")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--bin-um", type=float, default=25.0)

    p = sub.add_parser("build", help="build dataset directories")
    p.add_argument("ids", nargs="*", help="dataset ids (default: all)")

    p = sub.add_parser("report", help="dataset inventory table")
    p.add_argument("--run-id", default="sanity_realx")

    args = ap.parse_args(argv)
    manifest = xc.load_manifest(args.manifest)
    manifest_retrieved = manifest.get("retrieved", "")

    return {
        "verify": cmd_verify,
        "fetch": cmd_fetch,
        "pick": cmd_pick,
        "build": cmd_build,
        "report": cmd_report,
    }[args.cmd](args, manifest)


if __name__ == "__main__":
    raise SystemExit(main())
