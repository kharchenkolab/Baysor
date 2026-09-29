#!/usr/bin/env python3
"""Fetch and crop public 10x Xenium datasets for the Baysor benchmark suite.

Subcommands (all read ``benchmarks/datasets/real_xenium.yaml`` by default):

``verify``
    HEAD-check every source zip URL and compare its size with the manifest,
    then verify the recorded sha256 of every cached member (``downloaded_members``)
    and of every built dataset's ``molecules.parquet``/``meta.json`` (``outputs``).
``fetch``
    Download the manifest-listed members of each source zip into the cache
    with :mod:`remotezip` (never the whole zip) through the shared
    :mod:`download` helper; size + sha256 of every member are recorded under
    ``downloaded_members`` and verified on reuse.
``pick``
    Choose crop boxes (``crop.bbox_um``) for entries that do not have one yet,
    from the molecule/cell density histograms of the source bundle, applying
    the per-crop ``crop.criteria`` (see README).  With
    ``--write`` the chosen boxes are saved into the manifest.
``build``
    Build ``$BAYSOR_BENCH_DATA/real/<id>/`` per the dataset contract:
    ``molecules.parquet``, ``meta.json``, ``reference/`` boundaries, cropped
    focus images where applicable, the nucleus-label image prior where
    applicable, and a provenance ``README.md``.  sha256 of the two contract
    files are recorded into the manifest (``outputs``).
``record-hashes``
    Record/refresh the ``outputs`` sha256 of already-built datasets without
    rebuilding them.
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


def _group_members(group: list[dict], include_images: bool = True) -> list[str]:
    """Union of the member paths of every dataset in a source group."""
    members: list[str] = []
    for ds in group:
        members += _source_members(ds, include_images=include_images)
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


def _vendor_cluster_hists(url: str, group: list[dict], cells, bounds,
                          bin_um: float) -> list | None:
    """Per-cluster vendor cell histograms (composition criterion) or None."""
    member = next(
        (ds["members"].get("clusters") for ds in group
         if ds["members"].get("clusters")),
        None,
    )
    if member is None:
        return None
    clusters = xc.load_cell_clusters(xc.member_cache_path(url, member), cells)
    labels = np.unique(clusters)
    labels = labels[labels >= 0]
    if len(labels) == 0:
        return None
    print(f"  vendor clusters: {len(labels)} (from {member})")
    return [
        xc.cells_histogram(cells, bounds, bin_um, mask=(clusters == k))
        for k in labels
    ]


# ---------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------


def cmd_verify(args, manifest) -> int:
    """Check source URLs, cached-member hashes and built-dataset hashes."""
    ok = True
    if not args.skip_urls:
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

    n_members = n_outputs = 0
    for zip_name, members in sorted(
        (manifest.get("downloaded_members") or {}).items()
    ):
        for member, spec in sorted(members.items()):
            path = xc.cache_root() / zip_name / member
            tag = f"member {zip_name}/{member}"
            if not path.exists():
                print(f"[skip] {tag}: not cached")
                continue
            n_members += 1
            size = path.stat().st_size
            if size != spec["bytes"]:
                print(f"[FAIL] {tag}: size {size} != {spec['bytes']}")
                ok = False
                continue
            digest = xc.download.sha256_file(path)
            if digest != spec["sha256"]:
                print(f"[FAIL] {tag}: sha256 {digest} != {spec['sha256']}")
                ok = False
            else:
                print(f"[ ok ] {tag}: sha256 {digest[:12]}… ({size} bytes)")

    root = Path(args.data_root) / "real" if args.data_root \
        else xc.data_root() / "real"
    for ds in manifest["datasets"]:
        for fname, want in sorted((ds.get("outputs") or {}).items()):
            path = root / ds["id"] / fname
            tag = f"output {ds['id']}/{fname}"
            if not path.exists():
                if args.require_built:
                    print(f"[FAIL] {tag}: missing")
                    ok = False
                else:
                    print(f"[skip] {tag}: not built")
                continue
            n_outputs += 1
            digest = xc.download.sha256_file(path)
            if digest != want:
                print(f"[FAIL] {tag}: sha256 {digest} != {want}")
                ok = False
            else:
                print(f"[ ok ] {tag}: sha256 {digest[:12]}…")
    print(f"verified {n_members} cached member(s), {n_outputs} output file(s)"
          + ("" if ok else " — FAILURES ABOVE"))
    return 0 if ok else 1


def cmd_fetch(args, manifest) -> int:
    datasets = _datasets(manifest, args.ids)
    hashes = manifest.setdefault("downloaded_members", {})
    changed = False
    total = 0
    for url, group in _group_by_source(datasets).items():
        print(f"source {url}", flush=True)
        zip_name = Path(url).name[:-4]
        record: dict = {}
        paths = xc.fetch_members(
            url, _group_members(group),
            expected=dict(hashes.get(zip_name, {})), record=record,
        )
        merged = dict(hashes.get(zip_name, {}))
        merged.update(record)
        if hashes.get(zip_name) != merged:
            hashes[zip_name] = merged
            changed = True
        total += sum(p.stat().st_size for p in set(paths.values()))
    if changed:
        xc.save_manifest(manifest, args.manifest)
        print(f"manifest member hashes updated: {args.manifest}")
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
        xc.fetch_members(url, _group_members(group, include_images=False))
        cells = _load_cells(_source_with_members(any_ds))
        if not needs_hist:
            continue
        mol, cel, bounds, bin_um = _hist_cache(
            _source_with_members(any_ds), cells, manifest["defaults"]["min_qv"], bin_um
        )
        cluster_hists = None
        if any(d["crop"].get("criteria", {}).get("min_clusters")
               for d in needs_hist):
            cluster_hists = _vendor_cluster_hists(url, group, cells, bounds, bin_um)
        for ds in needs_hist:
            crop = ds["crop"]
            criteria = dict(crop.get("criteria") or {})
            within = None
            if crop.get("within"):
                ref = crop["within"]
                within = bbox_by_id.get(ref)
                if within is None:
                    raise SystemExit(
                        f"{ds['id']}: within={ref} has no bbox yet; list it first"
                    )
            if criteria.get("min_clusters") and cluster_hists is None:
                raise SystemExit(
                    f"{ds['id']}: criteria.min_clusters needs a 'clusters' "
                    f"member for this source"
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
                cluster_hists=(cluster_hists if criteria.get("min_clusters")
                               else None),
                criteria=criteria,
            )
            print(
                f"{ds['id']}: bbox={res['bbox_um']} mols={res['mols']} "
                f"cells={res['cells']} cov={res['coverage']:.2f} "
                f"dens={res['cells_per_mm2']:.0f}/mm2 hint={res['density_hint']}"
                + (f" ncl={res['n_clusters']}" if "n_clusters" in res else "")
                + (f" empty_border={res['empty_border']}"
                   if "empty_border" in res else "")
                + (" RELAXED" if res["relaxed_density"] else "")
            )
            crop["bbox_um"] = res["bbox_um"]
            crop["pick"] = {
                k: res[k]
                for k in ("mols", "cells", "coverage", "cells_per_mm2",
                          "density_hint", "relaxed_density", "n_clusters",
                          "empty_border")
                if k in res
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
    out_base = Path(args.out_root) if args.out_root \
        else xc.data_root() / "real"
    record = not args.no_record_hashes and not args.out_root
    built = []
    changed = False
    for url, group in _group_by_source(datasets).items():
        any_ds = group[0]
        print(f"source {url}", flush=True)
        xc.fetch_members(url, _group_members(group))
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
        with_z = any(ds["crop"].get("keep_z") for ds in group)
        print(f"  reading {len(bboxes)} crop(s) from transcripts"
              + (" (with z)" if with_z else "") + " ...", flush=True)
        crops = xc.read_transcript_crops(tp, bboxes, min_qv, with_z=with_z)
        for ds, crop_tbl in zip(group, crops):
            built.append(_build_one(ds, crop_tbl, cells, defaults,
                                    out_base=out_base))
            if record:
                ds_dir = out_base / ds["id"]
                outs = {
                    fname: xc.download.sha256_file(ds_dir / fname)
                    for fname in ("molecules.parquet", "meta.json")
                }
                if ds.get("outputs") != outs:
                    ds["outputs"] = outs
                    changed = True
    print(f"built {len(built)} dataset(s):")
    for line in built:
        print("  " + line)
    if changed:
        xc.save_manifest(manifest, args.manifest)
        print(f"manifest output hashes updated: {args.manifest}")
    elif args.out_root and not args.no_record_hashes:
        print("note: --out-root given, output hashes not recorded")
    return 0


def cmd_record_hashes(args, manifest) -> int:
    """Record sha256 of each built dataset's molecules.parquet and meta.json."""
    root = Path(args.data_root) / "real" if args.data_root \
        else xc.data_root() / "real"
    ids = set(args.ids)
    changed = False
    recorded = missing = 0
    for ds in manifest["datasets"]:
        if ids and ds["id"] not in ids:
            continue
        files = {name: root / ds["id"] / name
                 for name in ("molecules.parquet", "meta.json")}
        absent = [name for name, path in files.items() if not path.is_file()]
        if absent:
            missing += 1
            msg = f"{ds['id']}: not built ({', '.join(absent)})"
            print(f"[skip] {msg}")
            if args.require_built:
                raise SystemExit(msg)
            continue
        outs = {name: xc.download.sha256_file(path)
                for name, path in files.items()}
        if ds.get("outputs") != outs:
            ds["outputs"] = outs
            changed = True
        recorded += 1
        print(f"[ ok ] {ds['id']}: "
              + " ".join(f"{n}={h[:12]}…" for n, h in outs.items()))
    if changed:
        xc.save_manifest(manifest, args.manifest)
        print(f"manifest output hashes updated: {args.manifest}")
    else:
        print("manifest output hashes unchanged")
    print(f"recorded {recorded} dataset(s), {missing} not built")
    return 0


def _build_one(ds: dict, crop_tbl, cells, defaults: dict,
               out_base: Path) -> str:
    import pyarrow.parquet as pq

    src = ds["source"]
    bbox = ds["crop"]["bbox_um"]
    ds_dir = out_base / ds["id"]
    ds_dir.mkdir(parents=True, exist_ok=True)

    keep_z = bool(ds["crop"].get("keep_z"))
    mol = xc.build_molecule_table(crop_tbl, bbox, keep_z=keep_z)
    z_range = z_std = None
    if keep_z and "z" in mol.column_names:
        z = np.asarray(mol["z"].to_numpy(zero_copy_only=False), dtype=np.float64)
        z = z[np.isfinite(z)]
        if len(z):
            z_range = (float(z.min()), float(z.max()))
            z_std = float(z.std())
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
        images = _build_images(ds, bbox, out_base)

    # image prior: vendor nucleus boundaries rasterised into the crop's
    # pixel frame (1 µm/px, absolute origin) - see rasterize_nucleus_labels
    if ds.get("prior_image"):
        spec_in = ds["prior_image"]
        member_key = spec_in.get("from", "nucleus_boundaries")
        rel = spec_in.get("file", "images/nuclei_labels.tif")
        bp = xc.member_cache_path(src["url"], ds["members"][member_key])
        raster = xc.rasterize_nucleus_labels(bp, bbox, ds_dir / rel)
        images.append({
            "name": "nuclei_labels",
            "file": rel,
            "pixel_size_um": raster["pixel_size_um"],
            "origin_um": raster["origin_um"],
            "source": (f"{member_key} of the source bundle rasterised at "
                       f"1 µm/px over the absolute coordinate origin; "
                       f"{raster['n_labels']} labels, "
                       f"{raster['n_labelled_pixels']} labelled pixels"),
        })
        print(f"  prior image {rel}: {raster['width']}x{raster['height']} "
              f"{raster['n_labels']} labels, "
              f"{raster['n_labelled_pixels']} labelled px")

    note = ds["crop"]["note"]
    if z_range is not None:
        note += (f"; keeps z: {z_range[0]:.2f}-{z_range[1]:.2f} µm "
                 f"(std {z_std:.2f} µm over {stats['n_molecules']} molecules)")
    baysor = {**defaults["baysor"], **(ds.get("baysor") or {})}
    meta = xc.build_meta(
        dataset_id=ds["id"],
        tier=ds["tier"],
        source=src,
        bbox=bbox,
        note=note,
        stats=stats,
        panel_genes=panel_genes,
        scale_um=scale_um,
        scale_method=scale_method,
        baysor=baysor,
        images=images,
        retrieved=src.get("retrieved", manifest_retrieved),
        difficulty_notes=ds.get("difficulty_notes", note),
        z_range=z_range,
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


def _build_images(ds: dict, bbox: Sequence, out_base: Path) -> list[dict]:
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
        xc.write_tif(out_base / ds["id"] / "images" / f"{name}.tif", arr)
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

    p = sub.add_parser("verify",
                       help="HEAD-check source URLs and verify recorded hashes")
    p.add_argument("--skip-urls", action="store_true",
                   help="skip the network URL checks (hash checks only)")
    p.add_argument("--require-built", action="store_true",
                   help="fail when a dataset with recorded outputs is not built")
    p.add_argument("--data-root", default=None,
                   help="data root containing real/<id>/ (default: $BAYSOR_BENCH_DATA)")

    p = sub.add_parser("fetch", help="download manifest members into the cache")
    p.add_argument("ids", nargs="*", help="dataset ids (default: all)")

    p = sub.add_parser("pick", help="choose crop bboxes")
    p.add_argument("ids", nargs="*", help="dataset ids (default: all)")
    p.add_argument("--write", action="store_true", help="save bboxes into the manifest")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--bin-um", type=float, default=25.0)

    p = sub.add_parser("build", help="build dataset directories")
    p.add_argument("ids", nargs="*", help="dataset ids (default: all)")
    p.add_argument("--out-root", default=None,
                   help="write datasets under this directory instead of "
                        "<data>/real (also disables output-hash recording)")
    p.add_argument("--no-record-hashes", action="store_true",
                   help="do not update output sha256 hashes in the manifest")

    p = sub.add_parser(
        "record-hashes",
        help="record sha256 of built molecules.parquet/meta.json in the manifest",
    )
    p.add_argument("ids", nargs="*", help="dataset ids (default: all)")
    p.add_argument("--data-root", default=None,
                   help="data root containing real/<id>/ (default: $BAYSOR_BENCH_DATA)")
    p.add_argument("--require-built", action="store_true",
                   help="fail instead of skipping datasets that are not built")

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
        "record-hashes": cmd_record_hashes,
        "report": cmd_report,
    }[args.cmd](args, manifest)


if __name__ == "__main__":
    raise SystemExit(main())
