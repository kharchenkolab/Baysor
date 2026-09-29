"""Regenerate every simulated dataset listed in ``benchmarks/datasets/sim.yaml``.

Deterministic: each dataset is rebuilt from its manifest entry (generator,
scenario, parameters, seed) and written to ``$BAYSOR_BENCH_DATA/sim/<id>/``.
Re-running on a fresh machine with the same manifest produces byte-identical
``molecules.parquet`` and ``meta.json``.

Three generators:

``trivial`` / ``strec``
    Built by :mod:`trivial` and :mod:`strec` from ``scenario`` / ``params``.
``noprior``
    Derived from another entry (``base:``): the *same molecules* with the
    ``prior`` column dropped and ``meta.baysor.prior = "none"``, so Baysor's
    default no-prior mode is measured against the same truth.

Every manifest entry carries the expected SHA-256 of both files under
``sha256:``; ``--verify-all`` checks the datasets on disk *and* the fresh
regeneration against those hashes (not only against each other), and
``--update-hashes`` (re)records them.

Usage::

    python benchmarks/simulate/generate_all.py                 # everything
    python benchmarks/simulate/generate_all.py --list          # show the manifest
    python benchmarks/simulate/generate_all.py --only sim_circles_gaps_g100
    python benchmarks/simulate/generate_all.py --verify sim_circles_gaps_g100
    python benchmarks/simulate/generate_all.py --verify-all    # vs manifest hashes
    python benchmarks/simulate/generate_all.py --update-hashes # record hashes
"""
from __future__ import annotations

import argparse
import copy
import json
import re
import sys
import tempfile
import time
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import common  # noqa: E402
import strec  # noqa: E402
import trivial  # noqa: E402

DEFAULT_MANIFEST = HERE.parent / "datasets" / "sim.yaml"
_HASH_KEYS = ("molecules.parquet", "meta.json")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")


def load_manifest(path: Path | str = DEFAULT_MANIFEST) -> dict:
    with open(path) as fh:
        manifest = yaml.safe_load(fh)
    if manifest.get("version") != 1:
        raise ValueError(f"unsupported manifest version in {path}")
    ids = [d["id"] for d in manifest["datasets"]]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate dataset ids in manifest")
    return manifest


def build(entry: dict, index: dict | None = None) -> tuple:
    """Build one manifest entry -> ``(molecules_df, meta)``.

    ``index`` maps id -> entry and is required for ``noprior`` entries, which
    rebuild their ``base`` and derive from it.
    """
    gen = entry["generator"]
    if gen == "noprior":
        if index is None:
            raise KeyError(f"{entry['id']}: noprior build needs the entry index")
        base = index.get(entry.get("base"))
        if base is None:
            raise KeyError(f"{entry['id']}: unknown base {entry.get('base')!r}")
        if base["generator"] == "noprior":
            raise ValueError(f"{entry['id']}: base must not be a noprior entry")
        df, meta = build(base, index)
        if "prior" in df.columns:
            df = df.drop(columns=["prior"])
        meta = copy.deepcopy(meta)
        meta["id"] = entry["id"]
        meta["tier"] = entry["tier"]
        meta["baysor"]["prior"] = "none"
        # dataset id recorded inside truth.params by the trivial generator
        tp = meta.get("truth", {}).get("params")
        if isinstance(tp, dict) and "dataset_id" in tp:
            tp["dataset_id"] = entry["id"]
            tp["tier"] = entry["tier"]
        meta["truth"]["prior_variant_of"] = base["id"]
        note = (f"noprior variant of {base['id']}: identical molecules (truth "
                "and positions), prior column dropped, baysor.prior = none")
        meta["source"]["note"] = f"{meta['source']['note']}; {note}"
        meta["difficulty"]["notes"] = (
            f"{meta['difficulty']['notes']}, no prior (noprior variant of "
            f"{base['id']})")
        return df, meta
    if gen == "trivial":
        return trivial.generate(entry["scenario"], seed=entry["seed"],
                                dataset_id=entry["id"], tier=entry["tier"],
                                params=entry.get("params") or {})
    if gen == "strec":
        p = entry.get("params") or {}
        z = p.get("z_slab_um")
        return strec.generate(dataset_id=entry["id"], tier=entry["tier"],
                              packing=float(p["packing"]), sigma=float(p["sigma"]),
                              model_kind=p["model"], seed=entry["seed"],
                              mean_tx=float(p.get("mean_tx", strec.MEAN_TX_DEFAULT)),
                              n_target=int(p.get("n_target", strec.N_TARGET_DEFAULT)),
                              geometry=p.get("geometry", "voronoi"),
                              bg_frac=float(p.get("bg_frac", 0.0)),
                              z_slab_um=(float(z) if z is not None else None),
                              prior_opts=p.get("prior_opts"))
    raise KeyError(f"unknown generator {gen!r}")


def validate_manifest(manifest: dict) -> dict:
    """Structural checks: base links, prior_opts provenance, hash format.

    Returns ``{id: entry}``.  Raises ``ValueError``/``KeyError`` on violation.
    """
    index = {d["id"]: d for d in manifest["datasets"]}
    for entry in manifest["datasets"]:
        eid = entry["id"]
        if entry["generator"] not in manifest["generators"]:
            raise KeyError(f"{eid}: unknown generator {entry['generator']!r}")
        if entry["tier"] not in ("quick", "full"):
            raise ValueError(f"{eid}: bad tier {entry['tier']!r}")
        if not isinstance(entry["seed"], int):
            raise ValueError(f"{eid}: seed must be an int")
        hashes = entry.get("sha256")
        if hashes is not None:
            if not isinstance(hashes, dict) or set(hashes) != set(_HASH_KEYS):
                raise ValueError(f"{eid}: sha256 must cover exactly {_HASH_KEYS}")
            for k, v in hashes.items():
                if not isinstance(v, str) or not _HEX64.match(v):
                    raise ValueError(f"{eid}: bad sha256 for {k}")
        if entry["generator"] == "noprior":
            base = index.get(entry.get("base"))
            if base is None:
                raise KeyError(f"{eid}: unknown base {entry.get('base')!r}")
            if base["generator"] == "noprior":
                raise ValueError(f"{eid}: base must not be a noprior entry")
            if base["tier"] != entry["tier"]:
                raise ValueError(f"{eid}: tier must match base {base['id']}")
            if base["seed"] != entry["seed"]:
                raise ValueError(f"{eid}: seed must match base {base['id']}")
        params = entry.get("params") or {}
        opts = params.get("prior_opts")
        if opts is not None:
            if opts.get("kind", "imperfect") != "imperfect":
                raise ValueError(f"{eid}: unknown prior_opts kind")
            if "seed" not in opts:
                raise ValueError(f"{eid}: prior_opts needs a seed")
            base = index.get(opts.get("base"))
            if base is None:
                raise KeyError(f"{eid}: prior_opts.base {opts.get('base')!r} not found")
            if base["generator"] != entry["generator"]:
                raise ValueError(f"{eid}: prior_opts.base generator mismatch")
            if base["seed"] != entry["seed"]:
                raise ValueError(f"{eid}: prior_opts.base seed mismatch")
            if entry["generator"] == "trivial" and base.get("scenario") != entry.get("scenario"):
                raise ValueError(f"{eid}: prior_opts.base scenario mismatch")
            own = {k: v for k, v in params.items() if k != "prior_opts"}
            their = {k: v for k, v in (base.get("params") or {}).items()
                     if k != "prior_opts"}
            if own != their:
                raise ValueError(f"{eid}: params differ from base {base['id']}: "
                                 f"{own} != {their}")
    return index


def generate_one(entry: dict, out_root: Path, index: dict | None = None) -> dict:
    t0 = time.time()
    df, meta = build(entry, index)
    out_dir = out_root / entry["id"]
    hashes = common.write_dataset(out_dir, df, meta)
    result = {"id": entry["id"], "n_molecules": int(len(df)),
              "n_genes": int(df["gene"].nunique()), "tier": entry["tier"],
              "seconds": round(time.time() - t0, 2), "out": str(out_dir),
              "sha256": hashes}
    expected = entry.get("sha256")
    if expected is not None:
        result["matches_manifest"] = all(hashes[k] == expected[k] for k in _HASH_KEYS)
    return result


def verify_one(entry: dict, out_root: Path, index: dict | None = None) -> dict:
    """Regenerate into a temp dir and check disk + regeneration against the
    manifest's ``sha256`` (the manifest is the source of truth)."""
    expected = entry.get("sha256") or {}
    problems: list[str] = []
    for k in _HASH_KEYS:
        if not expected.get(k):
            problems.append(f"manifest entry has no sha256 for {k}")
    with tempfile.TemporaryDirectory(prefix=f"verify_{entry['id']}_") as tmp:
        df, meta = build(entry, index)
        regen = common.write_dataset(Path(tmp) / entry["id"], df, meta)
    for k in _HASH_KEYS:
        if expected.get(k) and regen[k] != expected[k]:
            problems.append(f"{k}: regeneration differs from manifest")
    existing = out_root / entry["id"]
    disk: dict[str, str] = {}
    if not existing.exists():
        problems.append("dataset not on disk")
    else:
        for k in _HASH_KEYS:
            f = existing / k
            if not f.exists():
                problems.append(f"{k}: missing on disk")
            else:
                disk[k] = common.sha256(f)
                ref = expected.get(k) or regen[k]
                if disk[k] != ref:
                    problems.append(f"{k}: on-disk hash differs from "
                                    + ("manifest" if expected.get(k) else "regeneration"))
    return {"id": entry["id"], "ok": not problems, "mismatch": problems,
            "sha256": regen}


# ---------------------------------------------------------------------------
# Manifest hash bookkeeping (line-based, to keep comments/formatting intact)
# ---------------------------------------------------------------------------

def _entry_blocks(lines: list[str]) -> dict[str, tuple[int, int]]:
    """``{id: (start, end)}`` line ranges of each ``  - id: ...`` block."""
    starts = [(i, m.group(1)) for i, ln in enumerate(lines)
              if (m := re.match(r"^  - id: (\S+)$", ln))]
    blocks: dict[str, tuple[int, int]] = {}
    for k, (i, eid) in enumerate(starts):
        end = starts[k + 1][0] if k + 1 < len(starts) else len(lines)
        blocks[eid] = (i, end)
    return blocks


def update_hashes(manifest_path: Path, hashes_by_id: dict[str, dict]) -> list[str]:
    """Insert/replace the ``sha256:`` line of each updated entry, in place.

    Returns the list of ids whose line changed.  Comments and formatting of
    the rest of the manifest are preserved (purely textual edit).
    """
    lines = manifest_path.read_text().splitlines()
    blocks = _entry_blocks(lines)
    changed: list[str] = []
    for eid in sorted(hashes_by_id, key=lambda i: -blocks[i][0] if i in blocks else 0):
        if eid not in blocks:
            raise KeyError(f"{eid}: not in {manifest_path}")
        h = hashes_by_id[eid]
        line = (f"    sha256: {{molecules.parquet: {h['molecules.parquet']}, "
                f"meta.json: {h['meta.json']}}}")
        start, end = blocks[eid]
        existing = [i for i in range(start, end)
                    if lines[i].startswith("    sha256: ")]
        if existing:
            if lines[existing[0]] != line:
                lines[existing[0]] = line
                changed.append(eid)
            continue
        j = end
        while j > start + 1 and (not lines[j - 1].strip()
                                 or lines[j - 1].lstrip().startswith("#")):
            j -= 1
        lines.insert(j, line)
        changed.append(eid)
    text = "\n".join(lines) + "\n"
    manifest_path.write_text(text)
    # round-trip: the edited manifest must still parse with the same entries
    reparsed = load_manifest(manifest_path)
    validate_manifest(reparsed)
    for eid, h in hashes_by_id.items():
        got = next(d.get("sha256") for d in reparsed["datasets"] if d["id"] == eid)
        if got != h:
            raise RuntimeError(f"{eid}: sha256 did not survive the manifest edit")
    return changed


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    p.add_argument("--out-root", default=None,
                   help="default: $BAYSOR_BENCH_DATA/sim")
    p.add_argument("--only", action="append", default=None,
                   help="dataset id (repeatable); default: all")
    p.add_argument("--list", action="store_true", help="list manifest entries")
    p.add_argument("--verify", action="append", default=None, metavar="ID",
                   help="regenerate and compare hashes with the manifest + disk")
    p.add_argument("--verify-all", action="store_true",
                   help="verify every dataset in the manifest against its "
                        "committed sha256 (disk and regeneration)")
    p.add_argument("--update-hashes", action="store_true",
                   help="(re)record sha256 of generated files in the manifest")
    args = p.parse_args(argv)

    manifest = load_manifest(args.manifest)
    index = validate_manifest(manifest)
    entries = manifest["datasets"]
    if args.only:
        wanted = set(args.only)
        unknown = wanted - {e["id"] for e in entries}
        if unknown:
            p.error(f"unknown ids: {sorted(unknown)}")
        entries = [e for e in entries if e["id"] in wanted]

    if args.list:
        for e in entries:
            params = " ".join(f"{k}={v}" for k, v in (e.get("params") or {}).items())
            has_hash = "hash" if e.get("sha256") else "NOHASH"
            print(f"{e['id']:34s} {e['generator']:9s} "
                  f"{e.get('scenario', ''):16s} tier={e['tier']:5s} "
                  f"seed={e['seed']} {has_hash} {params}")
        return 0

    out_root = Path(args.out_root) if args.out_root else common.data_root() / "sim"

    if args.verify or args.verify_all:
        targets = entries if args.verify_all else [
            e for e in entries if e["id"] in set(args.verify)]
        results = [verify_one(e, out_root, index) for e in targets]
        print(json.dumps(results, indent=2))
        return 0 if all(r["ok"] for r in results) else 1

    results = []
    for entry in entries:
        info = generate_one(entry, out_root, index)
        results.append(info)
        flag = ""
        if "matches_manifest" in info:
            flag = "  OK" if info["matches_manifest"] else "  MISMATCH vs manifest"
        print(f"{info['id']:34s} {info['n_molecules']:8d} molecules "
              f"{info['n_genes']:5d} genes  {info['seconds']:7.2f}s{flag}")

    stale = [r["id"] for r in results if r.get("matches_manifest") is False]
    if args.update_hashes:
        changed = update_hashes(Path(args.manifest),
                                {r["id"]: r["sha256"] for r in results})
        print(f"recorded sha256 for {len(results)} datasets "
              f"({len(changed)} newly written or changed)")
        if changed and not args.only:
            # full regeneration: everything else must have matched
            pass
        stale = [i for i in stale if i not in changed] if args.only else []
    elif stale:
        print(f"ERROR: regenerated files differ from the manifest hashes for: "
              f"{stale}; run --update-hashes if the change is intended")

    summary = {"n_datasets": len(results),
               "n_molecules_total": sum(r["n_molecules"] for r in results),
               "out_root": str(out_root)}
    print(json.dumps(summary, indent=2))
    if stale:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
