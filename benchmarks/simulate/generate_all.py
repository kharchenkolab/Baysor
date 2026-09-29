"""Regenerate every simulated dataset listed in ``benchmarks/datasets/sim.yaml``.

Deterministic: each dataset is rebuilt from its manifest entry (generator,
scenario, parameters, seed) and written to ``$BAYSOR_BENCH_DATA/sim/<id>/``.
Re-running on a fresh machine with the same manifest produces byte-identical
``molecules.parquet`` and ``meta.json`` (verified with ``--verify``).

Usage::

    python benchmarks/simulate/generate_all.py                 # everything
    python benchmarks/simulate/generate_all.py --list          # show the manifest
    python benchmarks/simulate/generate_all.py --only sim_circles_gaps_g100
    python benchmarks/simulate/generate_all.py --verify sim_circles_gaps_g100
    python benchmarks/simulate/generate_all.py --verify-all    # hash-compare all
"""
from __future__ import annotations

import argparse
import json
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


def load_manifest(path: Path | str = DEFAULT_MANIFEST) -> dict:
    with open(path) as fh:
        manifest = yaml.safe_load(fh)
    if manifest.get("version") != 1:
        raise ValueError(f"unsupported manifest version in {path}")
    ids = [d["id"] for d in manifest["datasets"]]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate dataset ids in manifest")
    return manifest


def build(entry: dict) -> tuple:
    """Build one manifest entry -> ``(molecules_df, meta)``."""
    gen = entry["generator"]
    if gen == "trivial":
        return trivial.generate(entry["scenario"], seed=entry["seed"],
                                dataset_id=entry["id"], tier=entry["tier"],
                                params=entry.get("params") or {})
    if gen == "strec":
        p = entry.get("params") or {}
        return strec.generate(dataset_id=entry["id"], tier=entry["tier"],
                              packing=float(p["packing"]), sigma=float(p["sigma"]),
                              model_kind=p["model"], seed=entry["seed"],
                              mean_tx=float(p.get("mean_tx", strec.MEAN_TX_DEFAULT)),
                              n_target=int(p.get("n_target", strec.N_TARGET_DEFAULT)))
    raise KeyError(f"unknown generator {gen!r}")


def generate_one(entry: dict, out_root: Path) -> dict:
    t0 = time.time()
    df, meta = build(entry)
    out_dir = out_root / entry["id"]
    hashes = common.write_dataset(out_dir, df, meta)
    return {"id": entry["id"], "n_molecules": int(len(df)),
            "n_genes": int(df["gene"].nunique()), "tier": entry["tier"],
            "seconds": round(time.time() - t0, 2), "out": str(out_dir),
            "sha256": hashes}


def verify_one(entry: dict, out_root: Path) -> dict:
    """Regenerate into a temp dir and compare hashes with what is on disk."""
    existing = out_root / entry["id"]
    if not existing.exists():
        return {"id": entry["id"], "ok": False, "error": "dataset not on disk"}
    with tempfile.TemporaryDirectory(prefix=f"verify_{entry['id']}_") as tmp:
        df, meta = build(entry)
        hashes = common.write_dataset(Path(tmp) / entry["id"], df, meta)
    ok, mismatch = True, []
    for name in ("molecules.parquet", "meta.json"):
        disk = existing / name
        if not disk.exists():
            ok, _ = False, mismatch.append(f"{name}: missing on disk")
        elif common.sha256(disk) != hashes[name]:
            ok = False
            mismatch.append(f"{name}: hash differs")
    return {"id": entry["id"], "ok": ok, "mismatch": mismatch, "sha256": hashes}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    p.add_argument("--out-root", default=None,
                   help="default: $BAYSOR_BENCH_DATA/sim")
    p.add_argument("--only", action="append", default=None,
                   help="dataset id (repeatable); default: all")
    p.add_argument("--list", action="store_true", help="list manifest entries")
    p.add_argument("--verify", action="append", default=None, metavar="ID",
                   help="regenerate into a temp dir and compare file hashes")
    p.add_argument("--verify-all", action="store_true",
                   help="verify every dataset in the manifest")
    args = p.parse_args(argv)

    manifest = load_manifest(args.manifest)
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
            print(f"{e['id']:30s} {e['generator']:8s} {e.get('scenario', ''):16s} "
                  f"tier={e['tier']:5s} seed={e['seed']} {params}")
        return 0

    out_root = Path(args.out_root) if args.out_root else common.data_root() / "sim"

    if args.verify or args.verify_all:
        targets = entries if args.verify_all else [
            e for e in entries if e["id"] in set(args.verify)]
        results = [verify_one(e, out_root) for e in targets]
        print(json.dumps(results, indent=2))
        return 0 if all(r["ok"] for r in results) else 1

    results = []
    for entry in entries:
        info = generate_one(entry, out_root)
        results.append(info)
        print(f"{info['id']:30s} {info['n_molecules']:8d} molecules "
              f"{info['n_genes']:5d} genes  {info['seconds']:7.2f}s")
    summary = {"n_datasets": len(results),
               "n_molecules_total": sum(r["n_molecules"] for r in results),
               "out_root": str(out_root)}
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
