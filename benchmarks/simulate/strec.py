"""Simulated datasets from the *st-recoverability* generator (Yamasaki, 2026).

Imports ``generator.py``, ``expression.py``, ``oracle.py`` and ``config.py``
from a pinned clone of https://github.com/RuiYamasaki/st-recoverability
(MIT, commit ``e85faeccfb58...``) via ``sys.path`` — the code is never vendored
into this repository.  The clone and the two small reference downloads live
under ``$BAYSOR_BENCH_DATA/cache/sim/``.

Pipeline per dataset (mirrors the upstream ``src/demo.py`` and
``src/methods_baysor.py``):

  1. ``build_field(packing, sigma, seed, model=...)``  — Voronoi field of 400
     cells at the requested packing and displacement sigma;
  2. ``generate_transcripts(field, mean_tx_per_cell, seed)``;
  3. ``build_oracle_maps`` + ``oracle_assign`` (Bayes-optimal) and
     ``naive_assign`` (nearest nucleus) — both accuracies are stored in
     ``meta.truth``;
  4. contract mapping: ``true_cell`` -> ``cell`` (+1, 0 = background),
     ``interior`` -> ``interior``, observed coordinates -> ``x``/``y`` (um —
     checked against the generator's documented units), nucleus prior (3 um
     around the true centre, as in ``methods_baysor.py``) -> ``prior``,
     type name -> ``celltype``.

Realistic expression models:
  * MERFISH: ``build_realistic_model_from_merfish`` on the squidpy
    ``MERFISH_0.24.h5ad`` mirror (figshare) that upstream ``realism.py``
    downloads;
  * Xenium: ``build_realistic_model_from_xenium`` on the cell-feature matrix
    of Xenium breast cancer Rep 1 (10x release 1.0.1); only
    ``cell_feature_matrix.h5`` and ``cells.parquet`` are fetched from the 9.9
    GB ``outs.zip`` with ``remotezip`` range requests.

CLI (one dataset)::

    python benchmarks/simulate/strec.py --id strec_dense_s2_merfish \\
        --packing 13625 --sigma 2.0 --model merfish --seed 910002 \\
        -o $BAYSOR_BENCH_DATA/sim/strec_dense_s2_merfish
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from . import common
except ImportError:  # running as a plain script
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import common

GENERATOR_VERSION = 1

REPO_URL = "https://github.com/RuiYamasaki/st-recoverability"
REPO_COMMIT = "e85faeccfb58e0bef05d76817e98c79abff6bf51"  # pin "e85faeccfb58"
REPO_MODULES = ["generator", "expression", "oracle", "config"]

MERFISH_H5AD_URL = "https://ndownloader.figshare.com/files/40038538"  # squidpy MERFISH_0.24.h5ad
XENIUM_OUTS_ZIP = (
    "https://cf.10xgenomics.com/samples/xenium/1.0.1/"
    "Xenium_FFPE_Human_Breast_Cancer_Rep1/"
    "Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip"
)

# upstream generator constants (recorded in truth)
N_TARGET_DEFAULT = 400          # generator.N_TARGET cells per field
MEAN_TX_DEFAULT = 160.0         # demo.DENSITY
NUCLEUS_RADIUS_UM = 3.0         # methods_baysor.NUCLEUS_RADIUS_UM
MODEL_SEED_DEFAULT = 20260618   # config.MASTER_SEED (Xenium k-means seed)


def data_root() -> Path:
    """``$BAYSOR_BENCH_DATA`` (default ``<repo>/.bench-data``)."""
    return common.data_root()


def cache_dir() -> Path:
    return common.cache_dir()


# ---------------------------------------------------------------------------
# External repo + reference data
# ---------------------------------------------------------------------------

def ensure_repo() -> Path:
    """Clone the pinned st-recoverability commit into the cache (idempotent)."""
    repo = cache_dir() / "st-recoverability"
    if not (repo / ".git").exists():
        repo.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "clone", REPO_URL, str(repo)], check=True)
    head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    if head != REPO_COMMIT:
        try:
            subprocess.run(["git", "-C", str(repo), "cat-file", "-e", REPO_COMMIT],
                           check=True, capture_output=True)
        except subprocess.CalledProcessError:
            subprocess.run(["git", "-C", str(repo), "fetch", "origin", REPO_COMMIT],
                           check=True)
        subprocess.run(["git", "-C", str(repo), "checkout", "--detach", REPO_COMMIT],
                       check=True, capture_output=True)
    return repo


def _import_strec():
    """Import the pinned upstream modules through ``sys.path``."""
    src = ensure_repo() / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
    import config
    import expression
    import generator
    import oracle
    return generator, expression, oracle, config


def ensure_merfish() -> Path:
    """Download the squidpy MERFISH_0.24.h5ad mirror (cached, ~4 MB)."""
    dest = cache_dir() / "merfish_moffitt.h5ad"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(MERFISH_H5AD_URL,
                                 headers={"User-Agent": "baysor-bench-sim/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r, open(dest, "wb") as fh:
        while True:
            block = r.read(1 << 20)
            if not block:
                break
            fh.write(block)
    return dest


def ensure_xenium() -> tuple[Path, Path]:
    """Fetch only ``cell_feature_matrix.h5`` + ``cells.parquet`` from the 10x
    Xenium breast Rep1 ``outs.zip`` using HTTP range requests (remotezip)."""
    from remotezip import RemoteZip

    h5 = cache_dir() / "xenium_breast_rep1_cell_feature_matrix.h5"
    cells = cache_dir() / "xenium_breast_rep1_cells.parquet"
    if h5.exists() and h5.stat().st_size > 0 and cells.exists() and cells.stat().st_size > 0:
        return h5, cells
    cache_dir().mkdir(parents=True, exist_ok=True)
    with RemoteZip(XENIUM_OUTS_ZIP) as z:
        for member, dest in [("outs/cell_feature_matrix.h5", h5),
                             ("outs/cells.parquet", cells)]:
            if not (dest.exists() and dest.stat().st_size > 0):
                z.extract(member, cache_dir())
                os.replace(cache_dir() / member, dest)
    return h5, cells


def build_model(kind: str):
    """Expression model: ``disjoint`` (toy), ``merfish`` or ``xenium``."""
    generator, expression, oracle, config = _import_strec()
    if kind == "disjoint":
        return expression.build_disjoint_model(), {}
    if kind == "merfish":
        path = ensure_merfish()
        return expression.build_realistic_model_from_merfish(h5ad_path=str(path)), {
            "reference": "squidpy MERFISH_0.24.h5ad (figshare 40038538)",
        }
    if kind == "xenium":
        h5, cells = ensure_xenium()
        model, _real = expression.build_realistic_model_from_xenium(
            h5_path=str(h5), cells_parquet=str(cells),
            seed=MODEL_SEED_DEFAULT)
        return model, {
            "reference": "Xenium_FFPE_Human_Breast_Cancer_Rep1 cell_feature_matrix.h5 "
                         "(10x release 1.0.1)",
            "model_seed": MODEL_SEED_DEFAULT,
        }
    raise KeyError(f"unknown model {kind!r}; have disjoint, merfish, xenium")


# ---------------------------------------------------------------------------
# Dataset generation
# ---------------------------------------------------------------------------

def _check_units(field) -> dict:
    """Guard that the upstream generator really works in micrometres."""
    if not (10.0 <= field.L_um <= 5000.0):
        raise ValueError(f"field side {field.L_um} is not in um (check units)")
    if not (1.0 <= field.r_mean_um <= 50.0):
        raise ValueError(f"mean radius {field.r_mean_um} is not in um (check units)")
    if not (0.1 <= field.sigma_um <= 10.0):
        raise ValueError(f"sigma {field.sigma_um} is not in um (check units)")
    return {"L_um": float(field.L_um), "r_mean_um": float(field.r_mean_um),
            "sigma_um": float(field.sigma_um), "units": "um (verified)"}


def _accuracy(true_cell, assigned, interior) -> float:
    return float((assigned[interior] == true_cell[interior]).mean())


def generate(*, dataset_id: str, tier: str, packing: float, sigma: float,
             model_kind: str, seed: int, mean_tx: float = MEAN_TX_DEFAULT,
             n_target: int = N_TARGET_DEFAULT) -> tuple[pd.DataFrame, dict]:
    """Build one st-recoverability field and map it to the dataset contract."""
    generator, expression, oracle, config = _import_strec()
    model, model_info = build_model(model_kind)

    seed_field, seed_tx = int(seed), int(seed) + 1
    field = generator.build_field(packing, sigma, seed_field, n_target=n_target,
                                  model=model)
    tx = generator.generate_transcripts(field, mean_tx, seed_tx)

    # units check + oracle / naive accuracies (interior molecules, as in demo.py)
    units = _check_units(field)
    dmax, argcell = oracle.build_oracle_maps(field)
    # upstream oracle takes log(p[gene|type]); zeros are guarded by `pos`
    with np.errstate(divide="ignore"):
        oracle_assign = oracle.oracle_assign(field, dmax, argcell, tx.obs_xy, tx.gene)
    naive_assign = oracle.naive_assign(field, tx.obs_xy)
    oracle_acc = _accuracy(tx.true_cell, oracle_assign, tx.interior)
    naive_acc = _accuracy(tx.true_cell, naive_assign, tx.interior)

    # nucleus prior: within 3 um of the true centre (methods_baysor.py logic)
    prior, _ = common.nucleus_prior(tx.obs_xy, field.centers, NUCLEUS_RADIUS_UM)

    type_names = [str(t) for t in model.type_names]
    df = pd.DataFrame({
        "x": tx.obs_xy[:, 0], "y": tx.obs_xy[:, 1],
        "gene": np.asarray(model.gene_names, dtype=object)[tx.gene],
        "prior": prior,
        "cell": (tx.true_cell + 1).astype(np.int32),
        "interior": tx.interior.astype(bool),
        "celltype": np.asarray(type_names, dtype=object)[tx.true_type],
    })
    df = common.make_molecules(df)

    area_um2 = field.L_um ** 2
    stats = common.make_stats(df, area_um2=area_um2, n_true_cells=field.n_cells)
    n_genes = int(model.n_genes)
    difficulty = {
        "cell_density": common.density_class(stats["true_cells_per_mm2"]),
        "gene_panel": common.gene_panel_class(n_genes),
        "notes": f"st-recoverability field: packing {packing:g} cells/mm2, "
                 f"sigma {sigma:g} um, model {model.name}",
    }
    params = {
        "packing_cells_per_mm2": packing,
        "sigma_um": sigma,
        "model": model_kind,
        "model_name": model.name,
        "model_info": model_info,
        "mean_tx_per_cell": mean_tx,
        "n_target": n_target,
        "geometry": "voronoi",
        "displacement": "gaussian",
        "nucleus_radius_um": NUCLEUS_RADIUS_UM,
    }
    meta = common.make_meta(
        id=dataset_id, tier=tier,
        source={
            "generator": "benchmarks/simulate/strec.py",
            "generator_version": GENERATOR_VERSION,
            "seed": seed,
            "note": "st-recoverability generator imported from a pinned clone; "
                    f"repo {REPO_URL} commit {REPO_COMMIT} (MIT)",
        },
        crop={
            "bbox_um": [0.0, 0.0, float(field.L_um), float(field.L_um)],
            "z_range_um": None,
            "note": "simulated field (2D)",
        },
        stats=stats,
        difficulty=difficulty,
        baysor=common.make_baysor(field.r_mean_um),
        truth={
            "generator": "benchmarks/simulate/strec.py",
            "generator_version": GENERATOR_VERSION,
            "external": {
                "name": "st-recoverability",
                "repo": REPO_URL,
                "commit": REPO_COMMIT,
                "license": "MIT",
                "imported_modules": REPO_MODULES,
            },
            "seed": {"field": seed_field, "transcripts": seed_tx},
            "params": params,
            "units": units,
            "oracle_accuracy": oracle_acc,
            "naive_accuracy": naive_acc,
            "accuracy_subset": "interior molecules (tx.interior)",
        },
    )
    return df, meta


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Generate one st-recoverability simulated dataset "
                    "(contract: benchmarks/README.md).")
    p.add_argument("--id", required=True, help="dataset id (directory name)")
    p.add_argument("--packing", type=float, required=True,
                   help="cell packing, cells/mm^2 (e.g. 2575, 6000, 13625)")
    p.add_argument("--sigma", type=float, required=True,
                   help="displacement sigma in um (e.g. 1.0, 2.0, 3.0)")
    p.add_argument("--model", required=True,
                   choices=["disjoint", "merfish", "xenium"])
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--mean-tx", type=float, default=MEAN_TX_DEFAULT)
    p.add_argument("--n-target", type=int, default=N_TARGET_DEFAULT)
    p.add_argument("--tier", choices=["quick", "full"], default="quick")
    p.add_argument("-o", "--out-dir", required=True)
    args = p.parse_args(argv)

    df, meta = generate(dataset_id=args.id, tier=args.tier, packing=args.packing,
                        sigma=args.sigma, model_kind=args.model, seed=args.seed,
                        mean_tx=args.mean_tx, n_target=args.n_target)
    hashes = common.write_dataset(args.out_dir, df, meta)
    t = meta["truth"]
    print(json.dumps({"id": args.id, "out": args.out_dir,
                      "n_molecules": int(len(df)),
                      "n_genes": int(df["gene"].nunique()),
                      "oracle_accuracy": t["oracle_accuracy"],
                      "naive_accuracy": t["naive_accuracy"],
                      "sha256": hashes}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
