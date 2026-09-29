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

Wrapper extensions (optional, defaults reproduce the upstream datasets
byte-for-byte):

  * ``geometry="aniso"`` — the upstream anisotropic label mode: elongated,
    irregular, randomly oriented cells instead of Voronoi polygons;
  * ``bg_frac`` — ambient background molecules (``cell = 0``) added uniformly
    over the domain by this wrapper (the upstream generator emits only cell
    transcripts);
  * ``z_slab_um`` — a wrapper-added third dimension with displacement (the
    upstream generator is 2D); the 2D oracle is then not claimed;
  * ``prior_opts`` — an imperfect nucleus prior (missed / shifted / merged
    nuclei) via ``common.imperfect_nucleus_prior``; truth unchanged.

Realistic expression models:
  * MERFISH: ``build_realistic_model_from_merfish`` on the squidpy
    ``MERFISH_0.24.h5ad`` mirror (figshare) that upstream ``realism.py``
    downloads;
  * Xenium: ``build_realistic_model_from_xenium`` on the cell-feature matrix
    of Xenium breast cancer Rep 1 (10x release 1.0.1); only
    ``cell_feature_matrix.h5`` and ``cells.parquet`` are fetched from the 9.9
    GB ``outs.zip`` with ``remotezip`` range requests;
  * Prime 5K (``prime5k1000`` / ``prime5k5000``):
    ``build_model_from_xenium_h5`` on the Xenium Prime 5K ovarian
    ``cell_feature_matrix.h5`` (10x release 3.0.0, fetched via remotezip),
    a seeded 30 k-cell subsample with the panel subsetted to the top-N genes
    by mean expression.

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

# Xenium Prime 5K ovarian (10x release 3.0.0): the large-panel reference for
# the prime5k expression models.  In this zip the members sit at the top
# level (no ``outs/`` prefix, unlike the 1.x bundles).
PRIME5K_OUTS_ZIP = (
    "https://cf.10xgenomics.com/samples/xenium/3.0.0/"
    "Xenium_Prime_Ovarian_Cancer_FFPE_XRrun/"
    "Xenium_Prime_Ovarian_Cancer_FFPE_XRrun_outs.zip"
)
PRIME5K_MATRIX_MEMBER = "cell_feature_matrix.h5"
PRIME5K_CELL_SUBSAMPLE = 30_000   # cells clustered into the expression model

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


def ensure_prime5k() -> Path:
    """Fetch only ``cell_feature_matrix.h5`` (~109 MB) of the Xenium Prime 5K
    ovarian ``outs.zip`` (26.7 GB) with remotezip range requests (cached)."""
    from remotezip import RemoteZip

    dest = cache_dir() / "xenium_prime5k_ovarian_cell_feature_matrix.h5"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    cache_dir().mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".h5.part")
    with RemoteZip(PRIME5K_OUTS_ZIP) as z:
        info = next((i for i in z.infolist() if i.filename == PRIME5K_MATRIX_MEMBER),
                    None)
        if info is None:
            raise KeyError(f"member {PRIME5K_MATRIX_MEMBER!r} not in {PRIME5K_OUTS_ZIP}")
        with z.open(PRIME5K_MATRIX_MEMBER) as src, open(tmp, "wb") as out:
            while True:
                block = src.read(8 << 20)
                if not block:
                    break
                out.write(block)
    if tmp.stat().st_size != info.file_size:
        got, want = tmp.stat().st_size, info.file_size
        tmp.unlink()
        raise IOError(f"Prime 5K matrix download truncated: {got} != {want} bytes")
    os.replace(tmp, dest)
    return dest


def build_model_from_xenium_h5(h5_path: str | os.PathLike, *,
                               n_genes: int | None = None,
                               n_types: int = 15,
                               cell_subsample: int = PRIME5K_CELL_SUBSAMPLE,
                               seed: int = MODEL_SEED_DEFAULT,
                               min_cells_per_type: int = 50,
                               name: str = "realistic_xenium_h5_kmeans"
                               ):
    """Build a realistic expression model from a 10x ``cell_feature_matrix.h5``.

    Deterministic pipeline (``seed`` fixes both the cell subsample and the
    clustering):

    1. read the CSC matrix (features x cells), keep Gene Expression rows;
    2. subsample up to ``cell_subsample`` cells (sorted, so order-stable);
    3. optionally subset the panel to the top ``n_genes`` by mean expression
       over the subsample (ties broken by gene name) — a 1000-gene panel out
       of the Prime 5K, say;
    4. median-total normalisation + ``log1p``, PCA to 50 components and
       ``KMeans`` into ``n_types`` clusters (both seeded by ``seed``;
       clusters below ``min_cells_per_type`` dropped);
    5. per-cluster mean profiles -> rows-normalised composition, cluster
       frequencies -> type proportions, upstream ``_exclusivity`` for markers.

    This mirrors upstream ``build_realistic_model_from_xenium`` but works on a
    cell subsample and supports gene subsetting, so the 5000-gene matrix never
    has to be densified in full.  (MiniBatchKMeans — the upstream choice —
    collapses to singleton clusters on this matrix; PCA + Lloyd KMeans does
    not, and is faster.)
    """
    import h5py
    import scipy.sparse as sp
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA

    generator, expression, oracle, config = _import_strec()
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        shape = tuple(int(v) for v in m["shape"][:])
        M = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=shape)
        ftype = np.array([x.decode() for x in m["features"]["feature_type"][:]])
        names = np.array([x.decode() for x in m["features"]["name"][:]])
    gmask = ftype == "Gene Expression"
    M = M[gmask, :]                                # (genes, cells)
    genes = [str(g) for g in names[gmask]]
    if M.shape[0] == 0:
        raise ValueError(f"no Gene Expression features in {h5_path}")

    n_cells = M.shape[1]
    rng = np.random.default_rng(seed)
    take = np.sort(rng.choice(n_cells,
                              size=min(int(cell_subsample), n_cells),
                              replace=False))
    Ms = M[:, take]                                # (genes, cells_subsample)

    if n_genes is not None and n_genes < Ms.shape[0]:
        # top-N by mean expression over the subsample, ties broken by name
        means = np.asarray(Ms.sum(axis=1)).ravel().astype(np.float64)
        order = sorted(range(len(genes)), key=lambda i: (-means[i], genes[i]))
        keep = np.sort(np.asarray(order[:n_genes], dtype=np.int64))
        Ms = Ms[keep, :]
        genes = [genes[i] for i in keep]

    X = Ms.T.tocsr().toarray().astype(np.float32)  # cells x genes, dense subsample
    totals = X.sum(axis=1).astype(np.float64)
    sf = np.where(totals > 0, totals, 1.0).astype(np.float32)
    Xn = np.log1p(X / sf[:, None] * np.float32(np.median(totals)))
    n_components = int(min(50, Xn.shape[1], Xn.shape[0] - 1))
    Z = PCA(n_components=n_components, random_state=seed,
            svd_solver="randomized").fit_transform(Xn)
    km = KMeans(n_clusters=n_types, random_state=seed, n_init=10)
    raw = km.fit_predict(Z).astype(np.int32)

    keep_clusters = [t for t in range(n_types)
                     if int((raw == t).sum()) >= min_cells_per_type]
    if not keep_clusters:
        raise ValueError("every cluster smaller than min_cells_per_type")
    remap = {t: i for i, t in enumerate(keep_clusters)}
    keep_cell = np.isin(raw, list(remap))          # drop cells of tiny clusters
    X, raw = X[keep_cell], raw[keep_cell]
    labels = np.array([remap[t] for t in raw], dtype=np.int32)
    k = len(keep_clusters)

    mean_expr = np.vstack([X[labels == t].mean(axis=0) for t in range(k)])
    proportions = np.array([(labels == t).mean() for t in range(k)], dtype=float)
    row = mean_expr.sum(axis=1, keepdims=True)
    composition = np.divide(mean_expr, row, out=np.zeros_like(mean_expr),
                            where=row > 0)
    owner = expression._exclusivity(composition, 0.7)
    return expression.ExpressionModel(
        name=name, n_types=k, n_genes=len(genes),
        type_names=[f"k{t}" for t in range(k)],
        proportions=proportions, composition=composition, gene_names=genes,
        excl_threshold=0.7, excl_owner=owner, mean_expr=mean_expr,
    )


def build_prime5k_model(n_genes: int):
    """Prime 5K ovarian realistic model subset to ``n_genes`` (1000 or 5000)."""
    h5 = ensure_prime5k()
    return build_model_from_xenium_h5(
        h5, n_genes=n_genes, name=f"realistic_prime5k_k{n_genes}")


def build_model(kind: str):
    """Expression model: ``disjoint`` (toy), ``merfish``, ``xenium`` or a
    ``prime5k<N>`` large panel from Xenium Prime 5K (N = 1000, 5000)."""
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
    if kind.startswith("prime5k"):
        n_genes = int(kind[len("prime5k"):])
        if not 100 <= n_genes <= 10_000:
            raise ValueError(f"prime5k panel size must be 100..10000 genes, "
                             f"got {n_genes}")
        model = build_prime5k_model(n_genes)
        return model, {
            "reference": "Xenium_Prime_Ovarian_Cancer_FFPE_XRrun cell_feature_matrix.h5 "
                         "(10x release 3.0.0)",
            "url": PRIME5K_OUTS_ZIP,
            "member": PRIME5K_MATRIX_MEMBER,
            "n_types": 15,
            "cell_subsample": PRIME5K_CELL_SUBSAMPLE,
            "model_seed": MODEL_SEED_DEFAULT,
            "gene_selection": "top-N by mean expression over the subsample "
                              "(ties by gene name)",
        }
    raise KeyError(f"unknown model {kind!r}; have disjoint, merfish, xenium, "
                   "prime5k1000, prime5k5000")


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


def _naive_assign_3d(obs3: np.ndarray, centers3: np.ndarray) -> np.ndarray:
    """Nearest nucleus in 3D (x, y, z) — the naive baseline for the z variant."""
    from scipy.spatial import cKDTree
    _, nn = cKDTree(centers3).query(obs3, k=1)
    return nn.astype(np.int32)


def generate(*, dataset_id: str, tier: str, packing: float, sigma: float,
             model_kind: str, seed: int, mean_tx: float = MEAN_TX_DEFAULT,
             n_target: int = N_TARGET_DEFAULT, geometry: str = "voronoi",
             bg_frac: float = 0.0, z_slab_um: float | None = None,
             prior_opts: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """Build one st-recoverability field and map it to the dataset contract.

    Optional extensions over the upstream demo (defaults reproduce the
    original datasets byte-for-byte):

    ``geometry``
        ``"voronoi"`` (upstream default) or ``"aniso"`` — the anisotropic
        label mode (elongated, irregular cells) of the upstream generator.
    ``bg_frac``
        share of ambient molecules added by this wrapper uniformly over the
        domain with ``cell = 0`` (genes from the tissue-average profile); the
        upstream generator only emits cell transcripts.
    ``z_slab_um``
        3D variant: every cell gets a nucleus centre ``z0 ~ U(0, slab)``, its
        molecules a uniform ``z`` within half the cell thickness around it,
        and the observed ``z`` the same Gaussian displacement as x/y.  The
        wrapper adds the ``z`` column; the upstream generator is 2D.
    ``prior_opts``
        ``kind: imperfect`` nucleus segmentation (missed / shifted / merged
        nuclei, see :func:`common.imperfect_nucleus_prior`); truth unchanged.
    """
    generator, expression, oracle, config = _import_strec()
    model, model_info = build_model(model_kind)

    seed_field, seed_tx = int(seed), int(seed) + 1
    field = generator.build_field(packing, sigma, seed_field, n_target=n_target,
                                  model=model, geometry=geometry)
    tx = generator.generate_transcripts(field, mean_tx, seed_tx)

    # units check + oracle / naive accuracies (interior molecules, as in demo.py)
    units = _check_units(field)
    dmax, argcell = oracle.build_oracle_maps(field)
    # upstream oracle takes log(p[gene|type]); zeros are guarded by `pos`
    with np.errstate(divide="ignore"):
        oracle_assign = oracle.oracle_assign(field, dmax, argcell, tx.obs_xy, tx.gene)

    obs_xy = tx.obs_xy
    z_obs = None
    prior_obs = obs_xy
    prior_centers = field.centers
    if z_slab_um is not None:
        # wrapper-added third dimension (the upstream generator is 2D)
        slab = float(z_slab_um)
        rng_z = np.random.default_rng(int(seed) + 3)
        z0 = rng_z.uniform(0.0, slab, size=field.n_cells)
        half = min(field.r_mean_um, slab / 2.0)
        z_true = np.clip(z0[tx.true_cell]
                         + rng_z.uniform(-half, half, size=len(tx.true_cell)),
                         0.0, slab)
        z_obs = z_true + rng_z.normal(0.0, field.sigma_um,
                                      size=len(tx.true_cell))
        prior_obs = np.column_stack([obs_xy, z_obs])
        prior_centers = np.column_stack([field.centers, z0])
        naive_assign = _naive_assign_3d(prior_obs, prior_centers)
        # The upstream oracle works on a 2D pixel grid and cannot use z, while
        # Baysor sees z here; do not claim a ceiling the method could exceed.
        oracle_acc = None
    else:
        naive_assign = oracle.naive_assign(field, obs_xy)
        oracle_acc = _accuracy(tx.true_cell, oracle_assign, tx.interior)
    naive_acc = _accuracy(tx.true_cell, naive_assign, tx.interior)

    # ambient background molecules (cell = 0), added by this wrapper
    n_bg = 0
    bg_xy = bg_gene = bg_z = None
    if bg_frac > 0.0:
        rng_bg = np.random.default_rng(int(seed) + 2)
        n_bg = int(round(bg_frac / (1.0 - bg_frac) * len(tx.true_cell)))
        bg_xy = rng_bg.uniform(0.0, field.L_um, size=(n_bg, 2))
        pooled = (np.asarray(model.proportions, dtype=np.float64)
                  @ np.asarray(model.composition, dtype=np.float64))
        pooled = pooled / pooled.sum()
        bg_gene = rng_bg.choice(model.n_genes, size=n_bg, p=pooled)
        if z_slab_um is not None:
            bg_z = rng_bg.uniform(0.0, float(z_slab_um), size=n_bg)

    # nucleus prior over every molecule: within 3 um of the true centre
    # (methods_baysor.py logic); a background molecule inside a nucleus keeps
    # a prior label, as in real data
    if n_bg:
        all_obs = (np.column_stack([np.vstack([obs_xy, bg_xy]),
                                    np.concatenate([z_obs, bg_z])])
                   if z_slab_um is not None else np.vstack([obs_xy, bg_xy]))
    else:
        all_obs = prior_obs
    prior_info = None
    if prior_opts:
        if prior_opts.get("kind", "imperfect") != "imperfect":
            raise ValueError(f"unknown prior_opts kind {prior_opts.get('kind')!r}")
        if "seed" not in prior_opts:
            raise ValueError("prior_opts requires an explicit seed")
        prior, prior_info = common.imperfect_nucleus_prior(
            all_obs, prior_centers, NUCLEUS_RADIUS_UM,
            seed=int(prior_opts["seed"]),
            miss_frac=float(prior_opts.get("miss_frac", 0.2)),
            shift_um=tuple(prior_opts.get("shift_um", (1.0, 2.0))),
            merge_frac=float(prior_opts.get("merge_frac", 0.05)))
        if prior_opts.get("base"):
            prior_info["base"] = prior_opts["base"]
    else:
        prior, _ = common.nucleus_prior(all_obs, prior_centers,
                                        NUCLEUS_RADIUS_UM)

    type_names = [str(t) for t in model.type_names]
    gene_names = np.asarray(model.gene_names, dtype=object)
    x = np.concatenate([obs_xy[:, 0], bg_xy[:, 0]]) if n_bg else obs_xy[:, 0]
    y = np.concatenate([obs_xy[:, 1], bg_xy[:, 1]]) if n_bg else obs_xy[:, 1]
    gene = (np.concatenate([gene_names[tx.gene], gene_names[bg_gene]])
            if n_bg else gene_names[tx.gene])
    cell = (np.concatenate([(tx.true_cell + 1).astype(np.int32),
                            np.zeros(n_bg, dtype=np.int32)])
            if n_bg else (tx.true_cell + 1).astype(np.int32))
    interior = tx.interior.astype(bool)
    if n_bg:
        margin = field.r_mean_um + getattr(generator, "INTERIOR_MARGIN_SIGMA",
                                           3.0) * field.sigma_um
        bg_interior = ((bg_xy[:, 0] >= margin) & (bg_xy[:, 0] <= field.L_um - margin)
                       & (bg_xy[:, 1] >= margin) & (bg_xy[:, 1] <= field.L_um - margin))
        interior = np.concatenate([interior, bg_interior]).astype(bool)
    celltype = (np.concatenate([np.asarray(type_names, dtype=object)[tx.true_type],
                                np.full(n_bg, "", dtype=object)])
                if n_bg else np.asarray(type_names, dtype=object)[tx.true_type])
    df = pd.DataFrame({
        "x": x, "y": y, "gene": gene,
        "prior": prior, "cell": cell, "interior": interior,
        "celltype": celltype,
    })
    if z_obs is not None:
        z = (np.concatenate([z_obs, bg_z]) if n_bg else z_obs)
        df.insert(2, "z", z.astype(np.float64))
    df = common.make_molecules(df)

    area_um2 = field.L_um ** 2
    stats = common.make_stats(df, area_um2=area_um2, n_true_cells=field.n_cells)
    if z_slab_um is not None:
        stats["volume_um3"] = float(area_um2 * float(z_slab_um))
    n_genes = int(model.n_genes)
    notes = (f"st-recoverability field: packing {packing:g} cells/mm2, "
             f"sigma {sigma:g} um, model {model.name}")
    if geometry != "voronoi":
        notes += f", {geometry} geometry (elongated, irregular cells)"
    if bg_frac > 0.0:
        notes += f", {bg_frac:.0%} ambient background (cell = 0)"
    if z_slab_um is not None:
        notes += f", 3D slab {z_slab_um:g} um (z added by the wrapper)"
    difficulty = {
        "cell_density": common.density_class(stats["true_cells_per_mm2"]),
        "gene_panel": common.gene_panel_class(n_genes),
        "notes": notes,
    }
    params = {
        "packing_cells_per_mm2": packing,
        "sigma_um": sigma,
        "model": model_kind,
        "model_name": model.name,
        "model_info": model_info,
        "mean_tx_per_cell": mean_tx,
        "n_target": n_target,
        "geometry": geometry,
        "displacement": "gaussian",
        "nucleus_radius_um": NUCLEUS_RADIUS_UM,
    }
    if bg_frac > 0.0:
        params["ambient_bg_frac"] = float(bg_frac)
    if z_slab_um is not None:
        params.update({"z_slab_um": float(z_slab_um), "z_sigma_um": float(sigma),
                       "naive_assignment": "nearest nucleus in 3D"})
    accuracy_subset = "interior molecules (tx.interior)"
    if n_bg:
        accuracy_subset += "; ambient background (cell = 0) excluded"
    truth_extra = {}
    if prior_info is not None:
        truth_extra["prior"] = prior_info
    if oracle_acc is None:
        truth_extra["oracle_accuracy_note"] = (
            "not computed: the upstream oracle is a 2D pixel grid and cannot "
            "use z, while this dataset exposes z to the method; the 3D naive "
            "baseline is recorded instead")
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
            "z_range_um": ([0.0, float(z_slab_um)] if z_slab_um is not None
                           else None),
            "note": ("simulated field (2D)" if z_slab_um is None else
                     f"simulated field (2D layout) with wrapper-added z in a "
                     f"{z_slab_um:g} um slab"),
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
            "accuracy_subset": accuracy_subset,
            **truth_extra,
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
                   help="cell packing, cells/mm^2 (e.g. 1000, 2575, 6000, 13625)")
    p.add_argument("--sigma", type=float, required=True,
                   help="displacement sigma in um (e.g. 1.0, 2.0, 3.0)")
    p.add_argument("--model", required=True,
                   choices=["disjoint", "merfish", "xenium",
                            "prime5k1000", "prime5k5000"])
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--mean-tx", type=float, default=MEAN_TX_DEFAULT)
    p.add_argument("--n-target", type=int, default=N_TARGET_DEFAULT)
    p.add_argument("--geometry", choices=["voronoi", "aniso"], default="voronoi",
                   help="cell geometry: upstream Voronoi or anisotropic (elongated)")
    p.add_argument("--bg-frac", type=float, default=0.0,
                   help="share of ambient background molecules (cell = 0)")
    p.add_argument("--z-slab-um", type=float, default=None,
                   help="add a wrapper z dimension: cells in a slab of this thickness")
    p.add_argument("--prior-opts", default=None,
                   help="JSON object for an imperfect prior, e.g. "
                        "'{\"kind\": \"imperfect\", \"seed\": 7301}'")
    p.add_argument("--tier", choices=["quick", "full"], default="quick")
    p.add_argument("-o", "--out-dir", required=True)
    args = p.parse_args(argv)

    prior_opts = json.loads(args.prior_opts) if args.prior_opts else None
    df, meta = generate(dataset_id=args.id, tier=args.tier, packing=args.packing,
                        sigma=args.sigma, model_kind=args.model, seed=args.seed,
                        mean_tx=args.mean_tx, n_target=args.n_target,
                        geometry=args.geometry, bg_frac=args.bg_frac,
                        z_slab_um=args.z_slab_um, prior_opts=prior_opts)
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
