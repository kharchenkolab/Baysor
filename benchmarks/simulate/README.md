# `benchmarks/simulate` — simulated datasets with exact ground truth

Generators for the `sim` dataset group of the Baysor benchmark suite.
Contract: [`benchmarks/README.md`](../README.md).  Manifest:
[`benchmarks/datasets/sim.yaml`](../datasets/sim.yaml).

* `common.py` — shared helpers: hex lattices, Delaunay adjacency + DSATUR
  graph colouring, expression profiles, molecule placement, nucleus prior
  (perfect and *imperfect*: missed / shifted / merged nuclei), contract
  assembly/IO, tier budgets.
* `trivial.py` — seven in-repo scenarios (`circles_gaps`, `tiled_distinct`,
  `tiled_same`, `mixed_sizes`, `sparse_noisy`, `circles_gaps_3d`,
  `elongated_gaps`).
* `strec.py` — wrapper around the external
  [st-recoverability](https://github.com/RuiYamasaki/st-recoverability)
  generator (Yamasaki, 2026; MIT), pinned commit
  `e85faeccfb58e0bef05d76817e98c79abff6bf51`, imported through `sys.path`
  from a cache clone — never vendored into this repository.  Wrapper
  extensions: ambient background molecules, anisotropic (elongated) cell
  geometry, a wrapper-added `z` dimension, imperfect priors and the large
  Prime 5K expression models.
* `generate_all.py` — regenerates every manifest entry deterministically;
  `--verify` / `--verify-all` check the datasets on disk **and** a fresh
  regeneration against the `sha256` blocks committed in `sim.yaml`;
  `--update-hashes` records them.
* `sanity.py` — runs the release Baysor binary on chosen datasets and records
  wall time, a majority-matching and a one-to-one (Hungarian) assignment
  accuracy next to the oracle (reports: [`sanity_check.json`](sanity_check.json),
  [`sanity_check_noprior.json`](sanity_check_noprior.json)).
* `tests/` — pytest suite, 76 tests
  (`python -m pytest benchmarks/simulate/tests`).

## Datasets (51 total, 5.81 M molecules, ~126 MB parquet)

| dataset | tier | molecules | genes | cells | density | panel | kB | prior | oracle | naive |
|---|---|---|---|---|---|---|---|---|---|---|
| sim_circles_gaps_3d_g100 | quick | 93765 | 100 | 621 | medium | small | 3028 | column | — | — |
| sim_circles_gaps_3d_g100_noprior | quick | 93765 | 100 | 621 | medium | small | 2985 | none | — | — |
| sim_circles_gaps_g100 | quick | 92296 | 100 | 621 | medium | small | 2091 | column | — | — |
| sim_circles_gaps_g1000 | quick | 92296 | 1000 | 621 | medium | large | 2129 | column | — | — |
| sim_circles_gaps_g1000_noprior | quick | 92296 | 1000 | 621 | medium | large | 2071 | none | — | — |
| sim_circles_gaps_g100_imprior | quick | 92296 | 100 | 621 | medium | small | 2083 | imperfect | — | — |
| sim_circles_gaps_g100_noprior | quick | 92296 | 100 | 621 | medium | small | 2034 | none | — | — |
| sim_circles_gaps_g5000 | quick | 92296 | 4963 | 621 | medium | huge | 2178 | column | — | — |
| sim_circles_gaps_g5000_noprior | quick | 92296 | 4963 | 621 | medium | huge | 2121 | none | — | — |
| sim_elongated_gaps_g100 | quick | 37962 | 100 | 255 | sparse | small | 846 | column | — | — |
| sim_elongated_gaps_g100_noprior | quick | 37962 | 100 | 255 | sparse | small | 823 | none | — | — |
| sim_mixed_sizes_g100 | quick | 94564 | 100 | 672 | medium | small | 2120 | column | — | — |
| sim_mixed_sizes_g100_imprior | quick | 94564 | 100 | 672 | medium | small | 2109 | imperfect | — | — |
| sim_mixed_sizes_g100_noprior | quick | 94564 | 100 | 672 | medium | small | 2069 | none | — | — |
| sim_sparse_noisy_g100 | quick | 48399 | 100 | 247 | sparse | small | 1059 | column | — | — |
| sim_sparse_noisy_g100_noprior | quick | 48399 | 100 | 247 | sparse | small | 1043 | none | — | — |
| sim_tiled_distinct_g100 | quick | 115904 | 99 | 780 | medium | small | 2624 | column | — | — |
| sim_tiled_distinct_g1000 | quick | 115904 | 974 | 780 | medium | large | 2670 | column | — | — |
| sim_tiled_distinct_g1000_full | full | 972976 | 999 | 6525 | medium | large | 19329 | column | — | — |
| sim_tiled_distinct_g1000_noprior | quick | 115904 | 974 | 780 | medium | large | 2612 | none | — | — |
| sim_tiled_distinct_g100_noprior | quick | 115904 | 99 | 780 | medium | small | 2566 | none | — | — |
| sim_tiled_distinct_g5000 | quick | 115904 | 4557 | 780 | medium | huge | 2726 | column | — | — |
| sim_tiled_distinct_g5000_noprior | quick | 115904 | 4557 | 780 | medium | huge | 2669 | none | — | — |
| sim_tiled_same_g100 | quick | 115904 | 100 | 780 | medium | small | 2585 | column | — | — |
| sim_tiled_same_g100_noprior | quick | 115904 | 100 | 780 | medium | small | 2527 | none | — | — |
| strec_dense_s1_disjoint | quick | 63999 | 20 | 400 | dense | tiny | 1434 | column | 0.853 | 0.808 |
| strec_dense_s1_disjoint_noprior | quick | 63999 | 20 | 400 | dense | tiny | 1381 | none | 0.853 | 0.808 |
| strec_dense_s2_disjoint | quick | 64441 | 20 | 400 | dense | tiny | 1440 | column | 0.723 | 0.648 |
| strec_dense_s2_disjoint_3d | quick | 64267 | 20 | 400 | dense | tiny | 2050 | column | — | 0.671 |
| strec_dense_s2_disjoint_ambient5 | quick | 66823 | 20 | 400 | dense | tiny | 1504 | column | 0.723 | 0.636 |
| strec_dense_s2_disjoint_noprior | quick | 64441 | 20 | 400 | dense | tiny | 1392 | none | 0.723 | 0.648 |
| strec_dense_s2_merfish | quick | 63889 | 155 | 400 | dense | small | 1470 | column | 0.723 | 0.635 |
| strec_dense_s2_merfish_ambient20 | quick | 80195 | 155 | 400 | dense | small | 1844 | column | 0.733 | 0.641 |
| strec_dense_s2_merfish_aniso | quick | 64352 | 155 | 400 | dense | small | 1478 | column | 0.710 | 0.575 |
| strec_dense_s2_merfish_imprior | quick | 63889 | 155 | 400 | dense | small | 1462 | imperfect | 0.723 | 0.635 |
| strec_dense_s2_merfish_noprior | quick | 63889 | 155 | 400 | dense | small | 1419 | none | 0.723 | 0.635 |
| strec_dense_s2_prime5k1000 | quick | 64027 | 1000 | 400 | dense | large | 1496 | column | 0.698 | 0.640 |
| strec_dense_s2_prime5k5000 | quick | 64098 | 4424 | 400 | dense | huge | 1539 | column | 0.701 | 0.644 |
| strec_dense_s2_xenium | quick | 64285 | 311 | 400 | dense | medium | 1489 | column | 0.737 | 0.642 |
| strec_dense_s2_xenium_full | full | 999729 | 313 | 2500 | dense | medium | 20031 | column | 0.735 | 0.642 |
| strec_dense_s3_disjoint | quick | 63626 | 20 | 400 | dense | tiny | 1425 | column | 0.602 | 0.496 |
| strec_medium_s2_disjoint | quick | 63982 | 20 | 400 | medium | tiny | 1414 | column | 0.812 | 0.753 |
| strec_medium_s3_xenium | quick | 64211 | 312 | 400 | medium | medium | 1470 | column | 0.740 | 0.644 |
| strec_medium_s3_xenium_imprior | quick | 64211 | 312 | 400 | medium | medium | 1465 | imperfect | 0.740 | 0.644 |
| strec_medium_s3_xenium_noprior | quick | 64211 | 312 | 400 | medium | medium | 1436 | none | 0.740 | 0.644 |
| strec_sparse_s1_disjoint | quick | 63743 | 20 | 400 | medium | tiny | 1397 | column | 0.926 | 0.904 |
| strec_sparse_s1_disjoint_noprior | quick | 63743 | 20 | 400 | medium | tiny | 1376 | none | 0.926 | 0.904 |
| strec_sparse_s2_disjoint | quick | 64191 | 20 | 400 | medium | tiny | 1406 | column | 0.870 | 0.832 |
| strec_sparse_s2_merfish | quick | 63582 | 155 | 400 | medium | small | 1432 | column | 0.873 | 0.831 |
| strec_sparse_s2_merfish_noprior | quick | 63582 | 155 | 400 | medium | small | 1412 | none | 0.873 | 0.831 |
| strec_ultrasparse_s2_disjoint | quick | 64157 | 20 | 400 | sparse | tiny | 1395 | column | 0.913 | 0.888 |

`genes` is the number of genes *observed* in `molecules.parquet` (the panel
size is recorded in `meta.truth.params.n_genes` / `model_info`; a handful of
low-expression background genes of a large panel can be unobserved — the
Prime 5K model subsets to 5000 of the 5101 Gene Expression features, 4424 of
them fire in the field).  `oracle`/`naive` are the st-recoverability
Bayes-optimal / nearest-nucleus accuracies on interior molecules
(`meta.truth`); trivial datasets have exact truth by construction (no
displacement), so no oracle is computed.  `prior` is the value of
`meta.baysor.prior`: `column` (perfect 3 µm nucleus prior), `none` (the
noprior variants: identical molecules, no prior column) or `imperfect`
(degraded prior, truth unchanged).

Notes:

* the st-recoverability "sparse" regime (2575 cells/mm²) falls into the
  contract's `medium` density class (sparse is < 2500); only
  `strec_ultrasparse_s2_disjoint` (1000 cells/mm²) is truly sparse;
* `strec_*_noprior` share their base's molecules (only the prior column is
  dropped); `<base>_imprior` share the base's geometry and truth;
* every new dataset is quick tier (≤ 150k molecules); the two full-tier
  datasets are `sim_tiled_distinct_g1000_full` (0.97 M) and
  `strec_dense_s2_xenium_full` (~1.00 M).

## Regenerating

```bash
export BAYSOR_BENCH_DATA=/path/to/bench-data      # default: <repo>/.bench-data
PY=.deps/bench/bin/python

$PY benchmarks/simulate/generate_all.py                    # all 51 datasets (~90 s)
$PY benchmarks/simulate/generate_all.py --list             # show the manifest
$PY benchmarks/simulate/generate_all.py --only sim_circles_gaps_g100
$PY benchmarks/simulate/generate_all.py --verify <id>      # disk + regeneration
                                                           # vs manifest sha256
$PY benchmarks/simulate/generate_all.py --verify-all       # every dataset
$PY benchmarks/simulate/generate_all.py --update-hashes    # record sha256
```

Every manifest entry carries the expected `sha256` of `molecules.parquet`
and `meta.json`; `--verify-all` fails if either the files on disk or a fresh
regeneration differ from the manifest, so a fresh machine can prove it
reproduces the committed bytes.

Single datasets can also be built directly:

```bash
$PY benchmarks/simulate/trivial.py circles_gaps --id sim_circles_gaps_g100 \
    --seed 1201 --n-genes 100 -o "$BAYSOR_BENCH_DATA/sim/sim_circles_gaps_g100"
$PY benchmarks/simulate/strec.py --id strec_dense_s2_merfish \
    --packing 13625 --sigma 2.0 --model merfish --seed 910009 \
    -o "$BAYSOR_BENCH_DATA/sim/strec_dense_s2_merfish"
$PY benchmarks/simulate/strec.py --id strec_dense_s2_disjoint_aniso \
    --packing 13625 --sigma 2.0 --model disjoint --seed 910015 --geometry aniso \
    -o "$BAYSOR_BENCH_DATA/sim/strec_dense_s2_disjoint_aniso"
```

Everything is seeded: re-running `generate_all.py` reproduces
`molecules.parquet` and `meta.json` byte-for-byte (verified by `--verify`,
also covered by `tests/test_manifest.py`).

## Design notes

* **Trivial scenarios** place molecules uniformly inside round cells or
  Voronoi regions (no displacement); cell types get a few marker genes over a
  shared lognormal background; `tiled_distinct` uses strictly disjoint gene
  blocks per DSATUR colour of the Delaunay (= Voronoi) adjacency, so neighbour
  cells always differ in type.  Molecules per cell are lognormal in [50, 300]
  (medians vary by scenario; `mixed_sizes` uses 70/190 for small/large).
  Background molecules are uniform with `cell = 0`.
* **Panel ablation**: the `g100/g1000/g5000` variants and the tiled pair share
  seed 1201, hence identical geometry, cell ids and molecule counts — only the
  gene panel changes.  `tiled_same` is the exact geometry of
  `tiled_distinct` with one type (composition carries no boundary signal).
* **`elongated_gaps`** (seed 1205): fibroblast-like (8 × 3 µm) and
  neuron-like (11 × 2.4 µm) ellipses with random orientation, conservative
  gaps (centre pairs ≥ 2 · max semi-major + 2 µm), round 3 µm nucleus prior;
  both a with-prior and a noprior dataset exist.
* **Noprior variants** (`generator: noprior` in the manifest) rebuild their
  base and drop the `prior` column: `meta.baysor.prior = "none"`, every other
  byte of the molecule table identical to the base.  They cover all eleven
  quick trivial datasets and six st-recoverability datasets (sparse/medium/
  dense × σ ∈ {1, 2, 3} × disjoint/merfish/xenium).  The full-tier
  `sim_tiled_distinct_g1000_full` has no noprior copy because new datasets
  must stay in the quick tier; the quick `g1000` noprior covers the scenario.
* **Imperfect-prior variants** (`<id>_imprior`, seeds 7301–7304): 20% of
  nuclei missed (no prior anywhere in those cells), every remaining nucleus
  shifted by a uniformly random 1–2 µm, and 5% of the nuclei merged onto
  their nearest unmerged neighbour (one prior label covers two nuclei).
  Implemented by `common.imperfect_nucleus_prior` with an RNG independent of
  the molecule streams, so the truth columns are byte-identical to the base;
  the applied defects are recorded in `meta.truth.prior`.
* **`meta.baysor`**: `scale_um` = the true mean cell radius (disc radius for
  round cells, `sqrt(area/cell / pi)` otherwise); `min_molecules_per_cell` 20
  and `prior_confidence` 0.5 follow the contract example / Xenium config.
  `config` is `null` for simulated datasets: all parameters are passed as
  explicit CLI flags by the runner.
* **st-recoverability datasets** are a design over packing ∈
  {1000, 2575, 6000, 13625} cells/mm², σ ∈ {1, 2, 3} µm and model ∈
  {disjoint, MERFISH-realistic, Xenium-realistic, Prime 5K 1000/5000} —
  each factor covered, the required dense + σ=2 + realistic cell present, not
  the full cartesian product.  Oracle accuracy falls from 0.93 (sparse, σ=1)
  to 0.60 (dense, σ=3), reproducing the upstream frontier (their dense σ=2
  demo value ≈ 0.72).
* **Ambient noise**: `strec_dense_s2_disjoint_ambient5` /
  `strec_dense_s2_merfish_ambient20` add 5% / 20% uniformly distributed
  background molecules (`cell = 0`, genes from the tissue-average profile) in
  the wrapper; the oracle/naive accuracies are computed over true-cell
  interior molecules only (`meta.truth.accuracy_subset` says so).
* **Anisotropic geometry**: `strec_dense_s2_merfish_aniso` uses the upstream
  `geometry="aniso"` label mode (`_aniso_label`, elongated randomly oriented
  Mahalanobis cells) — irregular, non-Voronoi boundaries with the oracle
  recomputed on the anisotropic label image.
* **3D**: `strec_dense_s2_disjoint_3d` gets a `z` column from the wrapper
  (nucleus z₀ ~ U(0, 6 µm), molecule z uniform within the cell thickness,
  observed z displaced by the same σ as x/y; prior computed in 3D).  The
  upstream oracle is 2D-only and is therefore *not* claimed (`oracle_accuracy`
  is null with a note); `naive_accuracy` is nearest-nucleus in 3D.
* **Large panels**: `prime5k1000` / `prime5k5000` expression models are
  built by `build_model_from_xenium_h5` from the Xenium Prime 5K ovarian
  `cell_feature_matrix.h5`: 30 k-cell subsample, top-N genes by mean
  expression (ties by name), median-normalised log1p, PCA(50) + seeded
  KMeans into 15 types, per-cluster mean profiles.  (Upstream
  MiniBatchKMeans collapses to singleton clusters on this matrix; PCA +
  Lloyd does not.)  Both prime5k datasets share seed 910014, hence the same
  field geometry — a panel-size ablation at fixed field.

## External cache (`$BAYSOR_BENCH_DATA/cache/sim`, ~135 MB total)

| file | source | size |
|---|---|---|
| `st-recoverability/` | github clone at pinned commit (only `src/{generator,expression,oracle,config}.py` are imported) | 13 MB |
| `merfish_moffitt.h5ad` | squidpy `MERFISH_0.24.h5ad` mirror, figshare file 40038538 (same URL as upstream `realism.py`) | 3.7 MB |
| `xenium_breast_rep1_cell_feature_matrix.h5` | only this member of the 10x `Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip` (9.9 GB), fetched with `remotezip` range requests | 12 MB |
| `xenium_breast_rep1_cells.parquet` | likewise, from the same zip | 3.3 MB |
| `xenium_prime5k_ovarian_cell_feature_matrix.h5` | only this member (top-level `cell_feature_matrix.h5`) of the 10x `Xenium_Prime_Ovarian_Cancer_FFPE_XRrun_outs.zip` (26.7 GB), fetched with `remotezip` range requests | 104 MB |

All downloads are re-fetchable by the scripts (`ensure_repo`, `ensure_merfish`,
`ensure_xenium`, `ensure_prime5k`); nothing external is committed.

## Sanity check (BENCH-SIM part C)

`python benchmarks/simulate/sanity.py` runs the release binary on one trivial,
one sparse and one dense st-rec dataset with `meta.baysor` + `:prior`
(≤ 6 threads), writes run outputs to
`$BAYSOR_BENCH_DATA/runs/bench-sim-sanity/<id>/` and the report to
`sanity_check.json`.  Reports record the majority-match accuracy and a
one-to-one (Hungarian) matched accuracy — the latter is computed by
`sanity.hungarian_match_accuracy`, a helper for this check only (the
benchmark metric lives in `harness/metrics.py`).  Latest with-prior run:

| dataset | exit | wall time | majority-match | one-to-one | oracle | naive |
|---|---|---|---|---|---|---|
| sim_circles_gaps_g100 | 0 | 30.1 s | **0.987** | — | ≈ 0.99 (1 − background) | — |
| strec_sparse_s2_disjoint | 0 | 26.9 s | **0.628** | — | 0.870 | 0.832 |
| strec_dense_s2_merfish | 0 | 31.5 s | **0.619** | — | 0.723 | 0.635 |

No-prior sanity run
(`sanity.py --ids ..._noprior --report sanity_check_noprior.json`, same
binary, no `:prior` flag):

| dataset | exit | wall time | majority-match | one-to-one | oracle | naive |
|---|---|---|---|---|---|---|
| strec_sparse_s1_disjoint_noprior | 0 | 23.6 s | 0.638 | **0.633** | 0.926 | 0.904 |
| strec_dense_s2_merfish_noprior | 0 | 28.4 s | 0.517 | **0.501** | 0.723 | 0.635 |
| sim_circles_gaps_g100_noprior | 0 | 26.3 s | 0.987 | **0.987** | ≈ 0.99 | — |

Without a prior, Baysor keeps the near-perfect score on the gap-separated
trivial field (the geometry alone resolves it) but drops ~0.09–0.22 on the
displaced st-recoverability fields — the no-prior mode is measurably weaker
than the prior mode and sits between the oracle and (on the sparse field)
below the nearest-nucleus baseline, which is what the harness now measures
properly.  Full metrics are BENCH-HARNESS's job; these numbers are a
plausibility check only.
