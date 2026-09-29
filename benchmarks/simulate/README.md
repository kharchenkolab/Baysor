# `benchmarks/simulate` — simulated datasets with exact ground truth

Generators for the `sim` dataset group of the Baysor benchmark suite.
Contract: [`benchmarks/README.md`](../README.md).  Manifest:
[`benchmarks/datasets/sim.yaml`](../datasets/sim.yaml).

* `common.py` — shared helpers: hex lattices, Delaunay adjacency + DSATUR
  graph colouring, expression profiles, molecule placement, nucleus prior,
  contract assembly/IO, tier budgets.
* `trivial.py` — six in-repo scenarios (`circles_gaps`, `tiled_distinct`,
  `tiled_same`, `mixed_sizes`, `sparse_noisy`, `circles_gaps_3d`).
* `strec.py` — wrapper around the external
  [st-recoverability](https://github.com/RuiYamasaki/st-recoverability)
  generator (Yamasaki, 2026; MIT), pinned commit
  `e85faeccfb58e0bef05d76817e98c79abff6bf51`, imported through `sys.path`
  from a cache clone — never vendored into this repository.
* `generate_all.py` — regenerates every manifest entry deterministically;
  `--verify` proves byte-identity by re-hashing.
* `sanity.py` — runs the release Baysor binary on three datasets and records
  wall time + a majority-matching assignment accuracy next to the oracle
  (report: [`sanity_check.json`](sanity_check.json)).
* `tests/` — pytest suite (`python -m pytest benchmarks/simulate/tests`).

## Datasets (21 total, 2.59 M molecules, ~58 MB parquet)

| dataset                        | tier  | molecules | genes | cells | density  | panel  | kB   | oracle | naive |
|--------------------------------|-------|-----------|-------|-------|----------|--------|------|--------|-------|
| sim_circles_gaps_g100          | quick |     92296 |   100 |   621 | medium   | small  | 2091 | —      | —     |
| sim_circles_gaps_g1000         | quick |     92296 |  1000 |   621 | medium   | large  | 2129 | —      | —     |
| sim_circles_gaps_g5000         | quick |     92296 |  4963 |   621 | medium   | huge   | 2178 | —      | —     |
| sim_tiled_distinct_g100        | quick |    115904 |    99 |   780 | medium   | small  | 2624 | —      | —     |
| sim_tiled_distinct_g1000       | quick |    115904 |   974 |   780 | medium   | large  | 2670 | —      | —     |
| sim_tiled_distinct_g5000       | quick |    115904 |  4557 |   780 | medium   | huge   | 2726 | —      | —     |
| sim_tiled_distinct_g1000_full  | full  |    972976 |   999 |  6525 | medium   | large  | 19329| —      | —     |
| sim_tiled_same_g100            | quick |    115904 |   100 |   780 | medium   | small  | 2585 | —      | —     |
| sim_mixed_sizes_g100           | quick |     94564 |   100 |   672 | medium   | small  | 2120 | —      | —     |
| sim_sparse_noisy_g100          | quick |     48399 |   100 |   247 | sparse   | small  | 1059 | —      | —     |
| sim_circles_gaps_3d_g100       | quick |     93765 |   100 |   621 | medium   | small  | 3028 | —      | —     |
| strec_sparse_s1_disjoint       | quick |     63743 |    20 |   400 | medium   | tiny   | 1397 | 0.926  | 0.904 |
| strec_sparse_s2_disjoint       | quick |     64191 |    20 |   400 | medium   | tiny   | 1406 | 0.870  | 0.832 |
| strec_sparse_s2_merfish        | quick |     63582 |   155 |   400 | medium   | small  | 1432 | 0.873  | 0.831 |
| strec_medium_s2_disjoint       | quick |     63982 |    20 |   400 | medium   | tiny   | 1414 | 0.812  | 0.753 |
| strec_medium_s3_xenium         | quick |     64211 |   312 |   400 | medium   | medium | 1470 | 0.740  | 0.644 |
| strec_dense_s1_disjoint        | quick |     63999 |    20 |   400 | dense    | tiny   | 1434 | 0.853  | 0.808 |
| strec_dense_s2_disjoint        | quick |     64441 |    20 |   400 | dense    | tiny   | 1440 | 0.723  | 0.648 |
| strec_dense_s3_disjoint        | quick |     63626 |    20 |   400 | dense    | tiny   | 1425 | 0.602  | 0.496 |
| strec_dense_s2_merfish         | quick |     63889 |   155 |   400 | dense    | small  | 1470 | 0.723  | 0.635 |
| strec_dense_s2_xenium          | quick |     64285 |   311 |   400 | dense    | medium | 1489 | 0.737  | 0.642 |

`genes` is the number of genes *observed* in `molecules.parquet` (the panel
size is recorded in `meta.truth.params.n_genes`; a handful of low-expression
background genes of a large panel can be unobserved).  `oracle`/`naive` are
the st-recoverability Bayes-optimal / nearest-nucleus accuracies on interior
molecules (`meta.truth`); trivial datasets have exact truth by construction
(no displacement), so no oracle is computed.

Note: the st-recoverability "sparse" regime (2575 cells/mm²) falls into the
contract's `medium` density class (sparse is < 2500).

## Regenerating

```bash
export BAYSOR_BENCH_DATA=/path/to/bench-data      # default: <repo>/.bench-data
PY=.deps/bench/bin/python

$PY benchmarks/simulate/generate_all.py                    # all 21 datasets (~30 s)
$PY benchmarks/simulate/generate_all.py --list             # show the manifest
$PY benchmarks/simulate/generate_all.py --only sim_circles_gaps_g100
$PY benchmarks/simulate/generate_all.py --verify <id>      # byte-identity check
$PY benchmarks/simulate/generate_all.py --verify-all
```

Single datasets can also be built directly:

```bash
$PY benchmarks/simulate/trivial.py circles_gaps --id sim_circles_gaps_g100 \
    --seed 1201 --n-genes 100 -o "$BAYSOR_BENCH_DATA/sim/sim_circles_gaps_g100"
$PY benchmarks/simulate/strec.py --id strec_dense_s2_merfish \
    --packing 13625 --sigma 2.0 --model merfish --seed 910009 \
    -o "$BAYSOR_BENCH_DATA/sim/strec_dense_s2_merfish"
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
* **Prior column**: molecules within `r_nucleus` (3 µm, 2 µm for the small
  cells of `mixed_sizes`) of the nearest *true* centre carry that cell's id
  (1-based); 0 = no prior.  The st-rec datasets use the 3 µm nucleus prior of
  the upstream `methods_baysor.py`.
* **`meta.baysor`**: `scale_um` = the true mean cell radius (disc radius for
  round cells, `sqrt(area/cell / pi)` otherwise); `min_molecules_per_cell` 20
  and `prior_confidence` 0.5 follow the contract example / Xenium config.
  `config` is `null` for simulated datasets: all parameters are passed as
  explicit CLI flags by the runner.
* **st-recoverability datasets** are a 10-point design over packing ∈
  {2575, 6000, 13625} cells/mm², σ ∈ {1, 2, 3} µm and model ∈ {disjoint,
  MERFISH-realistic, Xenium-realistic} — each factor covered, the required
  dense + σ=2 + realistic cell present, not the full cartesian product.
  Oracle accuracy falls from 0.93 (sparse, σ=1) to 0.60 (dense, σ=3),
  reproducing the upstream frontier (their dense σ=2 demo value ≈ 0.72).

## External cache (`$BAYSOR_BENCH_DATA/cache/sim`, ~32 MB total)

| file | source | size |
|---|---|---|
| `st-recoverability/` | github clone at pinned commit (only `src/{generator,expression,oracle,config}.py` are imported) | 13 MB |
| `merfish_moffitt.h5ad` | squidpy `MERFISH_0.24.h5ad` mirror, figshare file 40038538 (same URL as upstream `realism.py`) | 3.9 MB |
| `xenium_breast_rep1_cell_feature_matrix.h5` | only this member of the 10x `Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip` (9.9 GB), fetched with `remotezip` range requests | 12 MB |
| `xenium_breast_rep1_cells.parquet` | likewise, from the same zip | 3.5 MB |

All four are re-downloadable by the scripts; nothing external is committed.

## Sanity check (BENCH-SIM part C)

`python benchmarks/simulate/sanity.py` runs the release binary on one trivial,
one sparse and one dense st-rec dataset with `meta.baysor` + `:prior`
(≤ 6 threads), writes run outputs to
`$BAYSOR_BENCH_DATA/runs/bench-sim-sanity/<id>/` and the report to
`sanity_check.json`.  Latest run:

| dataset | exit | wall time | majority-match accuracy | oracle | naive |
|---|---|---|---|---|---|
| sim_circles_gaps_g100 | 0 | 30.1 s | **0.987** | ≈ 0.99 (1 − background) | — |
| strec_sparse_s2_disjoint | 0 | 26.9 s | **0.628** | 0.870 | 0.832 |
| strec_dense_s2_merfish | 0 | 31.5 s | **0.619** | 0.723 | 0.635 |

All three runs finish and produce `segmentation*.csv`/loom/polygons output.
The trivial dataset is near-perfect as expected.  On the st-rec fields Baysor
sits below the oracle ceiling and around the naive nearest-nucleus baseline —
the same relationship the st-recoverability authors report for Baysor
(segmentation is merged: ~309/400 cells predicted in the sparse field).
Full metrics are BENCH-HARNESS's job; these numbers are a plausibility check
only.
