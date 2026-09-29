# cellAdmix admixture audit as a Baysor benchmark metric

Runs the [cellAdmix-core](https://github.com/kharchenkolab/cellAdmix-core)
admixture audit (`fit.audit_admixture()`) on contract-format benchmark
molecules and writes a JSON result that scores a segmentation: **lower total
admixture = cleaner segmentation**. The audit estimates leaked molecules per
ordered cell-type pair from the spatial exposure gradient (target cells
stratified by the number of source-type neighbours among their 15 nearest
cells, source-marker excess over the zero-exposure baseline, extrapolated by
the markers' transcriptome share). It needs no reference annotation.

```
total_admixture_rate = Σ_pairs A_{S→T} / M          (over detected pairs)
```

## Files

| file | what |
|---|---|
| `install.sh`, `INSTALL.md` | build/install of the pinned cellAdmix Python bindings (+3 committed patches) |
| `fetch_pancreas.py` | download the 10x Xenium pancreas FFPE members (remotezip, ~155 MB) and cut contract-format crops |
| `audit.py` | **the CLI**: cellAdmix audit → `result.json` |
| `degrade.py` | deliberate degradation (border reassignment, dilation) → assignment parquet |
| `transfer.py` | transfer baseline cell types to another segmentation by molecule overlap |
| `align.py` | row-align a segmentation's per-molecule labels (e.g. Baysor output) → assignment parquet |
| `store.py` | `TabularCellAdmix` dataset wrapper + quick clustering |
| `summarize.py` | validation summary table + monotonicity checks |
| `run_validation.sh` | end-to-end validation (vendor vs degraded, seeds, runtime) |
| `run_baysor_validation.sh` | Baysor segmentation of the same crop → audit + comparison |
| `tests/` | pytest tests for the non-trivial logic |
| `results/` | committed baselines: `summary.{json,md}` + one audit JSON per variant |

## Quick start

```bash
./install.sh                       # build cellAdmix bindings (see INSTALL.md)
export BAYSOR_BENCH_DATA=/home/vpetukhov/Projects/Baysor/.bench-data
PY=.deps/bench/bin/python          # bench env python

# 1. dataset (cached after the first run)
python fetch_pancreas.py --crop-id pancreas_crop_quick --side-um 625 --max-molecules 150000

# 2. audit the vendor segmentation (typing = cellAdmix quick clustering)
python audit.py --molecules $BAYSOR_BENCH_DATA/cache/celladmix/datasets/pancreas_crop_quick/molecules.parquet \
  --cell-column cell_vendor --out vendor.json --threads 6 --seed 1

# 3. full validation + results/
./run_validation.sh && ./run_baysor_validation.sh

# 4. tests
python -m pytest tests/
```

## `audit.py`

```
--molecules molecules.parquet        contract table (x, y, gene, ...)
--assignment assignment.parquet      per-molecule cell int, 0 = unassigned,
                                     same row order; optional cell_label column
    XOR
--cell-column cell_vendor            segmentation column inside --molecules
--celltypes celltypes.parquet        fixed typing (else quick clustering)
--image membrane.tif                 recorded only; not used by the audit
--seed / --threads / --out
--cluster-resolution 1.0             quick-clustering resolution (= number of types)
--min-target-cells 200 --min-reference-cells 100 --n-pool 60 --neighbor-k 15
                                     audit parameters (defaults shown)
--save-celltypes out.parquet         write the typing actually used
```

Output JSON: `total_admixture_rate`, `total_admixture_molecules`,
`pairs_top` (top-20 detected pairs by estimated admixed molecules, with
per-pair `rate`, `coverage`, `q_value`, exposure/reference counts), counts
(molecules input/written/used, cells, genes, cell types, pairs
evaluated/detected), `filters`, all parameters (seed, threads, typing,
cluster resolution, fit parameters from the run manifest, audit parameters),
the pinned cellAdmix commit (+ verification against the source clone), input
sha256, and per-stage runtimes. Unassigned molecules are dropped exactly like
`keep_unassigned = FALSE` (the loader's token set — `""`, `"0"`, `"NA"`,
`"NaN"`, `"null"`, `"UNASSIGNED"`, `"unassigned"`, `"cell_0"` — is mirrored
in `audit.py`).

### `n_pool` default

cellAdmix's constructor default is `n_pool=20`; this benchmark defaults to
`n_pool=60`. On small panels (the crop has 376 genes) many genes crowd at
marker contrast ≈ 2, so the 20-gene pool cutoff reshuffles under small
perturbations and the pool's transcriptome coverage — the extrapolation
factor — swings between segmentations (observed: 0.500 → 0.018 for one pair,
which dropped it out of detection and inverted the metric's ordering).
With `n_pool=60` the dominant high-share markers stay in the pool and the
metric is stable; `--n-pool 20` is still available.

## Cell typing and comparability

The audit needs a cell type per cell, so **typing depends on the
segmentation** — this is the main comparability trap:

* **Re-clustering each segmentation** (quick clustering, fixed seed) types
  whatever cells that segmentation produced. Cell ids differ between
  segmentations, so nothing is shared: both the type partition *and* the
  cell universe change, and the metric moves for typing reasons as well as
  segmentation reasons. Measured on identical border-30%-degraded molecules:
  re-cluster typing gives 0.1171 vs 0.0946 for transferred typing — +24% on
  the *same* segmentation.
* **Type the baseline once and transfer** (recommended): cluster the vendor
  segmentation once (`--save-celltypes`), then give every other segmentation
  the same labels. When cell ids are shared (the degradation variants reuse
  vendor cell ids) the labels apply directly; when they differ (Baysor),
  `transfer.py` votes the baseline type by molecule overlap (majority over
  the baseline types of the molecules a target cell contains, ties by
  lexicographic order, cells without a single typed molecule dropped from the
  annotation — their molecules are excluded by `audit.py` and counted in
  `filters`).

**Recommendation: transfer.** It holds the typing constant so score
differences reflect segmentation differences only, which is the quantity a
segmentation benchmark must measure. Re-clustering each segmentation is
useful only as a diagnostic (this repo reports it as
`border30_recluster`), never as the comparison basis.

Two further typing details:

* `audit.py` clusters with `min_molecules=1`, `min_genes=1`,
  `cells_max=-1` so *every* cell gets a label — cells without a type would
  otherwise collapse into a spurious `"None"` pseudo-type inside the audit.
* The C++ clustering's parallel graph step is racy: on the same store and
  seed, 1 of 6 multi-threaded runs returned a different partition (total rate
  moved by 0.0009 ≈ 1.1%). `store.quick_cluster()` therefore runs
  single-threaded (clustering takes < 1 s on the crop), which made every
  repeated run bit-identical.

## Validation on the real pancreas crop

Dataset: `cache/celladmix/datasets/pancreas_crop_quick` — Xenium V1 FFPE
human pancreas, 450×450 µm densest-window crop, **141,351 molecules, 376
genes, 2,864 vendor cells** (qv ≥ 20, control probes/codewords removed,
rows sorted by (y, x)); 120,011 molecules carry a vendor cell label.

Degradations of the vendor segmentation (`degrade.py`):

* `border10/30` — the 10% / 30% of assigned molecules with the smallest
  distance to a foreign-labelled molecule are re-assigned to the cell owning
  that nearest neighbour (prefix ranking ⇒ 30% ⊇ 10% by construction);
  12,001 / 36,003 molecules moved, max cross-cell distance 0.91 / 1.52 µm.
* `dilate2` — every cell's convex hull is dilated by 2 µm: molecules within
  2 µm of another cell's hull flip to the nearest such cell (63,911 from
  neighbours), background molecules within 2 µm are absorbed (16,064).

Results (`results/summary.md`, all with the *same transferred typing* except
where noted):

| variant | total_admixture_rate | admixed molecules | detected/evaluated pairs | runtime |
|---|---|---|---|---|
| vendor | **0.075851** | 9,103 | 9/38 | 4.7 s |
| vendor (seed 2) | 0.075851 | 9,103 | 9/38 | 4.7 s |
| vendor (same seed, rerun) | 0.075851 | 9,103 | 9/38 | 4.7 s |
| border10 | 0.089800 | 10,777 | 10/39 | 4.1 s |
| border30 | 0.094625 | 11,356 | 10/39 | 4.1 s |
| dilate 2 µm | 0.143193 | 19,485 | 10/39 | 4.2 s |
| border30, re-clustered typing | 0.117148 | 14,059 | 8/37 | 4.6 s |
| Baysor (transferred typing) | 0.048090 | 6,748 | 5/27 | 3.9 s |

* **Monotonicity**: vendor (0.0759) < border10 (0.0898) < border30 (0.0946)
  and vendor < dilate2 (0.1432) — the audit reports monotonically higher
  admixture for worse segmentations (`summarize.py` checks PASS in
  `results/summary.md`).
* **Stochasticity**: three full runs (seed 1, seed 2, seed 1 again) are
  bit-identical — measured seed tolerance **0.000000**; the seed reaches only
  the NMF fit, which the audit does not use, and quick clustering is
  deterministic single-threaded. With multi-threaded clustering the
  observed run-to-run spread was ±0.0009; when typings *differ* (two
  different clustering draws) totals moved by up to ~0.006, so cross-typing
  comparisons need a tolerance of ≈ 0.005 absolute while
  same-typing comparisons are exact.
* **Runtime on the quick crop**: 3.9–4.7 s per audit (store build + typing +
  NMF fit + audit; 6 threads); the complete `run_validation.sh` (7 audits +
  3 degradations + summaries) takes ≈ 80 s.
* **Baysor vs vendor**: the Baysor run (35.7 s wall, `OMP_NUM_THREADS=6`,
  parameters from the crop's `meta.json`: scale 5 µm, scale-std 25%,
  min-molecules-per-cell 20, no prior) produced 1,783 cells (743 noise
  molecules unassigned). Its audit, typed by transferring the vendor
  clustering over molecule overlap (1,762 of 1,783 cells typed), scores
  **0.0481 — 37% lower admixture than the vendor segmentation (0.0759)**:
  on this crop Baysor's segmentation is markedly cleaner by the audit's own
  (reference-free) measure.

## Regenerating everything

```bash
./install.sh                                        # cellAdmix bindings
python fetch_pancreas.py                            # raw zip members + crop
./run_validation.sh                                 # degradations + audits + results/
./run_baysor_validation.sh                          # Baysor variant + final summary
python -m pytest tests/                             # 34 tests
```

All randomness is seeded (`--seed`, default 1); crop selection is
deterministic (densest window on a fixed grid, budget-driven shrink);
degradations are deterministic prefix/geometric rules.

## Data locations and sizes (this task's footprint)

```
$BAYSOR_BENCH_DATA/cache/celladmix/
  raw_pancreas/            148 MB   transcripts.parquet (8,073,840 rows) + cells.parquet,
                                    fetched with remotezip range requests
  datasets/pancreas_crop_quick/ 2.3 MB   molecules.parquet (141,351 rows, 376 genes)
                                    + meta.json
  work/validation/          90 MB   input stores, runs, assignment parquets, baysor_seg/
  src/                     1.1 GB   pinned cellAdmix-core clone + C++ test build tree
```

Total ≈ 1.35 GB (task budget: 150 GB). The committed `results/` files are
≈ 45 KB of JSON/Markdown; no data is committed.

## Known issues found along the way

1. **Upstream segfault**: parquet *tabular* store builds crash
   (use-after-free of the parquet reader; pure-C++ repro). Fixed by
   `patches/0002-keepalive-tabular-parquet-reader.patch` — worth reporting to
   kharchenkolab/cellAdmix-core.
2. **Racy parallel clustering**: see "Cell typing" above; worked around by
   single-threaded clustering.
3. **Pool-cutoff instability** on small panels: see `n_pool` above.
4. `tests/test_bridge.cpp` missing `#include <algorithm>` (gcc 15), patched
   locally as patch 0003 so the C++ tests build.
