# Baysor benchmark suite

Regression and quality benchmarks for the Baysor segmentation algorithm on
cropped real datasets and simulated datasets with known ground truth.

The suite answers two questions:

1. **No algorithm change** (refactoring, performance work, bug fixes that must
   not alter results): simulated-data metrics stay within the noise floor of the
   stored baseline, and real-data segmentations stay essentially the same as
   the stored baseline segmentation.
2. **Algorithm change**: accuracy on simulated data goes up, and the cellAdmix
   admixture audit on real data is not worse than the baseline.

Data is never committed. Code, dataset manifests and small baseline metric
files are.

## Layout

```
benchmarks/
  README.md              this file: layout and the dataset contract
  DATASETS.md            generated inventory of every dataset + coverage matrix
  environment.yml        Python environment for the suite
  datasets/              manifests: one YAML per dataset group (sim, real_xenium,
                         real_other) + suites.yaml (the regular/release suites)
  simulate/              generators for simulated datasets
  fetch/                 download + crop scripts for real datasets
  harness/               runner, metrics, baseline comparison, reports
  celladmix/             cellAdmix admixture audit on a segmentation
  baselines/             committed baseline metrics (small JSON/CSV per dataset)
```

## Dataset inventory

[`DATASETS.md`](DATASETS.md) lists every dataset (id, kind, tier,
platform/generator, tissue/scenario, genes, molecules, area, cells/mm²,
density and gene-panel class, 2D/3D, prior, images, `admixture_capable`,
source, notes) grouped by kind and platform, plus the density × gene-panel
coverage matrix per kind. It also carries the measured **resource columns**
(6-thread CPU time ± SD, 6-thread wall and peak RAM, 1-thread wall/RAM,
cellAdmix audit time — from the committed
[`baselines/bugfixes-35e8a7e/resources.csv`](baselines/bugfixes-35e8a7e/resources.csv),
regenerated with [`harness/resources.py`](harness/resources.py) from the
existing runs; `TODO` = never measured, never guessed) and a **Suite**
column (membership in [`datasets/suites.yaml`](datasets/suites.yaml)).
Regenerate it with:

```bash
$PY benchmarks/harness/inventory.py            # writes benchmarks/DATASETS.md
$PY benchmarks/harness/inventory.py --check    # exit 1 when stale
```

Validate the datasets against this contract (columns/dtypes/sort, meta
fields and enums, stats consistency, image/prior/config references,
manifest sha256) with:

```bash
$PY benchmarks/harness/validate_datasets.py    # exit 1 on metadata errors
$PY benchmarks/harness/validate_datasets.py --strict   # data findings too
```

## Workflow

See [`harness/README.md`](harness/README.md) for the benchmark workflow:
running datasets (`run.py`), metric definitions, baselines and the measured
noise floor (including Baysor's determinism findings), and the
`--expect same` / `--expect improved` comparison (`compare.py`, one-shot
`bench.sh`).

## How to test a change

Official baselines of the current algorithm:

| baseline | flavor | contents |
|---|---|---|
| [`baselines/bugfixes-35e8a7e-t1`](baselines/bugfixes-35e8a7e-t1/) | `identical` | quick tier, 1 thread, 1 replicate, no cellAdmix |
| [`baselines/bugfixes-35e8a7e`](baselines/bugfixes-35e8a7e/) | noise floor | quick + full tier, 6 threads, 3 replicates (full tier: see its README), cellAdmix audit with stable typing |

Setup used by every command below:

```bash
export BAYSOR_BENCH_DATA=/home/vpetukhov/Projects/Baysor/.bench-data
PY=.deps/bench/bin/python
B=/path/to/baysor                 # your build of the same sources
```

### The suites (`datasets/suites.yaml`)

Two committed suites (schema in [`harness/suites.py`](harness/suites.py),
resolution via `run.py --suite` / `compare.py --suite`). The times are
estimates from
[`baselines/bugfixes-35e8a7e/resources.csv`](baselines/bugfixes-35e8a7e/resources.csv)
(measured Baysor wall/CPU × replicates + cellAdmix audit; reproduce with
`$PY benchmarks/harness/suites.py --suite <name>`):

| suite | steps (compare mode) | coverage | est. wall | est. CPU |
|---|---|---|---|---|
| `regular` | `exact`: 1 thr × 1 rep → `identical` vs `-t1`; `noise`: 6 thr × 1 rep + audit → `same` vs `bugfixes-35e8a7e` | 4-dataset bitwise subset; 23-dataset coverage list | 3.8 + 14.4 = **18.1 min** core (+ ~1–2 min metrics/typing ≈ **~20 min**) | ~45 CPU-min |
| `release` | `quick6`+`full6`: 6 thr × 3 rep + audit → `same` (one run folder); `quick1`: 1 thr × 1 rep → `identical` | every dataset (78 = 65 quick + 13 full); the quick tier again at 1 thread | 156.7 + 232.4 + 87.9 = **477 min ≈ 8 h** (+ metrics bookkeeping; the historical `benchbase-b` quick+full passes observed ≈ 7 h against the 389 min core) | ~24 CPU-h |

* `regular` runs with `--no-ami`: AMI is informational (no gate reads it)
  but costs ~30–60 s of metrics time per sim replicate; `release` computes
  AMI so regenerated baselines keep their current contents.
* `regular` coverage: every `trivial.py` scenario **with and without
  prior** (12 datasets); st-recoverability sparse and dense; one 3D
  simulation; every real platform at quick tier (Xenium, CosMx, MERFISH,
  ISS, osmFISH, STARmap); gene-panel classes `tiny`–`huge` through cheap
  quick crops (`strec_*` tiny, `sim_circles_gaps_g100`/ISS small, pancreas
  medium, CosMx/STARmap large, `xenium_prime5k_ovarian_quick` huge).
* **The `*_admix` crop: yes, one fits** — `xenium_lung_cancer_admix`
  (~97 s/rep Baysor + ~6 s audit ≈ 103 s of the ~20 min budget) joins the
  noise step, so the cellAdmix audit and the admixture gate
  (`compare --expect improved`, `--celltypes-from bugfixes-35e8a7e` fixed
  typing/pairs) are exercised on a full-size admixture crop in every
  regular run; the quick Xenium crops in the list exercise the audit on
  quick data too. Without it the gate would still evaluate on
  `admixture_capable` quick crops, but never on the crop class the audit
  was calibrated for.
* Steps sharing a `group` run in the same folder `runs/<id>`; the group
  holding the suite's `identical` step keeps the bare `--run-id`
  (1-thread output-path-length sensitivity: keep it ≤ 17 characters),
  other groups get `<id>-<group>` (`<id>-noise`).

Resolve both suites without running Baysor (validates dataset ids,
baselines, run-ids and prints the estimates):

```bash
$PY benchmarks/harness/run.py --suite regular --run-id dry --dry-run
$PY benchmarks/harness/run.py --suite release --run-id dry --dry-run
$PY benchmarks/harness/suites.py --suite regular     # plan + estimates only
```

### Step 0: the C++ unit tests (part of `regular`)

```bash
cmake --preset tests && cmake --build --preset tests
ctest --test-dir build/tests --output-on-failure
# ~8 s (Release) / ~90 s (coverage build)
```

### Every change/PR: the `regular` suite (~20 min)

```bash
benchmarks/harness/bench.sh --baysor $B --preset regular --run-id reg-1
```

Runs both steps and compares every group (`exact` → `identical`, `noise`
→ `same`); exit 0 = pass, 1 = regression, 2 = setup error. `same` keeps
its strict *unchanged-binary* rule per check, but the **suite verdict
downgrades a `same` group whose only failing checks are the
`binary_sha256` provenance rows** (the normal result of a rebuilt binary
whose metrics all stayed within the noise floor) and prints a note — the
bitwise `identical` gate is then the verdict on behaviour.

Variants (override the non-bitwise group's mode):

* **algorithm change** (assignments are *supposed* to change):
  `bench.sh --preset regular --expect improved` — judges only the noise
  group (mean gain > noise + no per-dataset regression + admixture gate
  on `xenium_lung_cancer_admix`); the bitwise group is skipped by design.
* **harness/data/config change, binary untouched**:
  `bench.sh --preset regular --expect same` — strict mode, the sha gate
  fails if the binary really changed.
* **refactoring/bug fix that must not change results**: the default
  invocation; `identical` must pass and the noise metrics must stay
  within tolerance.

The legacy single-step presets `refactor` (1 thr × 1 rep, `identical`) and
`algorithm` (6 thr × 3 rep, `improved`) still exist for ad-hoc runs — see
[`harness/README.md`](harness/README.md).

### Before a release: the `release` suite (~8 h)

```bash
# algorithm release: full noise-floor run judged on improvement
benchmarks/harness/bench.sh --baysor $B --preset release --run-id rel-1 \
    --expect improved
# unchanged-algorithm release: default (--expect same) + bitwise gate
benchmarks/harness/bench.sh --baysor $B --preset release --run-id rel-2
```

Everything: all quick and full datasets at 6 threads × 3 replicates with
the cellAdmix audit (on the `*_admix` crops as on every real dataset),
plus 1-thread `identical` over the whole quick tier. After it passes,
freeze the new baselines (below).

Validate suite resolution without running Baysor (also shown above):
`run.py --suite <name> --dry-run` prints every step's datasets, threads,
replicates, timeouts and estimated time.

### Updating the baselines after an accepted change

```bash
# rerun both configurations with the new binary
$PY benchmarks/harness/run.py --baysor $B --datasets quick --run-id benchbase-b \
    --replicates 3 --threads 6 --timeout 1800 --skip-existing \
    --celltypes-from bugfixes-35e8a7e --label <new-sha>
$PY benchmarks/harness/run.py --baysor $B --datasets quick --run-id benchbase-t1 \
    --replicates 1 --threads 1 --timeout 1800 --no-celladmix \
    --skip-existing --label <new-sha>
# freeze them (adds --allow-incomplete/--identical as needed; see baseline.py)
$PY benchmarks/harness/baseline.py create --run-id benchbase-b \
    --name bugfixes-<new-sha> --force
$PY benchmarks/harness/baseline.py create --run-id benchbase-t1 \
    --name bugfixes-<new-sha>-t1 --force --identical
$PY benchmarks/harness/baseline_summary.py --baseline bugfixes-<new-sha>
```

(`--skip-existing` only reuses replicates whose binary sha256, thread
count and scale factor match, so a new binary reruns everything; the full
tier runs are appended by repeating the command with `--datasets full`
`--timeout 5400`.) See [`baselines/bugfixes-35e8a7e/README.md`](baselines/bugfixes-35e8a7e/README.md)
for the exact official-baseline invocations of this binary.

## Data location

All data lives under `$BAYSOR_BENCH_DATA` (default: `<repo>/.bench-data`,
gitignored), in one directory per dataset:

```
$BAYSOR_BENCH_DATA/
  sim/<dataset_id>/...
  real/<dataset_id>/...
  runs/<run_id>/<dataset_id>/...     harness outputs (Baysor results, metrics)
  cache/                             raw downloads kept for re-cropping (may be deleted)
```

## Dataset contract

Every dataset directory contains:

### `molecules.parquet` (required)

| column | type | required | meaning |
|---|---|---|---|
| `x`, `y` | float64 | yes | coordinates in µm |
| `z` | float64 | no | µm; present only for 3D datasets |
| `gene` | string (dictionary ok) | yes | gene name; control probes and blank codewords already removed |
| `qv` | float32 | no | vendor quality value, if available (already filtered by the fetch script) |
| `prior` | int32 | no | prior segmentation label for Baysor's `:prior` option (0 = no prior), e.g. the vendor nucleus id |
| `cell_vendor` | string | real only, if available | the vendor's cell assignment (empty = unassigned), kept for reference comparisons |
| `cell` | int32 | sim only | ground-truth cell id (0 = background / noise molecule) |
| `interior` | bool | sim only, optional | molecule counts toward accuracy metrics; excludes edge effects. Defaults to all true |
| `celltype` | string | sim only, optional | ground-truth cell type of the true cell |

Rows are sorted by (`y`, `x`). No other columns are required. Generators may
add extra columns prefixed with `aux_`.

### `meta.json` (required)

```json
{
  "id": "xenium_pancreas_377_crop1",
  "kind": "real",
  "tier": "quick",
  "platform": "Xenium",
  "source": {"url": "...", "doi": "...", "license": "...", "original_dataset": "...", "retrieved": "2026-09-29"},
  "crop": {"bbox_um": [x0, y0, x1, y1], "z_range_um": null, "note": "why this region"},
  "stats": {"n_molecules": 0, "n_genes": 0, "area_um2": 0.0, "molecules_per_um2": 0.0,
            "n_vendor_cells": 0, "vendor_cells_per_mm2": 0.0},
  "difficulty": {"cell_density": "sparse|medium|dense", "gene_panel": "tiny|small|medium|large|huge", "notes": ""},
  "baysor": {"scale_um": 5.0, "scale_std": "25%", "min_molecules_per_cell": 20,
             "prior": "none" ,
             "prior_confidence": 0.5,
             "config": "configs/xenium.toml",
             "extra_args": []},
  "images": [{"name": "dapi", "file": "images/dapi.tif", "pixel_size_um": 0.2125, "origin_um": [x0, y0]}],
  "truth": null
}
```

- `kind` is either `real` or `sim`.
- `tier` is either `quick` (at most about 150k molecules, runs in about a
  minute) or `full` (at most about 3M molecules).
- `baysor.prior` is one of:
  - `"none"`
  - `"column"`: use the `prior` column as `:prior`
  - `"image:<relative path>"`: a label TIFF in the dataset directory, same
    pixel frame as `images`
- `images` is optional. It holds cropped DAPI and membrane/boundary stains, when
  the source has them. They are used by cellAdmix membrane scoring and by future
  image-aware methods.
- For `sim`, `truth` holds the generator name and version or commit, all
  parameters, the seed and, when available, the oracle (best achievable)
  assignment accuracy. For `real` it stays `null`.
- `gene_panel` classes:

  | class | genes |
  |---|---|
  | `tiny` | < 50 |
  | `small` | 50–250 |
  | `medium` | 250–700 |
  | `large` | 700–2000 |
  | `huge` | > 2000 |

- `cell_density` classes, by nucleus or cell density:

  | class | cells/mm² |
  |---|---|
  | `sparse` | < 2500 |
  | `medium` | 2500–7000 |
  | `dense` | > 7000 |

### Optional files

- `images/*.tif`: cropped stains, single-channel, uint16 or uint8.
- `reference/`: vendor cell or nucleus boundaries for the crop (parquet), or a
  reference annotation.
- `README.md`: provenance notes specific to the dataset.

## Environment

The Python environment used by the whole suite is created by:

```bash
micromamba create -p .deps/bench -f benchmarks/environment.yml
```

On the development machine it is already created at
`/home/vpetukhov/Projects/Baysor/.deps/bench`.
