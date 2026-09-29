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
  environment.yml        Python environment for the suite
  datasets/              manifests: one YAML per dataset group (sim, real_xenium, real_other)
  simulate/              generators for simulated datasets
  fetch/                 download + crop scripts for real datasets
  harness/               runner, metrics, baseline comparison, reports
  celladmix/             cellAdmix admixture audit on a segmentation
  baselines/             committed baseline metrics (small JSON/CSV per dataset)
```

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
