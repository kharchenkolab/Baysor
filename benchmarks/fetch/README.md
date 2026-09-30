# BENCH-REALX: real Xenium datasets

Fetch + crop pipeline for the real Xenium part of the benchmark suite.
The dataset contract (files, `meta.json`, tiers) is defined in
[`../README.md`](../README.md); the dataset definitions live in
[`../datasets/real_xenium.yaml`](../datasets/real_xenium.yaml).

## Files

| file | purpose |
|---|---|
| `xenium.py` | argparse CLI: `verify`, `fetch`, `pick`, `build`, `record-hashes`, `report` |
| `xenium_common.py` | library: partial zip downloads, streaming transcript reads, crop selection, image-prior rasterisation, meta assembly |
| `download.py` | shared download helper (retries/backoff/Retry-After, Range resume, sha256, atomic rename); also usable by other fetchers |
| `sanity_run.sh` | run the Release Baysor binary on one dataset (command assembled from `meta.json`), record wall/RSS into `runs/<run_id>/<id>/timing.json` |
| `tests/` | pytest coverage of the non-trivial logic (`conftest.py` hosts the local HTTP server fixture) |

## Rebuild from scratch

```bash
export BAYSOR_BENCH_DATA=/path/to/.bench-data   # shared suite data root
PY=.deps/bench/bin/python                       # suite Python env
$PY benchmarks/fetch/xenium.py verify           # HEAD-check URLs + verify recorded hashes
$PY benchmarks/fetch/xenium.py fetch            # stream needed zip members to cache
$PY benchmarks/fetch/xenium.py pick --write     # only needed if a crop.bbox_um is null
$PY benchmarks/fetch/xenium.py build            # write real/<id>/ dataset dirs (+ output hashes)
$PY benchmarks/fetch/xenium.py record-hashes    # refresh output hashes without rebuilding
$PY benchmarks/fetch/xenium.py report           # inventory table
$PY -m pytest benchmarks/fetch/tests            # unit tests
```

`build` groups datasets by source bundle and reads each `transcripts.parquet`
once for all of its crops: row groups are iterated with pyarrow, only the
needed columns are materialised (`z_location` additionally when a crop has
`crop.keep_z`), and rows are filtered by `qv >= 20`, by the real-gene rule
(`is_gene` / `codeword_category` when present, otherwise the `NegControl*`,
`UnassignedCodeword*`, `DeprecatedCodeword*`, `Intergenic_Region*`,
`BLANK*` name prefixes) and by the crop box.  Nothing larger than a row
group plus the cropped result is ever held in memory, so the 2.2 GB Prime 5K
table is safe on a laptop-sized machine.  Builds are byte-reproducible:
`build <id> --out-root <tmp>` produces the same `molecules.parquet`/`meta.json`
sha256 as the committed manifest (`--out-root` disables hash recording).

Raw members are cached under `$BAYSOR_BENCH_DATA/cache/realx/` (~5.6 GB for
the current six bundles: five original sources plus breast Rep2) so
re-cropping never re-downloads.

## Downloads: retries, resume, hashes

All member downloads go through the shared [`download.py`](download.py):

* exponential-backoff retries for connection errors, timeouts and HTTP
  429/5xx, honouring `Retry-After` (capped at 5 minutes);
* HTTP `Range` resume of interrupted `.part` files (direct files;
  DEFLATE zip members restart from the member start because a deflate
  stream cannot be decompressed from an arbitrary offset - resuming the
  compressed range would need the full compressed prefix anyway);
* write to `<name>.part`, verify, then atomic rename;
* sha256 per member: streamed while downloading, recorded in the manifest
  under `downloaded_members` and re-verified whenever the cache is reused
  (a same-size corrupt file is detected and refetched); sizes are checked
  against the zip's advertised member size and the response
  `Content-Length`/`Content-Range`.

Other fetchers (e.g. the other-platforms `fetch/other.py`) can adopt the
helper with a few lines: `download.fetch_url(url, dst, expected_size=...,
expected_sha256=...)` does everything above for a direct URL.

## Dataset hashes

`build` (and `record-hashes` for already-built datasets) records the sha256
of each dataset's `molecules.parquet` and `meta.json` into
`real_xenium.yaml` under `outputs:`.  `verify` checks them (plus the cached
member hashes and the source URLs); `verify --skip-urls` is offline,
`verify --require-built` also fails on datasets that are missing on disk.
This is what pins the dataset contents the baselines (kept locally in
`$BAYSOR_BENCH_DATA/baselines/`, never committed) were computed on.

## Crop selection criteria

Existing crops keep the bboxes they were picked with; the criteria below
are documented here and applied to the **new** crops only (the `criteria:`
block in `real_xenium.yaml`, measured values recorded under
`crop.pick`).  A candidate box must satisfy all of them in addition to the
molecule/cell/coverage limits; criteria are never relaxed (only the
`density_hint` falls back to "densest quartile" as before):

* **Cell-type composition diversity** - `criteria.min_clusters`: the number
  of distinct vendor clusters (graphclust `analysis/.../clusters.csv` from
  the bundle's analysis outputs, when the bundle ships one; breast Rep1,
  Rep2 and mouse brain do, lung/ovarian/pancreas do not) whose cell
  centroids fall inside the box.  Requires a `clusters` entry in the
  dataset's `members`.
* **Tissue edge or empty space** - `criteria.max_coverage`: the fraction of
  occupied histogram bins in the box must not exceed this (e.g. 0.8), and
  `criteria.min_empty_border`: at least this fraction of the box's
  2-bin-wide border ring must be empty, i.e. the empty space touches the
  crop edge and is a real tissue boundary, not an interior hole.
* **Density** - the existing `density_hint` (`any` / `dense` >= 7000
  cells/mm² / `stromal` 1000-7000 cells/mm²).

## Datasets (18)

The original ten crops (`xenium_pancreas_377_{full,quick}`,
`xenium_breast_rep1_{dense_full,dense_quick,stroma_quick}`,
`xenium_mouse_brain_ff_quick`, `xenium_prime5k_ovarian_{full,quick}`,
`xenium_lung_cancer_quick`) are unchanged.  New ids:

| id | tier | what it adds |
|---|---|---|
| `xenium_breast_rep1_z_quick` | quick | keeps `z` (25-45 µm sub-section depth); `meta.baysor.extra_args` carries `-z z`, Baysor runs `BmmData (3D)` |
| `xenium_breast_rep1_imageprior_quick` | quick | image prior: vendor nucleus boundaries rasterised to `images/nuclei_labels.tif`, `prior: "image:images/nuclei_labels.tif"` |
| `xenium_mouse_brain_ff_edge_quick` | quick | tissue edge / empty space (coverage 0.66, 40% of the border ring empty) |
| `xenium_breast_rep2_dense_quick` | quick | second section of the same tumour (Rep2, same panel), dense + composition criteria |
| `xenium_breast_rep1_dense_admix` | full | admixture gate: 3559 vendor cells / 550k molecules |
| `xenium_breast_rep1_stroma_admix` | full | admixture gate: 2210 vendor cells / 550k molecules |
| `xenium_mouse_brain_ff_admix` | full | admixture gate: 2645 vendor cells / 550k molecules |
| `xenium_prime5k_ovarian_admix` | full | admixture gate: 4482 vendor cells / 550k molecules (large panel, larger box) |
| `xenium_lung_cancer_admix` | full | admixture gate: 7560 vendor cells / 550k molecules |

Every `*_admix` crop carries >= 2000 vendor cells at <= ~600k molecules,
which is what the cellAdmix audit needs for statistical power.

## Image prior frame

Baysor maps a molecule at `(x, y)` µm to pixel `(round(x) - 1, round(y) - 1)`
(`load_prior_from_image` in `src/data_loading/prior_segmentation.cpp`), so the
label raster must be a **1 µm/pixel grid aligned with the absolute coordinate
origin**: pixel `(c, r)` covers `x in [c+0.5, c+1.5)` and the image spans
`[0, ceil(x1)) x [0, ceil(y1))` pixels.  `rasterize_nucleus_labels` tests
polygon/pixel-box intersections exactly with shapely, so every molecule inside
a nucleus lands on that nucleus's pixel (~99.7% measured on the Rep1 crop;
the rest are vendor `overlaps_nucleus=0` edge cases).  Labels are dense
`1..K` per crop (uint16, 0 = background), deterministic across rebuilds.

The label TIFF is written **uncompressed on purpose**: Baysor reads it row by
row starting at the first row of the molecule window (usually mid-strip), and
the bundled libtiff refuses sub-strip scanline reads of Deflate data
(`Compression algorithm does not support random access`).  Uncompressed strips
seek arbitrarily; the 1325x2450 Rep1 frame is 6.5 MB.

## Notes and gotchas

* The four `morphology_focus/morphology_focus_NNNN.ome.tif` files of a bundle
  form **one multi-file OME image**: channel `N` lives in file `NNNN`
  (0 = DAPI, 1 = boundary stain), and each file's OME XML references its
  siblings, so tifffile resolves every channel from any member.  Windowed
  reads go through `series.aszarr()` and touch only the needed tiles.
* `configs/xenium.toml` maps vendor column names (`x_location`, ...) and sets
  `unassigned_label = "UNASSIGNED"`.  Our `molecules.parquet` uses contract
  names and a numeric `prior` column (0 = no prior), so `meta.baysor.extra_args`
  carries explicit `-x x -y y -g gene --qv-column qv --unassigned-prior-label 0`.
  Without the last flag Baysor would treat label `0` as a real prior segment.
* Breast Rep1 transcripts contain eight `antisense_*` features that are not in
  the panel gene list.  They are kept (the fetch filter drops exactly the
  prefixes listed above) and are excluded at run time by
  `configs/xenium.toml`'s `exclude_genes`.
* `scale_um` = 1.5 × √(median vendor `nucleus_area` / π) over the cells whose
  centroid is in the crop; the exact method string is stored in
  `meta.baysor.scale_um_method`.
* Sanity runs: `benchmarks/fetch/sanity_run.sh <dataset_id>` (6 threads,
  `/usr/bin/time -v`, outputs in `$BAYSOR_BENCH_DATA/runs/sanity_realx/<id>/`).
  The command line is assembled from `meta.json` exactly like the harness
  does it (config, `extra_args`, scale, prior column/image/none, prior
  confidence), so 3D and image-prior datasets run through the same script.
