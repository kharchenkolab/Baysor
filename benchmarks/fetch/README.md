# BENCH-REALX: real Xenium datasets

Fetch + crop pipeline for the real Xenium part of the benchmark suite.
The dataset contract (files, `meta.json`, tiers) is defined in
[`../README.md`](../README.md); the dataset definitions live in
[`../datasets/real_xenium.yaml`](../datasets/real_xenium.yaml).

## Files

| file | purpose |
|---|---|
| `xenium.py` | argparse CLI: `verify`, `fetch`, `pick`, `build`, `report` |
| `xenium_common.py` | library: partial zip downloads, streaming transcript reads, crop selection, meta assembly |
| `sanity_run.sh` | run the Release Baysor binary on one dataset, record wall/RSS into `runs/<run_id>/<id>/timing.json` |
| `tests/test_xenium_common.py` | pytest coverage of the non-trivial logic |

## Rebuild from scratch

```bash
export BAYSOR_BENCH_DATA=/path/to/.bench-data   # shared suite data root
PY=.deps/bench/bin/python                       # suite Python env
$PY benchmarks/fetch/xenium.py verify           # HEAD-check source URLs/sizes
$PY benchmarks/fetch/xenium.py fetch            # stream needed zip members to cache
$PY benchmarks/fetch/xenium.py pick --write     # only needed if a crop.bbox_um is null
$PY benchmarks/fetch/xenium.py build            # write real/<id>/ dataset dirs
$PY benchmarks/fetch/xenium.py report           # inventory table
$PY -m pytest benchmarks/fetch/tests            # unit tests
```

`build` groups datasets by source bundle and reads each `transcripts.parquet`
once for all of its crops: row groups are iterated with pyarrow, only the
needed columns are materialised, and rows are filtered by `qv >= 20`, by the
real-gene rule (`is_gene` / `codeword_category` when present, otherwise the
`NegControl*`, `UnassignedCodeword*`, `DeprecatedCodeword*`,
`Intergenic_Region*`, `BLANK*` name prefixes) and by the crop box.  Nothing
larger than a row group plus the cropped result is ever held in memory, so the
2.2 GB Prime 5K table is safe on a laptop-sized machine.

Raw members are cached under `$BAYSOR_BENCH_DATA/cache/realx/` (~5 GB for the
current five bundles) so re-cropping never re-downloads.

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
  prefixes listed in the task) and are excluded at run time by
  `configs/xenium.toml`'s `exclude_genes`.
* `scale_um` = 1.5 × √(median vendor `nucleus_area` / π) over the cells whose
  centroid is in the crop; the exact method string is stored in
  `meta.baysor.scale_um_method`.
* Sanity runs: `benchmarks/fetch/sanity_run.sh <dataset_id>` (6 threads,
  `/usr/bin/time -v`, outputs in `$BAYSOR_BENCH_DATA/runs/sanity_realx/<id>/`).
