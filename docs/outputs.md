# Outputs

`baysor run` writes one output directory (`-o/--output`, default
`segmentation`) in one of two styles (`--output-style`):

- `legacy` (default) — CSV / GeoJSON / Loom, compatible with the classic
  Baysor outputs and `xeniumranger import-segmentation`
- `parquet` — Parquet / GeoParquet / 10x-style HDF5 for downstream analysis
  in Python / R / DuckDB

```bash
baysor run --output-style legacy -o out ...
baysor run --output-style parquet -o out ...
```

For file-by-file definitions, including exact column order and storage layout,
see [Output files](output_files.md).

## Legacy output

- `segmentation.csv` — per-molecule segmentation table
- `segmentation_cell_stats.csv` — per-cell statistics
- `segmentation_polygons_2d.json` — joined cell polygons, GeoJSON (for 3D
  runs: all molecules pooled across the z-stack)
- `segmentation_polygons_3d.json` — per-layer cell polygons, GeoJSON (3D runs
  only)
- `segmentation_counts.loom` or `segmentation_counts.tsv` — count matrix
  (`--count-matrix-format`)
- `segmentation_params.dump.toml` — resolved run parameters
- `segmentation_log.log` — run log

With `--plot`, two extra HTML files are written: `diagnostic_report.html` and
`segmentation_plot.html`.

### Xenium compatibility

For Xenium-origin inputs (started from `experiment.xenium`), the `legacy`
bundle automatically adds the fields needed by
`xeniumranger import-segmentation`:

- `segmentation.csv` includes `transcript_id` and writes `is_noise` as
  `true` / `false`
- `segmentation_polygons_2d.json` uses a GeoJSON `FeatureCollection` with
  `properties.cell`; `--polygon-format GeometryCollectionLegacy` switches to
  the integer-id `GeometryCollection` layout required by Xenium Ranger 3.1
  and earlier

The two files contain exactly the same cells, so the import does not fail
with `EmptyCellsError` or `MissingCellPolygon`.

Use `legacy` when the result will be handed off to Xenium Ranger / Xenium
Explorer.

## Parquet output

- `molecules.parquet` — per-molecule segmentation table
- `cells.parquet` — per-cell statistics
- `cell_boundaries.parquet` — joined 2D cell polygons, GeoParquet
- `cell_boundaries_3d.parquet` — per-layer cell polygons, GeoParquet (3D runs
  only)
- `feature_matrix.h5` — 10x-style HDF5 feature-barcode matrix
- `run_params.toml` — resolved run parameters
- `run.log` — run log

With `--plot`, the same two HTML files as in `legacy` style are written.

## Legacy-only flags

`--polygon-format` and `--count-matrix-format` only affect the `legacy`
bundle; in `parquet` style they are ignored with a warning.

## Choosing between styles

Use `legacy` when you want the classic Baysor outputs, GeoJSON polygons, or
Xenium Ranger compatibility. Use `parquet` when you want Parquet / GeoParquet
tables and a 10x-style HDF5 matrix for Python / R / DuckDB analysis.
