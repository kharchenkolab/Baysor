# Outputs

`baysor run` writes to the directory selected by `-o` (default `segmentation`).
The default output style is `legacy`; select Parquet output with:

```bash
baysor run -m 30 -s 8 --output-style parquet -o out molecules.csv
```

## Files

| Contents | Legacy (default) | Parquet |
| --- | --- | --- |
| Molecule assignments | `segmentation.csv` | `molecules.parquet` |
| Cell statistics | `segmentation_cell_stats.csv` | `cells.parquet` |
| Joined 2D polygons | `segmentation_polygons_2d.json` | `cell_boundaries.parquet` |
| Per-layer polygons (3D only) | `segmentation_polygons_3d.json` | `cell_boundaries_3d.parquet` |
| Gene × cell count matrix | `segmentation_counts.loom` or `segmentation_counts.tsv` | `feature_matrix.h5` |
| Config values and invocation | `segmentation_params.dump.toml` | `run_params.toml` |
| Log | `segmentation_log.log` | `run.log` |

For 3D data, joined polygons pool each cell's molecules across the z-stack.
In molecule tables, cell `0` means noise / unassigned; other cells are named
`cell_1`, `cell_2`, etc. See [Output file formats](output_files.md) for columns
and storage layouts.

Add `-p` / `--plot` for two self-contained browser reports in either style:

- `diagnostic_report.html` — quality checks and run diagnostics
- `segmentation_plot.html` — molecule and cell visualization

## Legacy output

Use legacy output for CSV / GeoJSON / Loom workflows and Xenium Ranger.
`--count-matrix-format tsv` changes the default Loom matrix to TSV.
`--polygon-format` selects `FeatureCollection` (default),
`GeometryCollection`, `GeometryCollectionLegacy` or `none` (no polygon files).

For [Xenium Ranger import](xenium.md#xenium-explorer-handoff), keep transcript
IDs and choose the polygon format for your Ranger version. Molecule
assignments and joined polygons use matching cell IDs; failed boundary
estimates get fallback rectangles.

## Parquet output

Use Parquet / GeoParquet tables and the 10x-style HDF5 matrix for downstream
Python, R or DuckDB analysis. Transcript IDs are not included in this bundle;
use legacy output for Xenium Ranger.

`--polygon-format` and `--count-matrix-format` are legacy-only and ignored for
Parquet output.
