# Output files

This page defines the files written by `baysor run`. For bundle-level
selection, see [Outputs](outputs.md).

## Conventions

- `N` = number of retained molecules after input filtering and cropping
- `C` = number of final cells in the segmentation
- `G` = number of genes after input filtering
- cell names are written as `cell_<1-based-id>`
- noise / unassigned molecules are written with cell name `0`

## Legacy bundle

### `segmentation.csv`

Shape: `N` rows, one per retained molecule.

Exact column order:
`[transcript_id,] cell, gene, x, y [,z] [,confidence] [,cluster] [,ncv_color] [,assignment_confidence], is_noise`

Columns:

- `transcript_id`
  - present only for Xenium-origin inputs that preserve source transcript IDs
  - string
- `cell`
  - assigned cell name
  - string
  - `0` means noise / unassigned
- `gene`
  - gene name written from the current input panel
  - string
- `x`, `y`
  - molecule coordinates
  - float
- `z`
  - present for 3D inputs only
  - float
- `confidence`
  - molecule confidence estimated during noise modeling
  - float
- `cluster`
  - molecule-cluster label from the molecule clustering stage
  - integer
- `ncv_color`
  - per-molecule NCV color as hex RGB, for example `#4F9BD5`
  - string
- `assignment_confidence`
  - per-molecule posterior assignment confidence
  - float
- `is_noise`
  - last column, always present
  - for Xenium-origin inputs: literal `true` / `false`
  - otherwise: `1` / `0`

### `segmentation_cell_stats.csv`

Shape: `C` rows, one per final cell.

Exact column order:
`cell, x, y [,z] [,cluster] n_transcripts, density, elongation, area, avg_confidence [,avg_assignment_confidence] [,max_cluster_frac] [,lifespan]`

Columns:

- `cell`
  - cell name, for example `cell_17`
  - string
- `x`, `y`, optional `z`
  - centroid of molecules assigned to the cell
  - float
- `cluster`
  - assigned cell-level cluster label when molecule clustering was run
  - float-valued column containing integer labels
- `n_transcripts`
  - number of molecules assigned to the cell
  - float-valued column containing counts
- `density`
  - `n_transcripts / area` when area is available
  - float
- `elongation`
  - ratio of principal covariance eigenvalues of the cell molecule cloud
  - float
- `area`
  - 2D convex-hull area of the cell molecules
  - float
- `avg_confidence`
  - mean molecule confidence within the cell
  - float
- `avg_assignment_confidence`
  - mean posterior assignment confidence within the cell
  - float
- `max_cluster_frac`
  - fraction of molecules in the most frequent molecule-cluster label within
  the cell
  - float
- `lifespan`
  - number of traced iterations for the component GUID when tracing is
  available
  - float-valued column containing integer values

### `segmentation_polygons_2d.json`

Joined 2D cell polygons: one polygon per cell with at least one assigned
molecule (for 3D data, the polygons of all molecules of a cell pooled across
the z-stack). The cell set of this file is exactly the set of assigned cells
in `segmentation.csv`: a cell whose free-form boundary estimation fails or
produces fewer than three vertices gets a fallback rectangle around its
molecules instead of being omitted.

When `--polygon-format FeatureCollection` (default):

- root object type: `FeatureCollection`
- one feature per cell with a valid 2D polygon
- each feature has `id` (cell name), `geometry.type = "Polygon"`,
  `geometry.coordinates` (one closed outer ring in data coordinates), and
  `properties.cell` (cell name)

When `--polygon-format GeometryCollection`:

- root object type: `GeometryCollection`
- one geometry object per cell, each with `type: "Polygon"`, `coordinates`,
  and `cell` (cell name as a string)

When `--polygon-format GeometryCollectionLegacy`:

- root object type: `GeometryCollection`
- same layout as `GeometryCollection`, but `cell` is the integer part of the
  cell name (`cell_17` → `17`); this is the Baysor v0.7.1 handoff format for
  Xenium Ranger 3.x

When `--polygon-format none`: the file is omitted.

All format names are matched case-insensitively. Unknown values are
rejected.

### `segmentation_polygons_3d.json`

Written for 3D runs only.

Shape: JSON object keyed by layer name (for example `"z_003"` or a z value).

Structure: each value is a 2D polygon collection in the same schema as
`segmentation_polygons_2d.json`, describing the cell polygons within that
z-layer.

### `segmentation_counts.loom`

HDF5 layout (Loom spec `3.0.0`):

- `/matrix`
  - dataset type: `float32`
  - shape: `(G, C)` — rows are genes, columns are cells
- `/attrs/LOOM_SPEC_VERSION`
  - variable-length UTF-8 string
  - value: `"3.0.0"`
- `/row_attrs/Name`
  - UTF-8 string array of length `G`
  - gene names
- `/col_attrs/Name`
  - UTF-8 string array of length `C`
  - cell names
- `/col_attrs/CellID`
  - `float64`
  - length `C`
  - values `1, 2, ..., C`

Matrix values are counts.

### `segmentation_counts.tsv`

Dense tab-separated matrix:

- first column header: `gene`
- remaining column headers: one per cell name
- one data row per gene
- matrix shape on disk: `G x (1 + C)`

### `segmentation_params.dump.toml`

TOML document with the resolved run parameters after config-file loading and
CLI overrides, plus the CLI invocation as a comment on the first line.

### `segmentation_log.log`

Line-oriented text log: stage progress, iteration summaries, convergence
messages, save-step messages.

### `diagnostic_report.html`

Self-contained HTML document with diagnostic plots and summary panels.
Written only with `--plot`.

### `segmentation_plot.html`

Self-contained HTML document with an interactive molecule / cell
visualization. Written only with `--plot`.

## Parquet bundle

### `molecules.parquet`

Shape: `N` rows, one per retained molecule.

Exact column order:
`cell, gene, x, y [,z] [,confidence] [,cluster] [,ncv_color] [,assignment_confidence], is_noise`

Column types:

- `cell`, `gene`: UTF-8 string
- `x`, `y`, optional `z`: float64
- `confidence`, `assignment_confidence`: float64
- `cluster`: int32
- `ncv_color`: UTF-8 string
- `is_noise`: boolean

Note: `transcript_id` is not written to the parquet bundle; use `legacy` style
when you need Xenium Ranger compatibility.

### `cells.parquet`

Shape: `C` rows, one per final cell.

Columns: `cell` (UTF-8 string) followed by the same cell-stat columns and
order as `segmentation_cell_stats.csv` (all float64).

### `cell_boundaries.parquet`

Shape: one row per cell with a valid 2D polygon (2D runs; for 3D runs the
pooled polygons of the `2d` layer).

Columns:

- `cell`: UTF-8 string
- `n_vertices`: int32, number of polygon vertices excluding the duplicated
  closing vertex
- `geometry`: binary WKB polygon

GeoParquet metadata: primary geometry column `geometry`, encoding `WKB`,
geometry type `Polygon`, CRS unset.

### `cell_boundaries_3d.parquet`

Written for 3D runs only. Shape: one row per `(cell, layer)` polygon.

Columns: `cell` (UTF-8 string), `layer` (UTF-8 string), `n_vertices`
(int32), `geometry` (binary WKB polygon), with the same GeoParquet metadata
as `cell_boundaries.parquet`.

### `feature_matrix.h5`

10x-style HDF5 layout:

- `/matrix/barcodes` — UTF-8 string array of length `C`; cell names
- `/matrix/data` — int32, length `nnz`; nonzero values of the sparse matrix
- `/matrix/indices` — int32, length `nnz`; row indices into the feature axis
- `/matrix/indptr` — int64, length `C + 1`; column pointer array
- `/matrix/shape` — int64, length `2`; value `[G, C]`
- `/matrix/features/id` — UTF-8 string array of length `G`
- `/matrix/features/name` — UTF-8 string array of length `G`
- `/matrix/features/feature_type` — UTF-8 string array of length `G`; the
  current value for every gene is `Gene Expression`
- `/matrix/features/genome` — UTF-8 string array of length `G`; the current
  value for every gene is the empty string

The matrix is stored as CSC of the `G x C` count matrix; `id` and `name` are
both set to the gene names.

### `run_params.toml`

TOML document with the resolved run parameters for the parquet bundle.

### `run.log`

Line-oriented text log for the parquet bundle.

### `diagnostic_report.html` / `segmentation_plot.html`

Same as in the legacy bundle; written only with `--plot`.
