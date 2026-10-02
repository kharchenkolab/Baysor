# Output file formats

Reference for the files written by `baysor run`. For choosing a bundle, see
[Outputs](outputs.md).

## Conventions

- `N` — retained molecules after filtering and cropping
- `C` — final cells; `G` — retained genes
- Cell names are `cell_<1-based-id>`; noise / unassigned molecules have cell
  name `0`.
- Bracketed columns below are conditional: z for 3D, cluster fields when
  clustering runs, NCV colors unless skipped, and confidence / history fields
  when available.

## Legacy bundle

### `segmentation.csv`

One row per retained molecule (`N` rows). Column order:

```text
[transcript_id,] cell, gene, x, y [,z] [,confidence] [,cluster] [,ncv_color] [,assignment_confidence], is_noise
```

| Column | Meaning |
| --- | --- |
| `transcript_id` | Numeric source transcript ID, present when retained from the input. |
| `cell` | Assigned cell name, or `0` for noise. |
| `gene` | Gene name. |
| `x`, `y`, `z` | Molecule coordinates. |
| `confidence` | Molecule confidence from the noise model. |
| `cluster` | Integer molecule-cluster label. |
| `ncv_color` | NCV color as hex RGB, e.g. `#4F9BD5`. |
| `assignment_confidence` | Assignment stability estimated from recent iteration history. |
| `is_noise` | Always last; `true` / `false` when source transcript IDs are retained, `1` / `0` otherwise. |

### `segmentation_cell_stats.csv`

One row per final cell (`C` rows). Column order:

```text
cell, x, y [,z] [,cluster], n_transcripts, density, elongation, area, avg_confidence [,avg_assignment_confidence] [,max_cluster_frac] [,lifespan]
```

| Column | Meaning |
| --- | --- |
| `cell` | Cell name. |
| `x`, `y`, `z` | Centroid of assigned molecules. |
| `cluster` | Cell-level cluster label. |
| `n_transcripts` | Number of assigned molecules. |
| `density` | Molecule count / area when area is available. |
| `elongation` | Ratio of the two x/y covariance eigenvalues. |
| `area` | 2D convex-hull area of the cell's molecules. |
| `avg_confidence` | Mean molecule confidence. |
| `avg_assignment_confidence` | Mean assignment stability. |
| `max_cluster_frac` | Fraction of molecules in the most frequent molecule cluster. |
| `lifespan` | Consecutive stored iterations in which the component exists; limited by the retained history window. |

All columns except `cell` are numeric; labels and counts are stored as
floating-point values. Low `avg_assignment_confidence` can flag unstable
assignments, while low `max_cluster_frac` can flag mixtures of molecule
clusters. These are diagnostics, not automatic cell-quality filters.

### `segmentation_polygons_2d.json`

One joined polygon per assigned cell. For 3D data, molecules are pooled
across z. A failed boundary estimate gets a fallback rectangle so that the
cell set matches `segmentation.csv`.

`--polygon-format` controls the GeoJSON schema:

| Value | Schema |
| --- | --- |
| `FeatureCollection` (default) | One feature per cell, with `id` and `properties.cell` set to the cell name; `geometry.type` is `Polygon`. |
| `GeometryCollection` | One polygon geometry per cell with a string `cell` field. |
| `GeometryCollectionLegacy` | Same collection, but `cell` is an integer (`cell_17` → `17`) for Xenium Ranger 3.x. |
| `none` | No polygon files written. |

Coordinates use data units and one closed outer ring. Format names are
case-insensitive; unknown names are rejected.

### `segmentation_polygons_3d.json`

Written for 3D runs only: a JSON object keyed by layer (e.g. `z_003`), with a
2D polygon collection per layer using the schema above.
[`[plotting] max_z_slices`](configuration.md#plotting) (default `10`) limits
the layer count; stacks with more distinct z values are binned first.

### `segmentation_counts.loom`

Counts in Loom 3.0.0 format:

| HDF5 path | Type and shape | Contents |
| --- | --- | --- |
| `/matrix` | float32, `(G, C)` | Gene × cell counts. |
| `/attrs/LOOM_SPEC_VERSION` | UTF-8 string array, length 1 | `3.0.0`. |
| `/row_attrs/Name` | UTF-8 string array, length `G` | Gene names. |
| `/col_attrs/Name` | UTF-8 string array, length `C` | Cell names. |
| `/col_attrs/CellID` | float64 array, length `C` | IDs `1` … `C`. |

### `segmentation_counts.tsv`

Dense tab-separated counts: the first column is `gene`, followed by one
column per cell name, with one row per gene. Shape: `G × (1 + C)`.

### `segmentation_params.dump.toml`

Config values used by the pipeline, with the CLI invocation in a comment on
the first line. Keep the original config and command for reproducibility:
not every setting is serialized (for example, `tol` is omitted).

### `segmentation_log.log`

Text log with stage progress, iteration summaries, convergence and save
messages.

### `diagnostic_report.html`

Self-contained HTML diagnostics, written only with `--plot`.

### `segmentation_plot.html`

Self-contained molecule / cell visualization, written only with `--plot`.

## Parquet bundle

### `molecules.parquet`

`N` rows, with the same column order as `segmentation.csv` except that
`transcript_id` is not written:

| Columns | Type |
| --- | --- |
| `cell`, `gene`, `ncv_color` | UTF-8 string |
| `x`, `y`, `z`, `confidence`, `assignment_confidence` | float64 |
| `cluster` | int32 |
| `is_noise` | boolean |

The same conditional-column rules apply. Use legacy output for Xenium
Ranger compatibility.

### `cells.parquet`

`C` rows, with the same columns and order as `segmentation_cell_stats.csv`:
`cell` is UTF-8 string; all statistics are float64.

### `cell_boundaries.parquet`

One row per joined 2D cell polygon (pooled across z for 3D data):

| Column | Type | Contents |
| --- | --- | --- |
| `cell` | UTF-8 string | Cell name. |
| `n_vertices` | int32 | Vertices, excluding the repeated closing vertex. |
| `geometry` | binary | WKB polygon. |

GeoParquet metadata sets the primary geometry column to `geometry`, encoding
`WKB`, geometry type `Polygon`, and no CRS.

### `cell_boundaries_3d.parquet`

3D runs only: one row per `(cell, layer)` polygon. Column order is `cell`,
`layer`, `n_vertices`, `geometry`; `layer` is UTF-8 string. Other types and
GeoParquet metadata match `cell_boundaries.parquet`.

### `feature_matrix.h5`

10x-style sparse CSC matrix, genes × cells (`nnz` = nonzero entries):

| HDF5 path | Type and shape | Contents |
| --- | --- | --- |
| `/matrix/barcodes` | UTF-8 strings, length `C` | Cell names. |
| `/matrix/data` | int32, length `nnz` | Counts. |
| `/matrix/indices` | int32, length `nnz` | Gene-row indices. |
| `/matrix/indptr` | int64, length `C + 1` | Column pointers. |
| `/matrix/shape` | int64, length `2` | `[G, C]`. |
| `/matrix/features/id`, `/matrix/features/name` | UTF-8 strings, length `G` | Gene names in both arrays. |
| `/matrix/features/feature_type` | UTF-8 strings, length `G` | `Gene Expression`. |
| `/matrix/features/genome` | UTF-8 strings, length `G` | Empty strings. |

### `run_params.toml`

Same config dump and limitations as `segmentation_params.dump.toml`.

### `run.log`

Text log for Parquet output.

### `diagnostic_report.html` / `segmentation_plot.html`

Same reports as in legacy output; written only with `--plot`.
