## [0.7.1] — 2024-11-19

### Added

- `--polygon-format=GeometryCollectionLegacy` format similar to `GeometryCollection`, but with integer cell IDs. Used for compatibility with Xenium Ranger. This option will be deprecated after Xenium Ranger is updated.

### Fixed

- CLI argument parsing

### Changed

- `polygons.json` now has polygons for all cells in the dataset. *Cells with 1 or 2 points now have polygons with 4 vertices.*
- Polygon vertices are now slightly shifted farther from the cell center to avoid molecules on the border.

## [0.7.0] — 2024-09-13

### Added

- New algorithm for NCV estimation based on Random Indexing. It is used by default now.
- Support of Parquet format for input molecules (can directly use `transcripts.parquet` from Xenium)
- Now, in case of continuous z-stack (like Xenium), 3D polygons are estimated by binning the z-stack into 20 slices.
- Added support for string IDs for prior segmentation. *To specify the label for non-assigned molecules, use `--config.segmentation.unassigned_prior_label`*.
- New plotting configuration options (see `example_config.toml`)

### Removed

- `no-ncv-estimation` was removed, as the NCV algorithm is fast and memory-efficient now
- The Dirichlet sampling was removed, as it slowed-down the algorithm significantly, but didn't improve the segmentation quality
    - *Correspondingly, parameters `new_component_weight` and `new_component_fraction` were removed*

### Fixed

- Improved multithreading
- Optimized the algorithm performance
- Improved algorithm for polygon estimation: fixed bugs and reduced overlaps
- Fixed a bug in the `split` step, improved cell continuity

### Changed

- Cell IDs in polygon GeoJSON are now strings to match segmentation.csv
- The whole polygon.json format was changed from `GeometryColection` to `FeatureCollection` to match [10x format](https://www.10xgenomics.com/support/software/xenium-ranger/1.7/analysis/inputs/XR-input-overview#compat-files).
- Polygons are now saved by default using the [10x FeatureCollection format](https://www.10xgenomics.com/support/software/xenium-ranger/1.7/analysis/inputs/XR-input-overview#compat-files). Parameter `save-polygons` is replaced with `polygon-format=FeatureCollection`. Set it to `GeometryCollection` to save polygons in the format from Baysor v0.6 or to `none` to prevent saving polygons.
- 2D and 3D polygons are now stored in separate files (`polygons_2d/3d.json`)

## [0.6.2] — 2023-08-15

### Fixed

- Fixed wrong parameter interpretation in `neighborhood_count_matrix()`
- Fixed incompatibility with `julia-1.9` ([#78](https://github.com/kharchenkolab/Baysor/issues/78))
- Fixed `boundary_polygons_from_grid` error ([#82](https://github.com/kharchenkolab/Baysor/issues/82))
- Fixed a problem with huge integers in cell ids ([#83](https://github.com/kharchenkolab/Baysor/issues/83))

## [0.6.1] — 2023-05-11

### Fixed

- Fixed a bug with prior segmentation loading
- Fixed a bug with compartment genes

## [0.6.0] — 2023-04-20

### Added

- New output cell QC parameters `avg_assignment_confidence`, `max_cluster_frac` and `lifespan`
- Segmented cells are now saved to loom instead of TSV. To return an old behavior, use `count-matrix-format="tsv"`
- Minimal multi-threading (see README)

### Removed

- `iters` and `n-cells-init` parameters were removed from the CLI shortcuts. To change them, use the config or `--config.segmentation.iters` and `--config.segmentation.n_cells_init` parameters *(see 'Advanced configuration section in the readme')*.

### Changed

- Breaking changes in config file structure and CLI
- Greatly improved responsiveness of the CLI and simplified installation process
- Major refactoring of the code
- Various performance improvements
- Faster and more precise algorithm for estimating boundary polygons. Now each cell has exactly one polygon in the output GeoJSON.
    - *For method details see [Awrangjeb, 2015](https://doi.org/10.1109/IVCNZ.2015.7761536), it's pretty similar.*
    - *Closes [#15](https://github.com/kharchenkolab/Baysor/issues/15) and [#41](https://github.com/kharchenkolab/Baysor/issues/41), potentially also [#32](https://github.com/kharchenkolab/Baysor/issues/32) and [#37](https://github.com/kharchenkolab/Baysor/issues/37).*
- Using sparse PCA for NCV estimation on large datasets
- `baysor segfree` output is now fully compatible with loom v3 format
- Cells and NCVs now have IDs in the format `{type}{run_id}-{cell_id}`, where `type` is `C` for cells and `V` for NCVs, and `run_id` is a unique ID of Baysor run
- `--save-polygons` now works regardless of `-p`

## [0.5.2] — 2022-06-29

- Fixed some package versions, dependencies should cause fewer bugs now
- Fixed some bugs
- Fixed random seed for all CLI runs
- Adjusted some CLI parameters
- Added `segfree` run option to extract NCV vectors

## [0.5.1] — 2021-12-01

### Changed

- Slightly optimized compilation time
- Minor updates of core formulas

### Added

- CLI parameter `no-ncv-estimation` to disable estimation of NCVs

## [0.5.0] — 2021-06-03

### Changed

- `scale-std` now can be specified from CLI parameters
- Several bugs fixed
- Allow missing genes in the input data
- Updates in the core algorithm
- All diagnostic plots were updated

### Added

- `exclude-genes` option that removes genes from the data frame before segmentation.
- **Segmentation of compartments based on the list of compartment-specific genes**
- Using information about compartment per molecule in the segmentation algorithm when available
- **3D segmentation**

### Removed

- The data can not be split by frames anymore. *This functionality didn't work well previously and was hard to maintan.*

## [0.4.3] — 2021-04-06

### Changed

- Fixed the CLI installation bug
- More information in logging
- Better initialization for molecule clustering. *This should improve cell separability a lot!*
- Some memory optimization
- Fixed plotting performance by moving to Makie.jl

### Added

- Estimating scale when prior segmentation is provided as a CSV column
- Added the option `--save-polygons=GeoJSON` to save cell boundary polygons in the GeoJSON format

## [0.4.2] — 2020-11-26

### Changed

- Fixed Makefile julia version
- Improved polygon visualization
- Regressed to Plots 1.6.0 because of the performance issues
- Fixed docker build
- Small bug fixes
- `min-molecules-per-segment` parameter is working now
- Fixed visualization of prior segmentation

## [0.4.1] — 2020-10-30

### Added

- Added travis config
- Saving NCV colors to the `ncv_color` field of *segmentation.csv*

### Removed

- Dropped support for julia < 1.5

### Changed

- Updated dependencies
- `find_grid_point_labels_kde` now preserves label ids
- Fixed docker build
- Added ImageMagick dependency to fix problems with DAPI prior
