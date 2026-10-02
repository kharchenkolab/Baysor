# Changelog

All notable changes to the C++ line of Baysor are documented here.

## [cpp-0.9.0] — 2026-10-02

### Added

- Prebuilt binaries for every release: Linux x86-64 (glibc 2.28+), macOS arm64
  (macOS 12+) and Windows x86-64, with `SHA256SUMS`. No extra packages are needed.
- Docker images for every release on GHCR and Docker Hub (`X.Y.Z`, and `latest`
  for the newest stable release).
- Versioned documentation site with pages rewritten for the C++ implementation,
  a "Migrating from Baysor.jl (v0.7.x)" page and a "Performance" section.
- `-t, --threads` for `run`, `preview` and `segfree` (and the `threads` config
  key). By default Baysor uses `OMP_NUM_THREADS`, then the number of physical
  CPU cores.
- `--version` for `baysor` and every subcommand; `baysor run --version` prints a
  bare version string, as expected by Sopa.
- `--polygon-format GeometryCollectionLegacy` (integer cell ids) for Xenium
  Ranger 3.x.
- `[plotting] max_z_slices` sets the number of z-layers used for 3D polygons
  (default 10, #169).

### Changed

- OpenMP is no longer used: all parallel code runs on Baysor's own thread pool.
  Runs are reproducible: repeated runs give identical results, and all
  multi-threaded runs give the same result whatever the number of threads
  (a single-threaded run uses its own random-number stream and can differ
  slightly).
- Faster segmentation with less memory than cpp-0.8.3: about 2–3× faster at
  6 threads and 20–30 % faster at 1 thread, with 10–40 % lower peak memory.
- Molecule clustering scales to large gene panels: above 3,000 genes the ICA
  initialisation uses a truncated sparse decomposition, which is up to 60×
  faster with 5× less memory. On such panels clusters and segmentation can
  differ slightly from cpp-0.8.3.
- `run --plot` and `preview` reports are faster, and `[plotting] max_plot_size`
  is now honoured.
- The default thread count is the number of physical cores instead of all
  logical CPUs.
- Invalid config values are reported as errors naming the key instead of
  silently falling back to the default. Integer keys accept `50.0`, `1e2` and
  `1_000`, as in Baysor.jl.
- `--polygon-format` is case-insensitive, and unknown values are an error.

### Fixed

- `--prior-segmentation-confidence 1` keeps all molecules of a prior cell in one
  final cell, including duplicated CosMx transcripts (#117).
- Every cell in `segmentation.csv` gets exactly one polygon, so Xenium Ranger
  no longer reports `MissingCellPolygon` (#165).
- Multi-threaded runs reused the same random numbers in every iteration, which
  made results depend on the thread count and differ between repeated runs
  (e.g. about 7 % more cells on an ISS dataset at 6 threads).
- Single-threaded results no longer depend on the output path.
- Crashes on very small inputs: `--n-cells-init 1` or too few molecules for
  `-m`, fewer genes than clusters in the ICA initialisation, and NCV colours
  with very few anchors.
- Gene names stored in Parquet binary-dictionary columns are decoded as text.
- CLI help shows the actual `--tol` default.

## [cpp-0.8.3] — 2026-07-31

### Changed

- Changed the default `baysor run --output` directory from `segmentation.csv`
  to `segmentation`.

### Fixed

- Fixed Loom matrix orientation so `/matrix` rows match `/row_attrs` genes and
  columns match `/col_attrs` cells.

## [cpp-0.8.2] — 2026-04-30

### Changed

- Improved prior-based scale estimation for large prior segmentations by using
  exact KD-tree nearest-neighbor queries instead of an all-pairs center scan.
- Reduced segmentation-loop overhead when `tol = 0` by avoiding unnecessary
  assignment snapshot work.
- Parallelized connected-component splitting in the segmentation loop.
- Parallelized boundary-polygon construction across cells.
- Parallelized graph-clustering resolution attempts for Louvain and Leiden.
- Updated generated CLI help pages and Windows/CMake build support.
- Documented Louvain with about 10 coarse clusters as the recommended starting
  point for very large high-gene-panel Xenium runs.

### Fixed

- Fixed `unassigned_prior_label` config handling and the
  `--unassigned-prior-label` CLI override for transcript-native prior labels.
- Fixed preservation of existing prior options when the positional
  `prior_segmentation` argument is parsed.
- Fixed Windows CI/build issues.

## [cpp-0.8.1] — 2026-04-22

### Added

- `legacy` and `parquet` output styles for `baysor run`.
- A documented Xenium workflow based on `experiment.xenium` input plus `xeniumranger import-segmentation`.
- User-facing documentation pages for installation, CLI usage, inputs, outputs, preview, segfree, configuration, examples, and Xenium workflows.
- A file-level output reference for the current output bundles.
- Additional molecule clustering modes for `baysor run`:
  - `mrf`
  - `louvain`
  - `leiden`
  - `none`

### Changed

- `run`, `preview`, and `segfree` now resolve `experiment.xenium` automatically.
- `legacy` output now emits Xenium Ranger-friendly fields automatically for Xenium-origin inputs.
- NCV color generation in `run` and `preview` now uses an anchor-first streaming path:
  - learn the NCV basis from anchors
  - fit UMAP on sampled anchors
  - stream exact low-dimensional NCV vectors for all molecules
- Default `n_cells_init` is now prior-aware when a prior segmentation is present.
- `segmentation_counts.loom` writing now uses a row-oriented path instead of a late sparse-matrix duplication step.
- For graph clustering methods, `n_clusters` now controls the final coarse cluster count after anchor communities are merged.
- NCV-based clustering, 2D report UMAPs, and 3D NCV color UMAPs now share a consistent separation between:
  - the spatial neighborhood used to compute NCVs
  - the graph neighborhood used in NCV space

### Fixed

- Reduced clustering memory by removing the per-thread dense gene-by-gene correlation matrix.
- Reduced segmentation memory by switching persistent component gene counts away from dense `double` storage.
- Reduced NCV memory by removing global all-molecule neighborhood materialization and the full all-molecule high-dimensional NCV matrix from the `run` and `preview` paths.
- Reduced Loom write time and memory by avoiding an extra sparse-matrix conversion during output.
- Fixed NCV color / report embedding consistency when `cluster_method=louvain` or `leiden`, so report UMAP geometry no longer depends on an unintended clustering-only NCV neighborhood.
- Fixed streamed NCV projection to use the same logged feature transform as anchor NCV fitting, which restores non-collapsed NCV colors in `run` and `preview`.

## [0.8.0] — 2026-04-17

### Added

- Native C++ root build and test layout driven by CMake.
- 3D boundary estimation and 3D polygon output for volumetric datasets such as STARmap.
- Labeled TIFF prior support in the C++ loader.
- Adaptive NCV anchor selection that no longer fails when no molecules exceed a hard `0.95` confidence cutoff.

### Changed

- Repository structure now reflects the C++ implementation as the primary codebase.
- Example READMEs, Dockerfiles, and CI workflows now target the C++ CLI.

### Fixed

- Multiple segmentation-loop parity issues in the C++ implementation, including prior handling, drop thresholds, clustering, and history refinement behavior.
- STARmap NCV postprocessing now degrades gracefully instead of failing or producing empty anchors.
