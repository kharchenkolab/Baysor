# Changelog

All notable changes to the C++ line of Baysor are documented here.

## Unreleased

### Added

- Prebuilt binaries for every published GitHub release, built by the
  `release` GitHub Actions workflow: `baysor-<version>-linux-x86_64.tar.gz`,
  `baysor-<version>-macos-arm64.tar.gz`, `baysor-<version>-windows-x86_64.zip`
  and `SHA256SUMS`. They need no extra packages and run on any CPU of their
  architecture: Linux x86_64 with glibc 2.28 or newer, macOS 12 or newer on
  Apple silicon, 64-bit Windows 10 or newer.
- `baysor --version` prints the version.
- `packaging/`: the release build scripts; the Linux binary can be rebuilt
  locally in Docker with `packaging/linux/build-in-docker.sh`. The release
  procedure is described in `RELEASING.md`.
- A versioned documentation site (MkDocs + Material, versioned with `mike`):
  `mkdocs.yml`, rewritten `docs/` pages, and the `docs` GitHub workflow that
  builds the site strictly on docs changes and deploys one site version per
  GitHub release. The `latest` alias only points at the newest stable
  release: pre-releases, backport releases, and `workflow_dispatch` redeploys
  of older tags are deployed without moving `latest`.
- Release-binary installation documentation (Linux x86-64, macOS arm64,
  Windows x86-64 archives with `SHA256SUMS`), plus a Docker section on the
  installation page.
- A "Migrating from Baysor.jl (v0.7.x)" page and developer docs (source
  builds, tests, coverage, benchmarks, releasing).
- `docs/tools/check_cli_docs.py`, which fails when the docs mention a CLI
  option or config key that the sources do not define.
- `docs/tools/migrate_gh_pages.py`, a one-time maintainer migration for the
  `gh-pages` branch: archives the Julia site as `0.7.1 (Julia)` and keeps old
  `/dev/...` links working via redirect stubs.
- A dry-run mode for the `release` workflow (`workflow_dispatch` with
  `dry_run=true` plus optional `ref` and `platforms` inputs): builds and
  smoke-tests the release archives for any branch or commit without needing a
  release and without uploading to one, so the release build can be verified
  before tagging. On build failure the workflow uploads vcpkg's per-port
  build logs as the `vcpkg-logs-<platform>` artifact.

### Changed

- The benchmark suite moved out of this repository into
  [baysor-benchmarks](https://github.com/VPetukhov/baysor-benchmarks)
  (former `benchmarks/`, history preserved); `docs/development.md` links to
  it. Datasets and baselines stay under the local `.bench-data/` directory.
- `README.md` now points to the documentation site and release binaries.
- The documentation pages were rewritten against the C++ implementation
  (required `--min-molecules-per-cell`, actual option defaults, complete
  config-key reference, corrected output-file descriptions).

### Fixed

- Release and CI builds on Windows: the `autoconf2.71` MSYS2 package pinned
  inside vcpkg's gmp port was dropped from the MSYS2 mirrors (404 on all of
  them), breaking every Windows build; `packaging/vcpkg-overlay-ports/gmp`
  backports the upstream vcpkg fix (microsoft/vcpkg#53437, in no release yet).
- Release build on macOS: thrift (an Arrow/Parquet dependency) needs a bison
  newer than the Apple one (2.3) to generate its parser, so the release build
  installs Homebrew's bison.
- CLI help: `--tol` now shows its actual default (`0`), and `preview`/`segfree`
  `-o` is described as an output file rather than a file or directory;
  `configs/example_config.toml` shows the actual `max_plot_size` default
  (3000).


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
