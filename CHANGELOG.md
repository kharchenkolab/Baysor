# Changelog

All notable changes to the C++ line of Baysor are documented here.

## Unreleased

### Added

- `--threads` / `-t` on `run`, `preview` and `segfree` (and a top-level
  `threads` config key) to set the number of worker threads. Defaults to
  `BAYSOR_NUM_THREADS`, then `OMP_NUM_THREADS` (backward compatibility for
  existing scripts and the benchmark harness), then the number of CPU cores.
  The effective thread count is logged at start-up.

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
- Docker images of every release, built from the release binary by the
  `docker` job of the `release` workflow (`packaging/docker/Dockerfile`, a
  `debian:12-slim` runtime image with the extracted Linux archive,
  non-root user, `WORKDIR /data`) and pushed to GHCR
  (`ghcr.io/<owner>/baysor`; the package is private until made public once in
  its settings) and to Docker Hub (`vpetukhov/baysor`, or the
  `DOCKERHUB_REPOSITORY` repository variable) when the
  `DOCKERHUB_USERNAME`/`DOCKERHUB_TOKEN` secrets exist — otherwise Docker Hub
  is skipped with a warning and GHCR still works. Tags are `X.Y.Z`, plus
  `latest` only for the newest stable release, decided by
  `packaging/is_latest_release.py` with the same rule as the docs site's
  `latest` (unit-tested by `packaging/is_latest_release_test.py`). The image
  is smoke-tested in the container (`--version`, `--help`, a synthetic
  `baysor run`) before any push; dry runs build and smoke-test it without
  pushing. See `RELEASING.md`, "Docker images", and the installation docs.

### Changed

- The molecule graph's CGAL Delaunay triangulation is built once per `run`
  instead of twice (REPORT.md 6.4, ISSUES.md 8b): the confidence estimate
  computes the unfiltered edges once (`compute_molecule_adjacency`, returned
  as `ConfidenceEstimationDetails::adjacency`) and the segmentation graph is
  built from them with `build_molecule_graph(..., precomputed_edges)`, which
  applies the long-edge filter (n_mads=2.0) and produces a bit-identical
  graph to the former rebuild. Slides with exact duplicate coordinates —
  where the two historical `normalize_points()` runs drew different jitter
  batches — fall back to recomputing (`AdjacencyResult::normalize_rng_draws`
  records the batch size), so 1-thread output is bitwise identical in all
  cases. The confidence kNN keeps its results in
  flat `n x k` arrays with an in-place tied-run index sort instead of
  per-molecule vectors and a full `stable_sort` (review-bmm.md C7-1/C7-2),
  and the confidence estimate (and `preview`) read the kth-neighbour distance
  from a block-wise streaming query (`knn_kth_distances`) instead of holding
  the whole-slide `n x k` result (3.0 GiB and 223 M allocations at 10.6M
  molecules).

- The benchmark suite moved out of this repository into
  [baysor-benchmarks](https://github.com/VPetukhov/baysor-benchmarks)
  (former `benchmarks/`, history preserved); `docs/development.md` links to
  it. Datasets and baselines stay under the local `.bench-data/` directory.
- OpenMP is no longer used or required (no `libomp`, `vcomp140.dll`, or
  `-fopenmp` anywhere): all parallelism runs on Baysor's own persistent
  `std::thread` pool (`include/baysor/utils/thread_pool.h`), which also backs
  the FetchContent dependencies (umappp, knncolle, irlba, CppKmeans) through
  subpar's custom-parallelization hooks. Multi-threaded runs are
  deterministic: the E-step RNG streams are keyed by (iteration, chunk), and
  parallel reductions merge in a fixed order, so repeated runs are
  byte-identical and results do not depend on thread scheduling. With
  `--threads 1` results are bitwise identical to the previous OpenMP build.
  Where the Eigen version supports it (>= 3.4.90), Eigen's own GEMM thread
  pool replaces OpenMP for dense matrix products; older Eigen versions run
  dense products single-threaded. The pool wakes workers individually (no
  thundering herd), uses a short bounded spin-then-block on job hand-off and
  region completion (tunable via `BAYSOR_POOL_SPIN_US`, 0 disables), and the
  default thread count is the number of physical CPU cores. The umappp/kNN
  neighborhood-graph construction for NCV color embedding now also runs on
  the pool.
- Molecule clustering scales to large gene panels:
  - The NCV neighbourhood k-NN searches (k = genes / 10, used by
    `--cluster-method louvain|leiden` and the NCV colours) keep the k best
    candidates in a heap instead of an insertion-sorted array; same
    neighbours in the same order. About half the instructions of a whole run
    on an 8,400-gene CosMx crop.
  - These searches run in blocks bounded to ~32 MiB instead of 32,768
    queries, and the neighbourhood count matrix is assembled in place: peak
    RSS of that crop 475 → 233 MiB. Output unchanged.
  - The MRF clustering (`--cluster-method mrf`, the default) computes its
    convergence check inside the parallel E-step and runs the M-step in
    parallel over genes. Output unchanged at any thread count.
  - ICA initialisation of the MRF clustering on panels above 3,000 genes
    builds the gene co-occurrence matrix sparse and computes only the
    `n_clusters` leading eigenvectors (truncated SVD, irlba) instead of a dense
    eigen-decomposition whose cost grows with the cube of the gene count.
    A 150k-molecule, 17,500-gene CosMx whole-transcriptome crop with the
    default method now finishes in 86 s at 8 threads (1.1 GB peak RSS); it
    was stopped after > 70 min at 9.8 GB before.
    Panels up to 3,000 genes are unchanged. Above 3,000 genes the eigenvector
    signs, and with them the FastICA start, can differ from before, so
    clusters and segmentation may change within the usual run-to-run range.

- `run --plot` and `preview` reports are much faster: PNG images are
  encoded with zlib (level 1) instead of stb_image_write and several images
  are encoded in parallel. zlib is now an explicit build dependency (it was
  already required by HDF5 and libtiff).
- `[plotting] max_plot_size` (default 3000) is now honoured: it sets the
  longer side, in pixels, of the molecule images of the `run --plot`
  segmentation report and of the `preview` report (previously always
  6000 px wide), which halves the HTML size of a square dataset.
- Faster BMM segmentation loop with less memory. The E-step accumulates
  neighbour weights per distinct cell before the Julia-order dictionary
  (only molecules with two or more adjacent cells use it), hoists
  loop-invariant `exp`/`pow` factors and evaluates one `exp` per candidate
  instead of two. Grouping molecules by cell, the connected-component split,
  dropping empty cells, the cluster mode of the M-step, the prior-segment
  bookkeeping and applying the E-step result no longer allocate per call
  and run in parallel. Each iteration runs as one persistent parallel region
  of the thread pool (`parallel_region`, with work-shared loops, barriers and
  single blocks) instead of waking the pool for every loop; per-worker
  scratch no longer shares cache lines between workers; barrier waiters
  sleep on a futex on Linux and do not spin when there are more threads than
  physical cores (`BAYSOR_POOL_SPIN_US` still overrides). The assignment
  history keeps the newest entry plus per-iteration changes instead of
  `depth` full copies (about 63 instead of 200 bytes per molecule at the
  default depth), the history vote runs in parallel, and the molecule graph
  is moved into the BMM data instead of being copied. Results are unchanged
  except for possible last-bit differences of the E-step densities from the
  fused `exp`.

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
