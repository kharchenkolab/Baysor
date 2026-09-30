# Migrating from Baysor.jl (v0.7.x)

This site documents the **C++ line** of Baysor. The Julia implementation
(Baysor.jl v0.7.x) is no longer developed; its documentation is archived as
the version **0.7.1 (Julia)** in the version selector and at
[0.7.1/](https://kharchenkolab.github.io/Baysor/0.7.1/). Old links under
`/dev/...` redirect there.

The C++ line follows the v0.7.x segmentation algorithm, so results are
directly comparable. What changes is the packaging, the CLI surface, and some
outputs.

## Installation

- Baysor.jl was installed as a Julia package (`Pkg.add` + `Pkg.build`, or a
  Julia-based Docker image).
- The C++ line is a single native binary: [release
  binaries](installation.md#release-binaries) for Linux x86-64, macOS arm64,
  and Windows x86-64, a [source build](installation.md#building-from-source)
  via `./configure.sh`, or [Docker](installation.md#docker).

## CLI

The three subcommands keep their names: `baysor run`, `baysor preview`,
`baysor segfree`.

- **No more dotted config overrides.** Julia-era commands like
  `--config.segmentation.nuclei-genes=Neat1` or
  `--config.data.exclude-genes='Blank*'` are gone. Use the flat CLI flags
  (`--exclude-genes`, …) or a TOML config file (`-c`).
- **`min-molecules-per-cell` (`-m`) must be set**, on the CLI or in the
  config; the C++ CLI does not fall back to a default for it.
- **Scale requirement is unchanged**: provide `--scale` or a prior
  segmentation from which the scale is estimated.
- **Prior masks**: TIFF masks (binary or integer-labeled) are supported;
  MATLAB `.mat` masks are not. Boundary priors can now be passed as
  CSV/Parquet vertex tables (`vertex_x`, `vertex_y`, `label_id`/`cell_id`).
  Column priors (`:column_name`) work as before.
- **New options** in the C++ CLI include `--output-style legacy|parquet`,
  `--polygon-format`, `--count-matrix-format`, `--cluster-method
  mrf|louvain|leiden|none` with `--cluster-resolution`, `--cluster-graph-k`,
  `--cluster-n-dims`, `--cluster-basis-sample-size`, a `--tol` convergence
  criterion (Julia always ran exactly `--iters` iterations), `--min-qv` and
  `--qv-column` for quality filtering, coordinate crop flags
  (`--x-min`, …), `--force-2d`, and `--skip-ncv-color`.
- **Not yet ported**: compartment segmentation (`--nuclei-genes` /
  `--cyto-genes`) exists in the CLI but is not implemented in the C++ line
  yet; setting these options makes `run` exit with an error.

## Config files

Existing config files keep working: `[data]`, `[segmentation]`, and
`[plotting]` keys are unchanged (`[data]` is an alias for the preferred
`[molecules]` section). New configs should put prior settings in the
`[prior]` section; `[segmentation].unassigned_prior_label` and
`[segmentation].estimate_scale_from_centers` remain accepted for
compatibility. `ncv_method`, `min_pixels_per_cell`, and `max_plot_size` are
accepted but currently unused by the C++ pipeline.

## Outputs

- The default `--output` directory is `segmentation` (it was `segmentation.csv`
  in early C++ builds and a file prefix in Julia).
- Renamed HTML outputs: `segmentation_diagnostics.html` →
  `diagnostic_report.html`, `segmentation_borders.html` →
  `segmentation_plot.html`.
- In `segmentation.csv`, unassigned molecules have cell name `0` (Julia wrote
  an empty string), and `is_noise` is `true`/`false` for Xenium-origin inputs
  (`1`/`0` otherwise).
- The Loom count matrix orientation is fixed to the Loom spec: `/matrix` is
  genes × cells.
- A whole new `parquet` bundle (`--output-style parquet`) with
  Parquet/GeoParquet tables and a 10x-style HDF5 feature matrix is available.

Everything else — file names of the `legacy` bundle, column semantics, and
parameter meaning — matches the Julia release; see
[Output files](output_files.md) for the exact contracts.
