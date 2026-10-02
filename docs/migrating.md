# Migrating from Julia

The C++ implementation uses the same Baysor method but ships as a native
binary. CLI syntax and some output formats differ; results are not promised
to be identical to Julia. The old documentation is archived at
[0.7.1 (Julia)](https://kharchenkolab.github.io/Baysor/0.7.1/).

## Installation

Use a [release binary or Docker](installation.md); Julia is not required.
If you need the last Julia implementation instead, pin its revision:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="v0.7.1"))
Pkg.build("Baysor")
```

Follow the [archived installation guide](https://kharchenkolab.github.io/Baysor/0.7.1/installation/)
for Julia requirements. The default branch is C++, so an unpinned Julia
package install is no longer appropriate.

## CLI changes

The subcommands remain `baysor run`, `baysor preview` and `baysor segfree`.
Start with [Cell segmentation](run.md) for a working command.

- Set a positive `-m` on the CLI or in a config; there is no fallback value.
- For segmentation, supply `-s` / `--scale` or a usable prior. Scale is not
  inferred from `-m` alone.
- Replace dotted Julia-era overrides such as
  `--config.data.exclude-genes='Blank*'` with flat flags such as
  `--exclude-genes 'Blank*'`, or edit the TOML config passed with `-c`.
- Use `--threads` or `OMP_NUM_THREADS`, not `JULIA_NUM_THREADS`.
- TIFF masks and `:column_name` priors remain supported; MATLAB `.mat` masks
  are not. CSV / Parquet [boundary tables](priors.md) are also accepted.
- Compartment options (`--nuclei-genes` / `--cyto-genes`) are not implemented
  in C++; setting either makes `run` exit with an error.

New controls include Parquet output, molecule-clustering alternatives,
quality filtering, coordinate crops and convergence tolerance. See the
[CLI reference](cli.md) rather than translating old commands flag by flag.

## Config files

Most existing configs can be reused, but check the
[C++ key reference](configuration.md#config-key-reference) for behavior:

- `[data]` is an alias for `[molecules]`.
- Prefer `[prior]` for prior settings. The old
  `[segmentation].unassigned_prior_label` and
  `[segmentation].estimate_scale_from_centers` remain accepted.
- `ncv_method` and `min_pixels_per_cell` are accepted but do not affect the
  C++ pipeline. `max_plot_size` controls molecule-image size in HTML reports.
- `max_z_slices` (default `10`) controls the layer count for 3D polygons.

CLI flags still override config values. Output-directory and format choices
are CLI-only, not config keys.

## Outputs

`-o` now names a directory (default `segmentation`), not a Julia-era file
prefix. Other changes:

- HTML reports are named `diagnostic_report.html` and `segmentation_plot.html`
  (Julia: `segmentation_diagnostics.html` and `segmentation_borders.html`).
- In `segmentation.csv`, noise has cell name `0`, not an empty string.
  `is_noise` is `true` / `false` when transcript IDs are retained, and `1` / `0`
  otherwise.
- Loom count matrices use genes × cells orientation.
- `--output-style parquet` writes Parquet / GeoParquet tables and a 10x-style
  HDF5 matrix. Use legacy output for [Xenium Ranger](xenium.md).

See [Outputs](outputs.md) for the file list and
[Output file formats](output_files.md) for exact schemas.
