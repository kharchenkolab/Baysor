# Cell segmentation

```bash
baysor run -m 30 -s 8 molecules.csv
```

| Parameter | Meaning |
| --- | --- |
| `-m` | Minimum molecules expected in a real cell. Required on the CLI or in a config. |
| `-s` / `--scale` | Approximate cell radius, in coordinate units. Or pass a [prior](priors.md) as the second input to estimate it; `--prior-segmentation-confidence` controls trust in the prior (default `0.2`). |
| `-c` | [TOML config](configuration.md); CLI flags override its values. |
| `-o` | Output directory (default `segmentation`). |
| `--threads` / `-t` | Worker threads; physical CPU cores by default. |

The input is a CSV / Parquet table with `x`, `y` and `gene` columns; optional
`z` enables 3D segmentation. See [Input data](inputs.md) for other column
names and [Xenium](xenium.md) for `experiment.xenium` inputs. Need the binary?
See [Installation](installation.md).

## Choose `-m` and scale

`30` and `8` above are examples, not defaults. `-m` depends on how many
molecules your protocol detects per cell: sparse protocols need smaller
values than dense ones. It sets neighborhood and noise-estimation defaults,
not a strict filter on the size of output cells. To remove small cells after
segmentation, use `n_transcripts` in the [cell statistics](output_files.md#segmentation_cell_statscsv).

Scale is the most sensitive parameter. Use an approximate cell radius in the
same units as the coordinates (e.g. microns or pixels), not the diameter.
A poor choice can lead to over- or under-segmentation. `--scale-std` describes
variation across cells: the default is `25%` of scale, or you can pass an
absolute value.

A run needs a positive `-m` and either an explicit scale or a usable prior.
Use a [preview](preview.md) or a small crop to check the data and parameters
before a full run.

## Use a prior segmentation

If the table contains prior cell labels in `cell_id`:

```bash
baysor run -m 30 --prior-segmentation-confidence 0.5 \
  -o out --threads 8 molecules.csv :cell_id
```

The prior can also be a TIFF mask or a CSV / Parquet boundary table. Without
`-s`, Baysor estimates scale from the prior. Confidence `0` removes the
prior's assignment penalty; `1` keeps each retained prior segment together,
but allows cells to merge or grow. See [Prior segmentation](priors.md) for
formats, unassigned labels and the confidence-1 guarantee.

## Use a config

Save your settings in [a config file](configuration.md), then run:

```bash
baysor run -c config.toml -o out --threads 8 molecules.csv
```

The config must provide `min_molecules_per_cell` and a scale or prior. Protocol
[presets](configuration.md#protocol-presets) provide starting points; they
are downloaded separately from the binary.

## Inspect the results

The default [output bundle](outputs.md) contains molecule assignments
(`segmentation.csv`), cell statistics (`segmentation_cell_stats.csv`), a count
matrix (`segmentation_counts.loom`) and cell polygons. Add `-p` / `--plot` to
write `diagnostic_report.html` and `segmentation_plot.html` for quality checks.
Use `--output-style parquet` for Parquet tables and a 10x-style HDF5 matrix.

## Clustering methods

Molecule clustering supplies a coarse cell-type prior; it is not the final
cell segmentation. The default is `mrf` with `4` clusters. Alternatives are
`louvain`, `leiden` and `none` (no clustering prior).

For very large, high-gene-panel runs, try `--cluster-method louvain
--n-clusters 10`. For graph methods, the cluster count is a target after
merging communities; for `mrf` it is the fitted count. See the
[CLI reference](cli.md#clustering) for advanced controls.

Compartment segmentation (`--nuclei-genes` / `--cyto-genes`) is not implemented
in C++; setting either option makes `run` exit with an error.

## Threading

`run`, `preview` and `segfree` use the same `--threads` / `-t` option and
`threads` config key. An explicit CLI value takes precedence over the config;
auto (`0`) uses `OMP_NUM_THREADS`, then the number of physical CPU cores.
The effective count is logged at startup; some steps remain serial.

Repeated runs are deterministic for the same data, parameters and thread
count. Multi-threaded outputs agree across multi-threaded counts, but a
single-threaded run can differ; see the [benchmarks](performance/benchmarks.md#accuracy-and-reproducibility).

For all options and defaults, see the [CLI reference](cli.md) or run
`baysor run --help`.
