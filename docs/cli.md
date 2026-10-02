# Command-line reference

Start with [Cell segmentation](run.md), [Preview](preview.md) or
[Segmentation-free analysis](segfree.md) for typical commands. These tables
list the CLI options and defaults without a config file. Config values
become defaults; explicit CLI flags override them.

```bash
baysor run --help
baysor preview --help
baysor segfree --help
```

## Input and threads (all subcommands)

Each subcommand takes one required input: a CSV / Parquet molecule table or
Xenium `experiment.xenium` manifest. See [Input data](inputs.md).

| Option | Default | Meaning |
| --- | --- | --- |
| `-c, --config` | — | TOML file; use the space-separated form `-c config.toml`. |
| `-m, --min-molecules-per-cell` | — | Minimum molecules expected in a real cell. Must be positive (CLI or config); drives derived defaults. |
| `-x, --x-column` | `x` | x-coordinate column. |
| `-y, --y-column` | `y` | y-coordinate column. |
| `-z, --z-column` | `z` | Optional z-coordinate column. |
| `-g, --gene-column` | `gene` | Gene-name column. |
| `--qv-column` | `qv` | Quality-value column used by `--min-qv`. |
| `--min-qv` | `-1` | Drop molecules below this quality value; disabled unless set to a nonnegative value. |
| `--force-2d` | off | Ignore the z column. |
| `--x-min`, `--x-max` | −∞, +∞ | Keep molecules within this x range. |
| `--y-min`, `--y-max` | −∞, +∞ | Keep molecules within this y range. |
| `--z-min`, `--z-max` | −∞, +∞ | Keep molecules within this z range. |
| `-t, --threads` | auto | Worker threads; `OMP_NUM_THREADS`, then physical CPU cores. |

`--help` lists a subcommand's options. `baysor --version` prints the release
version.

## Cell segmentation (`run`)

The optional second input is a [prior segmentation](priors.md): a TIFF mask,
a boundary table or `:column_name` from the molecule table.

### Filtering

| Option | Default | Meaning |
| --- | --- | --- |
| `--min-molecules-per-gene` | `1` | Drop genes with fewer molecules. |
| `--exclude-genes` | empty | Comma-separated gene names or glob patterns (`*`, `?`); quote patterns, e.g. `--exclude-genes 'Blank*,MALAT1'`. |

These flags are `run`-only. The corresponding [config keys](configuration.md#molecules-data)
also apply to `preview` and `segfree`.

### Scale and convergence

| Option | Default | Meaning |
| --- | --- | --- |
| `-s, --scale` | — | Approximate cell radius. A positive scale disables scale estimation from the prior. |
| `--scale-std` | `25%` | Scale variation across cells: absolute value or percentage of scale. Replaced by the prior estimate when scale is estimated from a prior. |
| `--iters` | `500` | Maximum segmentation iterations. |
| `--tol` | `0` | Stop when fewer than this fraction of molecules change assignment over 20 consecutive iterations; `0` runs all iterations. |
| `--n-cells-init` | auto | Initial cell count; derived from molecule count and `-m`, reduced when a prior is available. |

### Clustering

| Option | Default | Meaning |
| --- | --- | --- |
| `--cluster-method` | `mrf` | `mrf`, `louvain`, `leiden` or `none`; `ica_mrf` is a legacy alias for `mrf`. |
| `--n-clusters` | auto | Exact count for `mrf` (default `4`); merged target for `louvain` / `leiden` (default `10`). |
| `--cluster-resolution` | `1.0` | Advanced overclustering resolution for graph methods. |
| `--cluster-graph-k` | `15` | NCV neighbors for graph clustering and NCV UMAPs. |
| `--cluster-n-dims` | `20` | NCV dimensions used by graph methods. |
| `--cluster-basis-sample-size` | `100000` | Maximum NCV basis anchors for graph methods. |

### Prior handling

| Option | Default | Meaning |
| --- | --- | --- |
| `--prior-segmentation-confidence` | `0.2` | Prior confidence in `[0, 1]`; see [Prior behavior](priors.md#prior-behavior). |
| `--unassigned-prior-label` | `0` | Column-prior value meaning no assignment. |

### Output

| Option | Default | Meaning |
| --- | --- | --- |
| `-o, --output` | `segmentation` | Output directory. |
| `--output-style` | `legacy` | `legacy` or `parquet`; see [Outputs](outputs.md). |
| `--polygon-format` | `FeatureCollection` | `FeatureCollection`, `GeometryCollection`, `GeometryCollectionLegacy` or `none` (case-insensitive). Legacy output only. |
| `--count-matrix-format` | `loom` | `loom` or `tsv`. Legacy output only. |
| `-p, --plot` | off | Write `diagnostic_report.html` and `segmentation_plot.html`. |
| `--skip-ncv-color` | off | Skip NCV color embedding; omits `ncv_color` from the molecule table. |

### Not implemented

`--nuclei-genes` and `--cyto-genes` accept comma-separated compartment-specific
genes (default empty), but compartment segmentation is not implemented in
C++. Setting either option makes `run` exit with an error.

## Dataset preview (`preview`)

In addition to the [shared options](#input-and-threads-all-subcommands):

| Option | Default | Meaning |
| --- | --- | --- |
| `-o, --output` | `preview.html` | Output HTML file, not a directory. |

## Segmentation-free analysis (`segfree`)

In addition to the [shared options](#input-and-threads-all-subcommands):

| Option | Default | Meaning |
| --- | --- | --- |
| `-k, --k-neighbors` | auto | Neighbors per NCV; `max(n_genes / 10, min_molecules_per_cell, 3)` with integer division. |
| `-o, --output` | `ncvs.loom` | Output Loom file, not a directory. |
