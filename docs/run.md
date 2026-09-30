# Cell segmentation (`baysor run`)

The `run` subcommand segments molecules into cells:

```bash
baysor run [OPTIONS] coordinates [prior_segmentation]
```

## Minimal forms

With an explicit scale and no prior:

```bash
baysor run -m 30 --scale 8 -o out molecules.csv
```

With a prior segmentation:

```bash
baysor run -m 30 -o out molecules.csv segmentation_mask.tif
```

With a config file (config values become defaults; CLI flags override them):

```bash
baysor run -c configs/xenium.toml -o out data/transcripts.parquet :cell_id
```

Two inputs are always required in some form:

- `min_molecules_per_cell` (`-m`) must be positive — set it on the CLI or in
  the config file;
- the cell scale must be known: either pass `--scale`, or pass a
  [prior segmentation](priors.md) from which the scale is estimated.

## Positionals

`coordinates` (required)
: Molecule table in CSV or Parquet form, or a Xenium
  `experiment.xenium` manifest (see [Input data](inputs.md)).

`prior_segmentation` (optional)
: One of an image mask (TIFF), a boundary table (CSV/Parquet), or
  `:column_name` for transcript-native prior labels in the molecule table.
  See [Prior segmentation](priors.md).

## Options

### Input columns and filtering

| Option | Default | Description |
| --- | --- | --- |
| `-c, --config` | — | TOML file with configuration (see [Configuration](configuration.md)) |
| `-x, --x-column` | `x` | Name of the x column |
| `-y, --y-column` | `y` | Name of the y column |
| `-z, --z-column` | `z` | Name of the z column |
| `-g, --gene-column` | `gene` | Name of the gene column |
| `--qv-column` | `qv` | Name of the quality-value column used by `--min-qv` |
| `-m, --min-molecules-per-cell` | — | Minimal number of molecules for a cell to be considered real. Required (CLI or config); drives several derived defaults |
| `--min-molecules-per-gene` | `1` | Minimal number of molecules per gene; genes below this are dropped |
| `--exclude-genes` | — | Comma-separated gene names or glob patterns (`*`, `?`) to drop, e.g. `Blank*,MALAT1` |
| `--min-qv` | `-1` | Drop molecules with quality value below this threshold (disabled unless set) |
| `--force-2d` | off | Ignore the z column in the data |
| `--x-min`, `--x-max` | ±∞ | Keep only molecules within this x range |
| `--y-min`, `--y-max` | ±∞ | Keep only molecules within this y range |
| `--z-min`, `--z-max` | ±∞ | Keep only molecules within this z range |

### Segmentation scale and convergence

| Option | Default | Description |
| --- | --- | --- |
| `-s, --scale` | — | Approximate cell radius. Explicit `--scale` disables scale estimation from the prior (`estimate_scale_from_prior`) |
| `--scale-std` | `25%` | Std of scale across cells: an absolute number, or `N%` of `--scale` |
| `--iters` | `500` | Maximum number of algorithm iterations |
| `--tol` | `0` | Convergence tolerance: stop once fewer than `tol` fraction of molecules change assignment over 20 consecutive iterations. `0` runs all `--iters` iterations |
| `--n-cells-init` | auto | Initial number of cells. Auto = `2 * n_molecules / min_molecules_per_cell`, reduced when a prior segmentation is present |

### Molecule clustering prior

| Option | Default | Description |
| --- | --- | --- |
| `--cluster-method` | `mrf` | Clustering prior: `mrf`, `louvain`, `leiden`, or `none` (legacy alias `ica_mrf` = `mrf`) |
| `--n-clusters` | auto | Target number of molecule clusters / major cell types: exact for `mrf` (default `4`), merged target for `louvain`/`leiden` (default `10`) |
| `--cluster-resolution` | `1.0` | Advanced overclustering resolution for `louvain`/`leiden` |
| `--cluster-graph-k` | `15` | NCV nearest neighbors used for graph clustering and NCV UMAPs |
| `--cluster-n-dims` | `20` | NCV dimensions used by `louvain`/`leiden` |
| `--cluster-basis-sample-size` | `100000` | Maximum number of basis anchors used by `louvain`/`leiden` |

### Prior handling

| Option | Default | Description |
| --- | --- | --- |
| `--prior-segmentation-confidence` | `0.2` | Confidence of the prior segmentation, in `[0, 1]` (see [Prior segmentation](priors.md)) |
| `--unassigned-prior-label` | `0` | Label for unassigned cells in transcript-native prior segmentation |

### Output

| Option | Default | Description |
| --- | --- | --- |
| `-o, --output` | `segmentation` | Output directory |
| `--output-style` | `legacy` | Output bundle style: `legacy` or `parquet` (see [Outputs](outputs.md)) |
| `--polygon-format` | `FeatureCollection` | GeoJSON root type for polygon output: `FeatureCollection`, `GeometryCollection`, or `none`. `legacy` style only |
| `--count-matrix-format` | `loom` | Count matrix format: `loom` or `tsv`. `legacy` style only |
| `-p, --plot` | off | Also write `diagnostic_report.html` and `segmentation_plot.html` |
| `--skip-ncv-color` | off | Skip the neighborhood-composition color embedding to speed up development runs |

`--nuclei-genes` and `--cyto-genes` exist in the CLI and config for
compartment-aware segmentation of intracellular structure, but compartment
segmentation is **not yet implemented** in the C++ line: setting either option
makes `run` exit with an error.

## Clustering methods

`run` supports four molecule-clustering modes, which provide a coarse
compatibility prior for segmentation:

- `mrf` — the legacy Baysor molecule clustering path; the default. `--n-clusters`
  is the exact number of molecule clusters to fit (default `4`).
- `louvain` — Louvain clustering of NCV basis anchors followed by label
  transfer to all molecules. `--n-clusters` is the target final coarse cluster
  count after anchor communities are merged (default `10`).
- `leiden` — the same anchor-based NCV workflow with Leiden refinement.
- `none` — disables the clustering prior entirely.

```bash
baysor run --cluster-method mrf --n-clusters 4 ...
baysor run --cluster-method louvain --n-clusters 10 ...
```

`--cluster-resolution`, `--cluster-graph-k`, `--cluster-n-dims`, and
`--cluster-basis-sample-size` are advanced controls for the graph methods; the
defaults are a reasonable starting point.

For very large runs, especially high-gene-panel Xenium runs such as 5K panels,
prefer Louvain clustering with about 10 coarse clusters:

```bash
baysor run --cluster-method louvain --n-clusters 10 ...
```

This keeps the clustering prior coarse enough for large-scale segmentation
while using the anchor-based NCV graph path instead of the legacy MRF
clustering model.

## Threading

Baysor uses OpenMP for parallel sections. Set the thread count with the
standard OpenMP environment variable:

```bash
OMP_NUM_THREADS=20 baysor run ...
```

Some phases are intentionally serial or only partially parallel, so CPU use
may not stay at the requested thread count for the full run.

## Config files

Protocol presets live in the repository under
[configs/](https://github.com/kharchenkolab/Baysor/tree/HEAD/configs):
[example_config.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/example_config.toml),
[xenium.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/xenium.toml),
[iss.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/iss.toml),
[starmap.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/starmap.toml),
[osm_fish.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/osm_fish.toml).

CLI flags override config values. See [Configuration](configuration.md) for
all config keys.

## Other subcommands

- [Dataset preview](preview.md)
- [Segmentation-free NCVs](segfree.md)

For the up-to-date option list, run `baysor run --help`.
