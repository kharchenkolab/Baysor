# Configuration

All options can be set in a TOML config file passed with `-c/--config`.
Config values become defaults; CLI flags override them. The resolved
parameters are written to the output directory
([`segmentation_params.dump.toml`](output_files.md#segmentation_paramsdumptoml)).

Protocol presets live in the repository under
[configs/](https://github.com/kharchenkolab/Baysor/tree/HEAD/configs):

- [example_config.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/example_config.toml)
- [xenium.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/xenium.toml)
- [iss.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/iss.toml)
- [starmap.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/starmap.toml)
- [osm_fish.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/osm_fish.toml)

```bash
baysor run -c configs/xenium.toml -o out data/transcripts.parquet :cell_id
```

## Value syntax

Values are parsed strictly by type: an unparsable or wrong-typed value aborts
the run with an error naming the key. Strings may be quoted (quotes are
stripped); integers accept digit separators (`1_000`) and integral floats
(`50.0`). Booleans are `true` / `false` (or `1` / `0`). `#` starts a comment.

## `[molecules]` / `[data]`

`[molecules]` is the preferred section name; `[data]` is accepted as an alias
for existing configs. If both are present, `[data]` is applied last.

| Key | Type | Default | Description |
| --- | --- | --- | --- |
| `x` | string | `x` | Name of the x column |
| `y` | string | `y` | Name of the y column |
| `z` | string | `z` | Name of the z column |
| `gene` | string | `gene` | Name of the gene column |
| `qv` | string | `qv` | Name of the quality-value column used by `min_qv` |
| `force_2d` | bool | `false` | Ignore the z column in the data |
| `min_molecules_per_gene` | int | `1` | Minimal number of molecules per gene |
| `exclude_genes` | string | empty | Comma-separated gene names or glob patterns (`*`, `?`) to drop |
| `min_molecules_per_cell` | int | `0` | Minimal number of molecules for a cell to be considered real; must be set (CLI or config) |
| `confidence_nn_id` | int | `0` | Nearest neighbors for confidence estimation. `0` = `max(min_molecules_per_cell / 2 + 1, 5)` |
| `min_qv` | float | `-1` | Drop molecules with quality value below this threshold; `-1` disables |
| `x_min`, `x_max` | float | ±∞ | Coordinate crop in x |
| `y_min`, `y_max` | float | ±∞ | Coordinate crop in y |
| `z_min`, `z_max` | float | ±∞ | Coordinate crop in z |
| `min_molecules_per_segment` | int | `0` | Minimal molecules for a prior segment to be kept. `0` = `max(min_molecules_per_cell / 4, 2)`. (Prior setting, kept in this section for compatibility; `[prior]` is preferred.) |

## `[segmentation]`

| Key | Type | Default | Description |
| --- | --- | --- | --- |
| `scale` | float | `-1` | Approximate cell radius; negative = estimate from the prior or `min_molecules_per_cell` |
| `scale_std` | string | `"25%"` | Std of scale across cells: absolute number, or `N%` of `scale` |
| `cluster_method` | string | `mrf` | Molecule clustering prior: `mrf`, `louvain`, `leiden`, `none` (legacy alias `ica_mrf`) |
| `n_clusters` | int | `0` | Number of molecule clusters / cell types; `0` = `4` for `mrf`, `10` for `louvain`/`leiden` |
| `cluster_resolution` | float | `1.0` | Overclustering resolution for `louvain`/`leiden` |
| `cluster_graph_k` | int | `15` | NCV neighbors for graph clustering and NCV UMAPs |
| `cluster_n_dims` | int | `20` | NCV dimensions used by `louvain`/`leiden` |
| `cluster_basis_sample_size` | int | `100000` | Maximum basis anchors used by `louvain`/`leiden` |
| `prior_segmentation_confidence` | float | `0.2` | Confidence of the prior segmentation, in `[0, 1]` |
| `iters` | int | `500` | Maximum number of algorithm iterations |
| `tol` | float | `0` | Convergence tolerance (fraction of molecules changing assignment over 20 iterations); `0` runs all `iters` |
| `n_cells_init` | int | `0` | Initial number of cells; `0` = auto (`2 * n_molecules / min_molecules_per_cell`, prior-aware) |
| `nuclei_genes` | string | empty | Nuclei-specific genes. **Not yet implemented** in the C++ line: setting it makes `run` exit with an error |
| `cyto_genes` | string | empty | Cytoplasm-specific genes. **Not yet implemented** in the C++ line |
| `estimate_scale_from_centers` | bool | `true` | Compatibility key for `[prior].estimate_scale_from_prior` |
| `unassigned_prior_label` | string | `0` | Compatibility key for `[prior].unassigned_label` |

## `[prior]`

Prior input specifics; see [Prior segmentation](priors.md) for semantics.

| Key | Type | Default | Description |
| --- | --- | --- | --- |
| `type` | string | `none` | Prior input type: `none`, `column`, `image`, `boundary` |
| `path` | string | empty | Path to the mask / boundary table |
| `column_name` | string | empty | Column name for `type = "column"` |
| `unassigned_label` | string | `0` | Label meaning "no prior assignment" for column priors |
| `min_molecules_per_segment` | int | `0` | Minimal molecules for a prior segment to be kept; `0` = `max(min_molecules_per_cell / 4, 2)` |
| `estimate_scale_from_prior` | bool | `true` | Estimate scale and scale_std from the prior segments |

## `[plotting]`

| Key | Type | Default | Description |
| --- | --- | --- | --- |
| `gene_composition_neighborhood` | int | `0` | Spatial neighborhood size used to compute NCVs. `0` = `max(n_genes / 10, min_molecules_per_cell, 3)` |
| `gene_composition_neigborhood` | int | — | Julia-era misspelling of the key above, still accepted |
| `min_pixels_per_cell` | int | `15` | Accepted for compatibility with Julia-era configs; currently unused by the C++ pipeline |
| `max_plot_size` | int | `3000` | Accepted for compatibility with Julia-era configs; currently unused by the C++ pipeline |
| `ncv_method` | string | `ri` | Accepted for compatibility (`ri`/`dense`/`sparse`); the C++ pipeline currently always uses random indexing |

## Clustering options

The pre-segmentation molecule clustering is configured through `[segmentation]`
(see also [run](run.md#clustering-methods)):

- `cluster_method = "mrf" | "louvain" | "leiden" | "none"`
- `n_clusters` — exact cluster count for `mrf`, merged target for `louvain`
  and `leiden`
- `cluster_resolution` — advanced overclustering resolution for `louvain` and
  `leiden`
- `cluster_graph_k` — neighbors in NCV space for the anchor graph and NCV
  UMAPs
- `cluster_n_dims`, `cluster_basis_sample_size` — anchor-based NCV controls

For very large, high-gene-panel Xenium runs (e.g. 5K panels), a good starting
point is:

```toml
[segmentation]
cluster_method = "louvain"
n_clusters = 10
```

Treat `10` as a practical starting point for the final coarse cluster count.

## Prior options

Prefer `[prior]` in new configs:

```toml
[prior]
unassigned_label = "UNASSIGNED"
estimate_scale_from_prior = true
```

`[segmentation].unassigned_prior_label` and
`[segmentation].estimate_scale_from_centers` remain accepted for older
configs. CLI flags (`--unassigned-prior-label`, `--scale`) override both
forms.

## CLI vs config

Config values set defaults; CLI flags override them:

```bash
baysor run \
  -c configs/xenium.toml \
  --x-min 0 --x-max 2000 \
  --y-min 0 --y-max 2000 \
  -o out \
  data/experiment.xenium :cell_id
```

## Protocol notes

### Xenium

[configs/xenium.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/xenium.toml)
maps the Xenium columns (`x_location`, `y_location`, `z_location`,
`feature_name`, `qv`), filters negative-control genes
(`NegControl*,BLANK_*,antisense_*`), sets `min_molecules_per_gene = 10`,
`min_molecules_per_cell = 50`, `prior_segmentation_confidence = 0.5`, and the
`UNASSIGNED` prior label.

### ISS / STARmap / osm-FISH

The preset configs mainly adjust `min_molecules_per_gene`,
`min_molecules_per_cell`, and the expected column layout.

For protocol-specific work, start from the closest preset and override only
the dataset-specific values.
