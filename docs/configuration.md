# Configuration

Save input and segmentation settings in a TOML file and pass it with `-c`.
Explicit CLI flags override the file:

```bash
baysor run -c config.toml -o out --threads 8 molecules.csv
```

A minimal `config.toml` for a table with `x`, `y` and `gene` columns is:

```toml
[molecules]
min_molecules_per_cell = 30

[segmentation]
scale = 8.0
```

These values are examples; choose `min_molecules_per_cell` and scale for your
[dataset](run.md#choose-m-and-scale). Output choices such as `--output-style`
and `--plot` are CLI flags, not config keys. Use the space-separated form
`-c config.toml` or `--config config.toml`.

## Protocol presets

Download a preset separately from the binary, then adjust it for your data:

| Preset | Use |
| --- | --- |
| [xenium.toml](https://github.com/kharchenkolab/Baysor/blob/cpp-0.9.0/configs/xenium.toml) | Xenium column names, control-gene filters, `min_molecules_per_cell = 50`, prior confidence `0.5` and unassigned label `UNASSIGNED`. |
| [iss.toml](https://github.com/kharchenkolab/Baysor/blob/cpp-0.9.0/configs/iss.toml) | Sparse ISS data. |
| [osm_fish.toml](https://github.com/kharchenkolab/Baysor/blob/cpp-0.9.0/configs/osm_fish.toml) | osm-FISH data. |
| [starmap.toml](https://github.com/kharchenkolab/Baysor/blob/cpp-0.9.0/configs/starmap.toml) | STARmap data. |

See the [Xenium workflow](xenium.md) for a download-and-run command, or
[example_config.toml](https://github.com/kharchenkolab/Baysor/blob/cpp-0.9.0/configs/example_config.toml)
for a longer template. The tables below describe the C++ defaults; some
comments in the example template still describe Julia-era behavior.

## Config key reference

Missing keys keep their defaults. Use `key = value` lines in the sections
below; quote strings and use `true` / `false` for booleans. `#` starts a
comment. Invalid numeric or boolean values cause an error; unknown keys are
ignored, so check spelling carefully.

Derived defaults use integer division and `m = max(min_molecules_per_cell, 3)`.

### Top-level keys

Place these before any `[section]`.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `threads` | int | `0` | Worker threads; auto uses `OMP_NUM_THREADS`, then physical CPU cores. |

## `[molecules]` / `[data]`

Use `[molecules]` in new configs. `[data]` is an alias; if both are present,
`[data]` is applied last.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `x`, `y`, `z` | string | `x`, `y`, `z` | Coordinate column names; z is optional. |
| `gene` | string | `gene` | Gene-name column. |
| `qv` | string | `qv` | Quality-value column used by `min_qv`. |
| `force_2d` | bool | `false` | Ignore the z column. |
| `min_molecules_per_gene` | int | `1` | Drop genes with fewer molecules. |
| `exclude_genes` | string | empty | Comma-separated gene names or glob patterns (`*`, `?`). |
| `min_molecules_per_cell` | int | `0` | Minimum molecules expected in a real cell; must be set to a positive value. |
| `confidence_nn_id` | int | `0` | Neighbors for noise-confidence estimation; `0` = `max(m / 2 + 1, 5)`. |
| `min_qv` | float | `-1` | Drop molecules below this quality value; negative disables. |
| `x_min`, `x_max` | float | −∞, +∞ | Crop in x. |
| `y_min`, `y_max` | float | −∞, +∞ | Crop in y. |
| `z_min`, `z_max` | float | −∞, +∞ | Crop in z. |
| `min_molecules_per_segment` | int | `0` | Prior segment threshold; `0` = `max(m / 4, 2)`. Prefer this key in `[prior]`. |

## `[segmentation]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `scale` | float | `-1` | Approximate cell radius; must be positive unless a usable prior supplies it. A positive value disables prior scale estimation. |
| `scale_std` | string | `"25%"` | Scale variation: absolute value or percentage of scale; replaced when scale is estimated from a prior. |
| `cluster_method` | string | `mrf` | `mrf`, `louvain`, `leiden` or `none`; legacy alias `ica_mrf`. |
| `n_clusters` | int | `0` | Auto = `4` for `mrf`, `10` for graph methods; exact count for `mrf`, merged target for `louvain` / `leiden`. |
| `cluster_resolution` | float | `1.0` | Advanced overclustering resolution for graph methods. |
| `cluster_graph_k` | int | `15` | NCV neighbors for graph clustering and NCV UMAPs. |
| `cluster_n_dims` | int | `20` | NCV dimensions used by graph methods. |
| `cluster_basis_sample_size` | int | `100000` | Maximum NCV basis anchors used by graph methods. |
| `prior_segmentation_confidence` | float | `0.2` | Prior confidence in `[0, 1]`; see [Prior segmentation](priors.md). |
| `iters` | int | `500` | Maximum segmentation iterations. |
| `tol` | float | `0` | Stop when fewer than this fraction of molecules change assignment over 20 consecutive iterations; `0` runs all iterations. |
| `n_cells_init` | int | `0` | Auto = `2 * floor(n_molecules / m)`, reduced when a prior is available. |
| `nuclei_genes`, `cyto_genes` | string | empty | Not implemented in C++; setting either makes `run` fail. |
| `estimate_scale_from_centers` | bool | `true` | Compatibility alias for `[prior].estimate_scale_from_prior`. |
| `unassigned_prior_label` | string | `0` | Compatibility alias for `[prior].unassigned_label`. |

## `[prior]`

These keys describe the [prior input](priors.md). The second CLI positional
argument overrides its type and source. `[prior]` keys take precedence over
the compatibility forms above; explicit CLI flags override both.

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `type` | string | `none` | `none`, `column`, `image` or `boundary`. |
| `path` | string | empty | TIFF mask or boundary-table path. |
| `column_name` | string | empty | Column name for a column prior, without the `:` prefix. |
| `unassigned_label` | string | `0` | Column-prior value meaning no assignment. |
| `min_molecules_per_segment` | int | `0` | Discard smaller prior segments; `0` = `max(m / 4, 2)`. |
| `estimate_scale_from_prior` | bool | `true` | Estimate scale and scale variation from retained prior segments. |

## `[plotting]`

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `gene_composition_neighborhood` | int | `0` | Spatial neighbors for NCVs; auto = `max(n_genes / 10, m, 3)`. |
| `gene_composition_neigborhood` | int | — | Julia-era misspelling, still accepted; the correctly spelled key takes precedence. |
| `min_pixels_per_cell` | int | `15` | Accepted for compatibility; unused by the C++ pipeline. |
| `max_plot_size` | int | `3000` | Longer side of molecule images in HTML reports, in pixels; values below 1 use the default. |
| `max_z_slices` | int | `10` | Maximum z-layers for 3D polygons; larger stacks are binned. Must be ≥ 1. |
| `ncv_method` | string | `ri` | Compatibility values `ri`, `dense`, `sparse`; the C++ pipeline uses random indexing regardless. |

A run records config values and its invocation in
[`segmentation_params.dump.toml`](output_files.md#segmentation_paramsdumptoml)
(or `run_params.toml` for Parquet output). Keep your original config as well:
the dump does not serialize every setting.
