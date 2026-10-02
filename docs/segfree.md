# Segmentation-free analysis

Analyze local gene composition without assigning molecules to cells:

```bash
baysor segfree -m 30 -k 100 -o ncvs.loom molecules.csv
```

Baysor computes one neighborhood composition vector (NCV) per molecule from
its nearest neighbors, then log-transforms it. `-k` / `--k-neighbors` sets the
neighborhood size: larger values emphasize broader spatial patterns. If
omitted, it is `max(n_genes / 10, min_molecules_per_cell, 3)` (integer division).

`-m` is required (CLI or config) for noise-estimation defaults even when `-k`
is set. No scale or prior is needed. CSV, Parquet and Xenium manifests are
accepted; use `-c config.toml` for column mappings and filters. See
[Input data](inputs.md) and [Configuration](configuration.md).

## Output

`-o` names a [Loom file](https://linnarssonlab.org/loompy/format/index.html)
(default `ncvs.loom`), not a directory:

- `/matrix` — genes × molecules, one log-transformed NCV per column
- `/col_attrs/Name` — molecule names `V1` … `VN`
- `/col_attrs/ncv_color` — local gene-composition color
- `/col_attrs/confidence` — molecule confidence from the noise model

These are overlapping neighborhoods, not segmented cells or a cell count
matrix. `--threads` works as for [run](run.md#threading). See the
[CLI reference](cli.md) or `baysor segfree --help` for all options.
