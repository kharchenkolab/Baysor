# Dataset preview

Check the data before choosing parameters for a full segmentation:

```bash
baysor preview -m 30 -o preview.html molecules.csv
```

Open `preview.html` in a browser. It shows molecule-confidence diagnostics,
local gene-composition colors and gene structure, without assigning molecules
to cells. `-m` must be positive (CLI or config); no scale or prior is needed.

For a large dataset, preview a crop first:

```bash
baysor preview -m 30 --x-min 0 --x-max 2000 --y-min 0 --y-max 2000 \
  -o preview.html molecules.csv
```

The bounds are in the input coordinate units. Whole-slide reports can be
large and slow to open.

CSV, Parquet and Xenium manifests are accepted; see [Input data](inputs.md).
Use `-c config.toml` for column mappings and filters, or `-x`, `-y`, `-z`,
`-g` for column names. The [Xenium preset](configuration.md#protocol-presets)
already maps Xenium columns.

`-o` names an HTML file (default `preview.html`), not a directory.
`--threads` controls worker threads as for [run](run.md#threading).
See the [CLI reference](cli.md) or `baysor preview --help` for all options.
