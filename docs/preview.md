# Dataset preview (`baysor preview`)

`preview` generates a self-contained HTML overview of a dataset without
running the full segmentation. It is useful for sanity checks and for
estimating parameters before a `run`:

```bash
baysor preview [OPTIONS] coordinates
```

## Typical use

```bash
baysor preview -c configs/xenium.toml -o preview.html data/transcripts.parquet
```

For Xenium-style columns without a config:

```bash
baysor preview -m 30 --qv-column qv -g feature_name -x x_location -y y_location \
  -o preview.html data/transcripts.parquet
```

## What it computes

- loads and filters the molecules
- estimates a molecule-confidence / noise model
- computes neighborhood-composition colors
- estimates a gene-structure embedding
- writes one HTML report with dataset diagnostics

The HTML output can be large for whole-slide datasets; run on a crop (see the
coordinate bounds below) first if in doubt.

## Options

| Option | Default | Description |
| --- | --- | --- |
| `coordinates` | — | required. CSV/Parquet molecule table, or a Xenium `experiment.xenium` manifest |
| `-c, --config` | — | TOML file with configuration |
| `-x, --x-column` | `x` | Name of the x column |
| `-y, --y-column` | `y` | Name of the y column |
| `-z, --z-column` | `z` | Name of the z column |
| `-g, --gene-column` | `gene` | Name of the gene column |
| `--qv-column` | `qv` | Name of the quality-value column used by `--min-qv` |
| `-m, --min-molecules-per-cell` | — | Minimal number of molecules for a cell to be considered real. Required (CLI or config) |
| `--min-qv` | `-1` | Drop molecules with quality value below this threshold |
| `--x-min`, `--x-max` | ±∞ | Keep only molecules within this x range |
| `--y-min`, `--y-max` | ±∞ | Keep only molecules within this y range |
| `--z-min`, `--z-max` | ±∞ | Keep only molecules within this z range |
| `-o, --output` | `preview.html` | Output HTML file |
| `--force-2d` | off | Ignore the z column in the data |
| `-t, --threads` | auto | Number of worker threads; auto = `BAYSOR_NUM_THREADS`, then `OMP_NUM_THREADS`, then CPU cores |

## Notes

- `preview` accepts a transcript table directly; for Xenium datasets it also
  accepts `experiment.xenium` and resolves the underlying transcript table
  automatically.
- the coordinate bounds make it cheap to preview a crop of a large dataset.
