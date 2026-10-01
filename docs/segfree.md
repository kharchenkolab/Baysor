# Segmentation-free analysis (`baysor segfree`)

`segfree` extracts neighborhood composition vectors (NCVs) without running
the cell segmentation algorithm — one vector per molecule, summarizing the
local transcriptional composition. Many analyses don't require segmentation
and can run on these local neighborhoods instead:

```bash
baysor segfree [OPTIONS] coordinates
```

## Typical use

```bash
baysor segfree -c configs/xenium.toml -k 100 -o ncvs.loom data/transcripts.parquet
```

## What it computes

- loads and filters the molecules
- builds per-molecule neighborhood-composition counts with `k` nearest
  neighbors and log-transforms them
- estimates per-molecule confidences (noise model)
- computes NCV colors
- writes a [loom](https://linnarssonlab.org/loompy/format/index.html) file
  with one NCV per molecule

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
| `-m, --min-molecules-per-cell` | — | Minimal number of molecules for a cell to be considered real. Required (CLI or config); used to derive `k` |
| `--min-qv` | `-1` | Drop molecules with quality value below this threshold |
| `--x-min`, `--x-max` | ±∞ | Keep only molecules within this x range |
| `--y-min`, `--y-max` | ±∞ | Keep only molecules within this y range |
| `--z-min`, `--z-max` | ±∞ | Keep only molecules within this z range |
| `-k, --k-neighbors` | auto | Number of neighbors per NCV. Auto = `max(n_genes / 10, min_molecules_per_cell, 3)` |
| `-o, --output` | `ncvs.loom` | Output Loom file |
| `--force-2d` | off | Ignore the z column in the data |
| `-t, --threads` | auto | Number of worker threads; auto = `BAYSOR_NUM_THREADS`, then `OMP_NUM_THREADS`, then physical CPU cores |

## Output

The Loom file contains:

- `/matrix` — the NCV matrix, one column per molecule
- `/col_attrs/ncv_color` — per-molecule NCV color
- `/col_attrs/confidence` — per-molecule confidence from the noise model

## Notes

- `segfree` accepts a transcript table directly; for Xenium datasets it also
  accepts `experiment.xenium` and resolves the underlying transcript table
  automatically.
- the output is segmentation-free: molecules are not grouped into cells, and
  every molecule is its own "cell" (`V1..VN`).
