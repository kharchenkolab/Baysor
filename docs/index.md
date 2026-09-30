# Baysor

**Bay**esian **s**egmentation **o**f imaging-based spatial t**r**anscriptomics data.

Baysor segments imaging-based spatial transcriptomics data using spatial
position, local gene composition, and optional prior segmentations. The
approach can take nuclear or cytoplasmic staining into account, but can also
segment based on the detected molecules alone. The method is described in the
[Nature Biotechnology paper](https://www.nature.com/articles/s41587-021-01044-w).

This site documents the **C++ line** of Baysor: a native C++17 implementation
of the Baysor segmentation algorithm, distributed as a single `baysor` binary
with three subcommands:

- [`baysor run`](run.md) — cell segmentation
- [`baysor preview`](preview.md) — quick dataset overview
- [`baysor segfree`](segfree.md) — segmentation-free neighborhood composition
  vectors (NCVs)

The documentation is versioned. The version selector in the header switches
between the docs of released versions; `latest` tracks the most recent release.
The old Julia (Baysor.jl v0.7.x) documentation is kept as an archived version —
see [Migrating from Baysor.jl](migrating.md).

## Quick start

Install a [release binary](installation.md#release-binaries) or
[build from source](installation.md#building-from-source), then:

```bash
baysor run -m 30 --scale 8 -o out molecules.csv
```

`-m/--min-molecules-per-cell` and either `--scale` or a
[prior segmentation](priors.md) are the only things you always need to decide
on. For Xenium data, start from [the Xenium workflow](xenium.md) instead.

## What is implemented

- CSV / Parquet molecule tables and Xenium `experiment.xenium` manifests as
  input
- prior segmentation from transcript columns (`:column_name`), TIFF masks, or
  boundary tables
- `legacy` (CSV/GeoJSON/Loom) and `parquet` (Parquet/GeoParquet/10x-HDF5)
  [output styles](outputs.md)
- 2D and 3D segmentation with per-layer polygon output
- HTML diagnostic reports (`run --plot`, `preview`)

## Where to go next

- [Installation](installation.md) — release binaries, source builds, Docker
- [Running Baysor](run.md) — the main `run` subcommand
- [Input data](inputs.md) and [Configuration](configuration.md) — formats and
  all options
- [Outputs](outputs.md) — what Baysor writes
- [Examples](examples.md) — runnable protocol-specific datasets
- [Development](development.md) — tests, coverage, benchmarks, releasing
- [Citation](citation.md)
