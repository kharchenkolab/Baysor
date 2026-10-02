# Baysor

Baysor finds cells in imaging-based spatial transcriptomics data using
molecule positions and gene composition, with or without a prior
segmentation. These docs cover the native C++ release **cpp-0.9.0**.

## Quick start

**Linux x86-64** — install and run on your molecule table:

```bash
curl -fLO https://github.com/kharchenkolab/Baysor/releases/download/cpp-0.9.0/baysor-0.9.0-linux-x86_64.tar.gz
tar -xzf baysor-0.9.0-linux-x86_64.tar.gz
./baysor-0.9.0-linux-x86_64/bin/baysor run -m 30 -s 8 molecules.csv
```

**Docker** — from the directory containing `molecules.csv`:

```bash
docker run --rm -v "$PWD:/data" ghcr.io/kharchenkolab/baysor:0.9.0 run -m 30 -s 8 molecules.csv
```

See [Installation](installation.md) for macOS / Windows binaries and system
requirements. The table needs `x`, `y` and `gene` columns; an optional `z`
column enables 3D segmentation. CSV and Parquet are supported.

| Setting | What to choose |
| --- | --- |
| `-m` | Minimum molecules expected in a real cell; choose for your protocol. |
| `-s` / `--scale` | Approximate cell radius in coordinate units. Or pass a [prior](priors.md) as the second input and set `--prior-segmentation-confidence` (default `0.2`). |
| `-c` | [TOML config](configuration.md); explicit CLI flags override it. |
| `-o` | Output directory (default `segmentation`). |
| `--threads` | Worker threads; physical CPU cores by default. |

`30` molecules and radius `8` are examples, not universal settings. For
Xenium data, start with the [Xenium workflow](xenium.md).

## Next steps

- [Cell segmentation](run.md) — choose parameters and inspect the results
- [Dataset preview](preview.md) — check the data before a full run
- [Segmentation-free analysis](segfree.md) — analyze local gene composition
  without assigning cells
- [Examples](examples.md) — workflows for Xenium, ISS, osm-FISH and STARmap
- [Outputs](outputs.md) — molecule assignments, count matrices and polygons
- [Performance](performance/benchmarks.md) — run time, memory and accuracy

The method is described in the [Nature Biotechnology paper](citation.md).
Questions? Start a [discussion](https://github.com/kharchenkolab/Baysor/discussions).

Use the version selector for other releases. For Baysor.jl, see
[Migrating from Julia](migrating.md) and the
[archived v0.7.1 docs](https://kharchenkolab.github.io/Baysor/0.7.1/).
