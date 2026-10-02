# Baysor

**Bay**esian **s**egmentation **o**f imaging-based spatial t**r**anscriptomics data.
Baysor finds cells from molecule positions and gene composition, with or
without a prior segmentation. This is the native C++ release **cpp-0.9.0**.

## Quick start

**Linux x86-64** — download the release binary and run on your molecule table:

```bash
curl -fLO https://github.com/kharchenkolab/Baysor/releases/download/cpp-0.9.0/baysor-0.9.0-linux-x86_64.tar.gz
tar -xzf baysor-0.9.0-linux-x86_64.tar.gz
./baysor-0.9.0-linux-x86_64/bin/baysor run -m 30 -s 8 molecules.csv
```

**Docker** — from the directory containing `molecules.csv`:

```bash
docker run --rm -v "$PWD:/data" ghcr.io/kharchenkolab/baysor:0.9.0 run -m 30 -s 8 molecules.csv
```

The table needs `x`, `y` and `gene` columns (optional `z` for 3D).
`30` molecules and radius `8` are examples, not universal settings.

| Setting | What to choose |
| --- | --- |
| `-m` | Minimum molecules expected in a real cell; choose for your protocol. |
| `-s` / `--scale` | Approximate cell radius in coordinate units. Alternatively, pass a prior as the second input and set `--prior-segmentation-confidence` (default `0.2`). |
| `-c` | TOML config file; explicit CLI flags override it. |
| `-o` | Output directory (default `segmentation`). |
| `--threads` | Worker threads; physical CPU cores by default. |

## Documentation

- [Installation](https://kharchenkolab.github.io/Baysor/latest/installation/) —
  macOS / Windows binaries, requirements, Docker and source builds
- [Cell segmentation](https://kharchenkolab.github.io/Baysor/latest/run/) —
  choosing parameters, using a prior and inspecting results
- [Xenium workflow](https://kharchenkolab.github.io/Baysor/latest/xenium/) and
  [examples](https://kharchenkolab.github.io/Baysor/latest/examples/)
- [Performance](https://kharchenkolab.github.io/Baysor/latest/performance/profiling/) —
  run time, memory and segmentation examples

The [documentation](https://kharchenkolab.github.io/Baysor/) is versioned per
release. For the old Julia implementation, see the
[migration guide](https://kharchenkolab.github.io/Baysor/latest/migrating/) and
[archived v0.7.1 docs](https://kharchenkolab.github.io/Baysor/0.7.1/).

## Citation

If you use Baysor in a publication, please cite:

```text
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
