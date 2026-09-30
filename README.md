# Baysor

**Bay**esian **s**egmentation **o**f imaging-based spatial t**r**anscriptomics data.

Baysor segments imaging-based spatial transcriptomics data using spatial
position, local gene composition, and optional prior segmentations. This
repository contains the native C++ implementation (the `cpp` line), a single
`baysor` binary with `run`, `preview`, and `segfree` subcommands.

- **Documentation:** [kharchenkolab.github.io/Baysor](https://kharchenkolab.github.io/Baysor/)
  (versioned per release; includes a
  [migration guide](https://kharchenkolab.github.io/Baysor/latest/migrating/)
  from Baysor.jl v0.7.x and the [archived Julia docs](https://kharchenkolab.github.io/Baysor/0.7.1/))
- **Release binaries:** [GitHub Releases](https://github.com/kharchenkolab/Baysor/releases)
  (Linux x86-64, macOS arm64, Windows x86-64, with `SHA256SUMS`)

## Quick start

```bash
baysor run -m 30 --scale 8 -o out molecules.csv
```

See the documentation for [installation](https://kharchenkolab.github.io/Baysor/latest/installation/)
(binaries, source builds, Docker) and the
[run reference](https://kharchenkolab.github.io/Baysor/latest/run/).

## Citation

```
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
