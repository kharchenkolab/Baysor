[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://kharchenkolab.github.io/Baysor/dev)

# Baysor

**Bay**esian **s**egmentation **o**f imaging-based spatial t**r**anscriptomics data

- [News (\[0.7.0\] — 2024-09-13)](#news-070--2024-09-13)
- [Overview](#overview)
- [Usage](#usage)
- [Installation](#installation)
  - [Binary download](#binary-download)
  - [Install as a Julia package](#install-as-a-julia-package)
- [Citation](#citation)

## News ([0.7.0] — 2024-09-13)

- Improved integration with 10x Xenium
- Better parallelism
- Optimized speed
- Improved NCV coloring
- Reworked polygon outputs

*See the [changelog](CHANGELOG.md) for more detalis.*

## Overview

Baysor is a tool for performing cell segmentation on imaging-based spatial transcriptomics data. It optimizes segmentation considering the likelihood of transcriptional composition, size and shape of the cell. The approach can take into account nuclear or cytoplasm staining, however, can also perform segmentation based on the detected molecules alone. The details of the method are described in the [paper](https://www.nature.com/articles/s41587-021-01044-w), or [pre-print](https://www.biorxiv.org/content/10.1101/2020.10.05.326777v1) (old version of the text). To reproduce the analysis from the paper see [BaysorAnalysis](https://github.com/kharchenkolab/BaysorAnalysis) repo.

## Usage

See [the documentation](https://kharchenkolab.github.io/Baysor/) for usage instructions.

## Installation

*For more details and alternative ways of installation see [the documentation](https://kharchenkolab.github.io/Baysor/)*

### Binary download

The easiest way to install Baysor on Linux is to download a binary from the [release section](https://github.com/kharchenkolab/Baysor/releases) (see *Assets*). There, you can use *bin/baysor* executable. For other platforms, "Install as a Julia package" is a recommended way.

### Install as a Julia package

[Install Julia](https://github.com/julialang/juliaup#installation):

```bash
curl -fsSL https://install.julialang.org | sh
```

Install Baysor:

```bash
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.build()'
```

## Citation

If you find Baysor useful for your publication, please cite:

```
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
