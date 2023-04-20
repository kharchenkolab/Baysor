[![Build Status](https://travis-ci.com/kharchenkolab/Baysor.svg?branch=master)](https://travis-ci.com/github/kharchenkolab/Baysor)
[![codecov](https://codecov.io/gh/kharchenkolab/Baysor/branch/master/graph/badge.svg?token=12AE1N3DY4)](undefined)

# Baysor

**Bay**esian **s**egmentation **o**f imaging-based spatial t**r**anscriptomics data

- [News (\[0.6.0\] — 2023-04-13)](#news-060--2023-04-13)
- [Overview](#overview)
- [Installation](#installation)
  - [Binary download](#binary-download)
  - [Install as a Julia package](#install-as-a-julia-package)
  - [Docker](#docker)
- [Run](#run)
  - [Dataset preview](#dataset-preview)
  - [Full run](#full-run)
    - [Normal run](#normal-run)
    - [Using a prior segmentation](#using-a-prior-segmentation)
      - [Segmenting stains](#segmenting-stains)
    - [Segmenting cells with pronounced intracellular structure](#segmenting-cells-with-pronounced-intracellular-structure)
    - [Outputs](#outputs)
    - [Choice of parameters](#choice-of-parameters)
  - [Extract NCVs (neighborhood composition vectors)](#extract-ncvs-neighborhood-composition-vectors)
  - [Multi-threading](#multi-threading)
- [Advanced configuration](#advanced-configuration)
- [Citation](#citation)

## News ([0.6.0] — 2023-04-13)

**A major update!**

- Simplified installation
- Fixed various issues: memory usage of NCVs, polygon output
- Finished loom support
- Updated CLI and configuration *(breaking changes in the config file structure!)*
- *And more*

*See the [changelog](CHANGELOG.md) for more detalis.*

## Overview

Baysor is a tool for performing cell segmentation on imaging-based spatial transcriptomics data. It optimizes segmentation considering the likelihood of transcriptional composition, size and shape of the cell. The approach can take into account nuclear or cytoplasm staining, however can also perform segmentation based on the detected molecules alone. The details of the method are described in the [paper](https://www.nature.com/articles/s41587-021-01044-w), or [pre-print](https://www.biorxiv.org/content/10.1101/2020.10.05.326777v1) (old version of the text). To reproduce analysis from the paper see [BaysorAnalysis](https://github.com/kharchenkolab/BaysorAnalysis) repo.

See the 16-min **[live-demo of Baysor](https://vimeo.com/558564804)** for an overview of the workflow!

**Do you have any question? Start [a discussion](https://github.com/kharchenkolab/Baysor/discussions)!**

## Installation

### Binary download

The easiest way to install Baysor on Linux or MacOS is to download a binary from the [release section](https://github.com/kharchenkolab/Baysor/releases) (see *Assets* for your platform). There, you can use *bin/baysor* executable.

### Install as a Julia package

If you need to install julia, [`juliaup`](https://github.com/julialang/juliaup#installation) is a recommended way. It's cross-platform and doesn't require admin privileges. *TL;DR: `curl -fsSL https://install.julialang.org | sh`* .

To install Baysor as a Julia package run the following command from your CLI *(it requires `gcc` or `clang` installed)*:

```bash
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.build()'
```

It will install all dependencies, compile the package and create an executable in `~/.julia/bin/baysor`. *This executable can be moved to any other place if you need it.*

### Docker

Alternatively, you can use Docker. It contains executable `baysor` to run Baysor from CLI, as well as IJulia installation to use Baysor with Jupyter.
The [repo](https://hub.docker.com/r/vpetukhov/baysor) also has images for older versions.

```bash
docker run -it --rm vpetukhov/baysor:latest
```

Build by hands:

```bash
docker pull julia:latest
cd Baysor/docker
docker build .
```

You can find more info about dockers at [Docker Cheat Sheet](https://github.com/wsargent/docker-cheat-sheet).

## Run

### Dataset preview

As full run takes some time, it can be useful to run a quick preview to get meaning from the data and to get some guesses about parameters of the full run.

```bash
baysor preview [-x X_COL -y Y_COL --gene GENE_COL -c config.toml -o OUTPUT_PATH] MOLECULES_CSV
```

Here:

- MOLECULES_CSV must contain information about *x* and *y* positions and gene assignment for each molecule
- Parameters X_COL, Y_COL and GENE_COL must specify the corresponding column names. Default values are "x", "y" and "gene" correspondingly
- OUTPUT_PATH determines path to the output html file with the diagnostic plots

### Full run

#### Normal run

To run the algorithm on your data, use

```bash
baysor run [-s SCALE -x X_COL -y Y_COL -z Z_COL --gene GENE_COL] -c config.toml MOLECULES_CSV [PRIOR_SEGMENTATION]
```

Here:

- MOLECULES_CSV must contain information about *x* and *y* positions and gene assignment for each molecule
- PRIOR_SEGMENTATION is optional molecule segmentation obtrained from other method (see [Using prior segmentation](#using-prior-segmentation))
- Parameters X_COL, Y_COL, Z_COL and GENE_COL must specify the corresponding column names. Default values are "x", "y", "z" and "gene" correspondingly.
- SCALE is a crutial parameter and must be approximately equal to the expected cell radius in the same units as "x", "y" and "z". In general, it's around 5-10 micrometers, and the preview run can be helpful to determine it for a specific dataset (by eyes, for now).

To see full list of command-line options run

```bash
baysor run --help
```

For more info see [examples](https://github.com/kharchenkolab/Baysor/tree/master/examples) (though probably are out-of-date). <!-- TODO: fix it -->

For the description of all config parameters, see [example_config.toml](https://github.com/kharchenkolab/Baysor/blob/master/configs/example_config.toml). <!-- TODO: check that it's up to date -->

#### Using a prior segmentation

In some cases, you may want to use another segmentation as a prior for Baysor. The most popular case is having a segmentation based on DAPI/poly-A stainings: such information helps to understand where nucleis are positioned, but it's often quite imprecise. To take this segmentation into account you can pass it as the second positional arguments to Baysor:

```bash
baysor run [ARGS] MOLECULES_CSV [PRIOR_SEGMENTATION]
```

Here, `PRIOR_SEGMENTATION` can be a path to a binary image with segmentation mask, an image with integer cell segmentation labels or a column name in the `MOLECULES_CSV` with integer cell assignment per molecule. In the later case, column name must have `:` prefix, e.g. for column `cell` you should use `baysor run [ARGS] molecules.csv :cell`. In case the image is too big to be stored in the tiff format, Baysor supports MATLAB '.mat' format: it should contain a single field with an integer matrix for either a binary mask or segmentation labels. When loading the segmentation, Baysor filters segments that have less than `min-molecules-per-segment` molecules. It can be set in the toml config, and the default value is `min-molecules-per-segment = min-molecules-per-cell / 4`. **Note:** only CSV column prior is currently supported for 3D segmentation.

To specify expected quality of the prior segmentation you may use `prior-segmentation-confidence` parameter. The value `0.0` makes the algorithm ignore the prior, while the value `1.0` restricts the algorithm from contradicting the prior. Prior segmentation is mainly needed for the cases where gene expression signal is not enough, e.g. with very sparse protocols (such as ISS or DARTFISH). Another potential usecase is high-quality data with visible sub-cellular structure. In these situations, setting `prior-segmentation-confidence > 0.7` is recommended. Otherwise, the default value `0.2` should work well.

##### Segmenting stains

If you have a non-segmented DAPI image, the simplest way to segment it would go through the following steps [ImageJ](https://fiji.sc/):

1. Open the image (File -> Open)
2. Go to Image -> Type and pick "8-bit"
3. Run Process -> Filters -> Gaussian Blur, using Sigma = 1.0. *The value can vary, depending on your DAPI, but 1.0 generally works fine.*
4. Run Image -> Adjust -> Auto Threshold, using Method = Default. *Different methods can give the best results for different cases. Often "Mean" also works well.*
5. Run Process -> Binary -> Watershed
6. Save the resulting image in the .tif

Another promising tool is [CellPose](https://github.com/mouseland/cellpose), however it may require some manual labeling to fine-tune the network.

#### Segmenting cells with pronounced intracellular structure

High-resolution protocols, such as MERFISH or seq-FISH, can capture intracellular structure. Most often, it would mean pronounced difference between nuclear and citoplasmic gene composition. By default, such differences would push Baysor to recognize compartments as different cells. However, if some compartment-specific genes are known, they may be used to mitigate the situation. These genes can be specified through `--config.segmentation.nuclei-genes` and `--config.segmentation.cyto-genes` options, *e.g.*:

```julia
baysor run -m 30 --n-clusters=1 -s 30 --scale-std=50% --config.segmentation.nuclei-genes=Neat1 --config.segmentation.cyto-genes=Apob,Net1,Slc5a1,Mptx2 --config.data.exclude-genes='Blank*' ./molecules.csv
```

Please, notice that it's highly recommended to set `--n-clusters=1`, so molecule clustering would not be affected by compartment differences.

**Note.** Currently, there is no automated way to determine such compartment-specific genes. So, the only way we can suggest is interactive explaration of data. In theory, it should be straightforward to infer such information from DAPI and poly-A stains, however it is not implemented yet. If you have particular need for such functionality, please submit an issue with description of your experimental set up.

#### Outputs

Segmentation results:
- ***segmentation_counts.loom*** or ***segmentation_counts.tsv*** (depends on `--count-matrix-format`): count matrix with segmented stats. In case of loom format, column attributes also contain the same info as ***segmentation_cell_stats.csv***.
- ***segmentation.csv***: segmentation info per molecule:
  - `confidence`: probability of a molecule to be real (i.e. not noise)
  - `cell`: id of the assigned cell. Value "" corresponds to noise.
  - `cluster`: id of molecule cluster
  - `assignment_confidence`: confidence that the molecule is assigned to a correct cell
  - `is_noise`: shows whether molecule was assigned to noise *(it's equal `true` if and only if `cell` == "")*
  - `ncv_color`: RGB code of the neighborhood composition coloring
- ***segmentation_cell_stats.csv***: diagnostic info about cells. The following parameters can be used to filter low-quality cells:
  - `area`: area of the convex hull around the cell molecules
  - `avg_confidence`: average confidence of the cell molecules
  - `density`: cell area divided by the number of molecules in cell
  - `elongation`: ration of the two eigenvalues of the cell covariance matrix
  - `n_transcripts`: number of molecules per cell
  - `avg_assignment_confidence`: average assignment confidence per cell. Cells with low `avg_assignment_confidence` have much higher chance of being an artifact.
  - `max_cluster_frac` *(only if `n-clusters > 1`)*: fraction of the molecules coming from the most popular cluster. Cells with low `max_cluster_frac` are often doublets.
  - `lifespan`: number of iterations the given component exists. The maximal `lifespan` is clipped proportionally to the total number of iterations. Components with short lifespan likely correspond to noise.

Visualization:
- ***segmentation_polygons.json***: polygons used for visualization in GeoJSON format. In case of 3D segmentation, it is an array of with GeoJSON polygons per z-plane, as well as "joint" polygons. *Shown only if `--save-polygons=geojson` is set*.
- ***segmentation_diagnostics.html***: visualization of the algoritm QC. *Shown only when `-p` is set.*
- ***segmentation_borders.html***: visualization of cell borders for the dataset colored by local gene expression composition (first part) and molecule clusters (second part). *Shown only when `-p` is set.*

Other:
- ***segmentation_params.dump.toml***: aggregated parameters from the config and CLI

#### Choice of parameters

Most important parameters:

- `scale` is the most sensitive parameter, which specifies expected radius of a cell. It doesn't have to be precise, but wrong set-up can lead to over- or under-segmentation. This parameter is inferred automatically if cell centers are provided.
- `min-molecules-per-cell` is the number of molecules, required for a cell to be considered as real. It really depends on a protocol. For instance, for ISS it's fine to set it to 3, while for MERFISH it can require hundreds of molecule.

Some other sensitive parameters (normally, shouldn't be changed):

- `new-component-weight` is proportional to probability of generating new cell for a molecule, instead of assigning it to one of existing cells. More precisely, probability to assign a molecule to a particular cell linearly depends on number of molecules, already assigned to this cell. And this parameter is used as number of molecules for a cell, which is just generated for this new molecule. The algorithm is robust to small changes in this parameter. And normally values in range 0.1-0.9 should work fine. Smaller values would lead to slower convergence of the algorithm, while larger values force emergence of large number of small cells on each iteration, which can produce noise in the result. In general, default value should work well.

Run parameters:

- `--config.segmentation.n-cells-init` expected number of cells in data. This parameter influence only convergence speed of the algorithm. It's better to set larger values than smaller ones.
- `--config.segmentation.iters` number of iterations for the algorithm. **At the moment, no convergence criteria is implemented, so it will work exactly `iters` iterations**. Thus, to small values would lead to non-convergence of the algorithm, while larger ones would just increase working time. Optimal values can be estimated by the convergence plots, produced among the results.


### Extract NCVs (neighborhood composition vectors)

To extract NCVs you need to run:

```bash
baysor segfree [-x X_COL -y Y_COL --gene GENE_COL -c config.toml -o OUTPUT_PATH] MOLECULES_CSV
```

The results will be stored in [loom format](https://linnarssonlab.org/loompy/format/index.html) with `/matrix` corresponding to the NCV matrix and `/col_attrs/ncv_color` showing the NCV color.

### Multi-threading

All running options support some basic multi-threading. To enable it, set `JULIA_NUM_THREADS` environment variable before running the script You can either do it globally by running `export JULIA_NUM_THREADS=13` or for an individual command:

```bash
JULIA_NUM_THREADS=13 baysor run -m 30 -s 30 ./molecules.csv
```

## Advanced configuration

The pipeline options are described in the CLI help. Run `baysor --help` or the corresponding command like `baysor run --help` for the list of main options.

However, there are additional parameters that can be spericied through the TOML config. See [example_config.toml](./configs/example_config.toml) for their description. All parameters from the config can also be passed through command line. For example, to set `exclude_genes` from the `data` section you need to pass `--config.data.exclude_genes='Blank*,MALAT1'` parameter.

## Preparing a release

```
export BAYSOR_VERSION=v0.6.0
docker build -t vpetukhov/baysor:latest -t "vpetukhov/baysor:$BAYSOR_VERSION" --build-arg CACHEBUST=$(date +%s) .
docker push vpetukhov/baysor

LazyModules_lazyload=false julia --project ./deps/build.jl app
zip -r "baysor-x86_x64-linux-${BAYSOR_VERSION}_build.zip" LICENSE README.md ./bin/baysor/*
```

## Citation

If you find Baysor useful for your publication, please cite:

```
Petukhov V, Xu RJ, Soldatov RA, Cadinu P, Khodosevich K, Moffitt JR & Kharchenko PV.
Cell segmentation in imaging-based spatial transcriptomics.
Nat Biotechnol (2021). https://doi.org/10.1038/s41587-021-01044-w
```
