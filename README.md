# Baysor

**B**a**y**esian **S**egmentation **o**f Spatial T**r**anscriptomics Data

- [Abstract](#abstract)
  - [Method description](#method-description)
- [Installation](#installation)
  - [Install as a Julia package](#install-as-a-julia-package)
  - [Build CLI application from source](#build-cli-application-from-source)
  - [Docker](#docker)
- [Run](#run)
  - [Dataset preview](#dataset-preview)
  - [Full run](#full-run)
    - [Outputs](#outputs)
    - [Choice of parameters](#choice-of-parameters)
    - [Multi-threading](#multi-threading)

## Abstract

Spatial transcriptomics is an emerging stack of technologies, which adds spatial dimension to conventional single-cell RNA-sequencing. New protocols, based on in-situ sequencing or RNA fluorescent in-situ hybridization (RNA-FISH) register positions of single molecules in fixed tissue slices. Analysis of such data at the level of individual cells, however, requires accurate identification of cell boundaries. While many existing methods are able to approximate cell center positions using nuclei stains, current protocols do not report robust signal on the cell membranes, making cell segmentation a key challenge in downstream analysis and interpretation of the data. To address this challenge, we developed a tool for **B**a**y**esian **S**egmentation **o**f **S**patial T**r**anscriptomics Data (Baysor), which optimizes segmentations considering the likelihood of transcriptional composition, size and shape of the cell. The Bayesian approach can take into likely cell center positions from DAPI or other staining, however can also perform segmentation without any additional information. We show that Baysor segmentation can nearly double the number of the identified cells, while reducing contamination and decreasing cell doublet rate. We demonstrate Baysor performance on data acquired using different spatially-resolved protocols, with up to 10^5 of cells and 10^3 genes.

### Method description

The most relevant description is in the [poster](https://yadi.sk/i/gZ7U0brphSH91g) from SCG'19 conference and the [minipresentation](https://slides.com/vpetukhov/baysor_dsc19) from Danish Single-Cell Symposium'19

## Installation

### Install as a Julia package

To install Baysor as a Julia package run the following code:

```julia
using Pkg
Pkg.add([PackageSpec(name="UMAP", rev="master"), PackageSpec(url="https://github.com/hms-dbmi/Baysor.git")])
import Baysor
```

After that you can create Baysor executable by running

```bash
echo "#! /usr/bin/env julia\nimport Baysor: run_cli\nrun_cli()" >> baysor && chmod +x baysor
```

### Build CLI application from source

On linux you can build the command-line tool from source. To do that you need to clone this package and run the Makefile:

```bash
git clone https://github.com/hms-dbmi/Baysor.git
cd Baysor/bin
make
```

The command takes ~10 minutes to complete. It creates executable file `./Baysor/bin/baysor`, which can be used for spatial data segmentation.

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

**NOTE: The algorithm is still in the alpha-version, so it can be unstable now and is continuously improved. Please, report any problems and suggestions to Issues.**

### Dataset preview

As full run takes some time, it can be useful to run a quick preview to get meaning from the data and to get some guesses about parameters of the full run.

```bash
baysor preview [-x X_COL -y Y_COL --gene GENE_COL -c config.toml -o OUTPUT_PATH] MOLECULES_CSV
```

Here:

- MOLECULES_CSV must contain information about *x* and *y* positions and gene assignment for each transcript
- Parameters X_COL, Y_COL and GENE_COL must specify the corresponding column names. Default values are "x", "y" and "gene" correspondingly
- OUTPUT_PATH determines path to the output html file with the diagnostic plots

<!-- TODO: describe output format -->

### Full run

To run the algorithm on your data, use

```bash
baysor run [-s SCALE -x X_COL -y Y_COL --gene GENE_COL] -c config.toml MOLECULES_CSV [CENTERS_CSV]
```

Here:

- MOLECULES_CSV must contain information about *x* and *y* positions and gene assignment for each transcript
- CENTERS_CSV must have the same column name for *x* and *y* coordinates as MOLECULES_CSV
- Parameters X_COL, Y_COL and GENE_COL must specify the corresponding column names. Default values are "x", "y" and "gene" correspondingly.
- SCALE is a crutial parameter and must be approximately equal to the expected cell radius in the same units as "x" and "y". In general it's around 5-10 micrometers, and the preview run can be helpful to determine it for a specific dataset (by eyes, for now).

To see full list of command-line options run

```bash
baysor run --help
```

For more info see [examples](https://github.com/hms-dbmi/Baysor/tree/master/examples) (though probably are out of date). <!-- TODO: fix it -->

For the description of all config parameters, see [example_config.toml](https://github.com/hms-dbmi/Baysor/blob/master/configs/example_config.toml). <!-- TODO: chech that it's up to date -->

#### Outputs

- segmentation_borders.html: visualization of cell borders for the dataset colored by local gene expression composition (first part) and molecule clusters (second part)
- segmentation_cell_stats.csv: diagnostic info about cells. Parameters "n_transcripts", "density", "elongation", "area" and "avg_confidence" can be used to filter low-quality cells.
- segmentation_config.toml: copy of the config to improve reproducibility
- segmentation_params.dump: aggregated parameters from the config and CLI
- segmentation.csv: segmentation info per molecule:
  - confidence: probability of a molecule to be real (i.e. not noise)
  - cell: id of the assigned cell. Value "0" correspond to noise.
  - cluster: id of molecule cluster
  - assignment_confidence: confidence that the molecule is assigned to a correct cell
  - is_noise: shows whether molecule was assigned to noise (it's equal "true" if and only if "cell" == 0)

#### Choice of parameters

Most important parameters:

- `scale` is the most sensitive parameter, which specifies expected radius of a cell. It doesn't have to be precise, but wrong set-up can lead to over- or under-segmentation. This parameter is inferred automatically if cell centers are provided.
- `min-molecules-per-cell` is the number of molecules, required for a cell to be considered as real. It really depends on a protocol. For instance, for ISS it's fine to set it to 3, while for MERFISH it can require hundreds of molecule.

Some other sensitive parameters (normally, shouldn't be changed):

- `new-component-weight` is proportional to probability of generating new cell for a molecule, instead of assigning it to one of existing cells. More precisely, probability to assign a molecule to a particular cell linearly depends on number of molecules, already assigned to this cell. And this parameter is used as number of molecules for a cell, which is just generated for this new molecule. The algorithm is robust to small changes in this parameter. And normally values in range 0.1-0.9 should work fine. Smaller values would lead to slower convergence of the algorithm, while larger values force emergence of large number of small cells on each iteration, which can produce noise in the result. In general, default value should work well.
- `center-component-weight` has the same meaning as `new-component-weight`, but works for cell centers, estimated from DAPI staining (i.e. from CENTERS_CSV). In general, it should be larger then `new-component-weight` and it reflects your confidence in the estimated centers. For instance, if you want all centers, which have enough molecules around, to be used, you just need to set some huge value here (e.g. 100000).

Run parameters:

- `num-cells-init` expected number of cells in data. This parameter influence only convergence speed of the algorithm. It's better to set larger values than smaller ones.
- `iters` number of iterations for the algorithm. **At the moment, no convergence criteria is implemented, so it will work exactly `iters` iterations**. Thus, to small values would lead to non-convergence of the algorithm, while larger ones would just increase working time. Optimal values can be estimated by the convergence plots, produced among the results.
- `n-frames` determines maximum parallelization level. To allow multi-threaded processing of the data, all the space is splitted in `n_frames` parts, which contain approximately equal number of molecules. Results can be a bit inaccurate on the borders of the frames, so it's better to avoid very large values.

#### Multi-threading

Currently, running multi-threaded version is only available by splitting dataset to independent frames, which creates artifacts on the frame borders. If you're fine with that, please use the command

```bash
JULIA_NUM_THREADS=10; baysor run -n $JULIA_NUM_THREADS [-s SCALE -x X_COL -y Y_COL --gene GENE_COL] -c config.toml MOLECULES_CSV [CENTERS_CSV]
```