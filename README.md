# Baysor

**B**a**y**esian **S**egmentation **o**f **S**patial T**r**anscriptomics Data

## Installation

### Install as a Julia package

To install Baysor as a Julia package run the following code:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/hms-dbmi/Baysor.git"))
import Baysor
```

After that you can create Baysor executable by running

```bash
echo "#! /usr/bin/env julia\nimport Baysor: run_cli\nrun_cli()" >> baysor && chmod +x baysor
```

### Build CLI application from source

To build the command-line tool you need to clone this package and run the Makefile:

```bash
git clone https://github.com/hms-dbmi/Baysor.git
cd Baysor/bin
make
```

The command takes ~10 minutes to complete. It creates executable file `./Baysor/bin/segment_data.jl`, which can be used for spatial data segmentation.

### Docker

Alternatevely, you can use Dockerfile. It contains executable `segment_data` to run Baysor from CLI, as well as IJulia installation to use Baysor with Jupyter.

Pre-built [image](https://hub.docker.com/r/vpetukhov/baysor):

```bash
docker run -it --rm vpetukhov/baysor:latest
```

Buil by hands:

```bash
cd Baysor/docker
docker build .
```

You can find more info about dockers at [Docker Cheat Sheet](https://github.com/wsargent/docker-cheat-sheet).

## Run

**NOTE: The algorithm is still in the alpha-version, so it can be unstable now and will definitely be improved. Please, report all found problems and suggestions to Issues.**

To run the algorithm on your data, use

```bash
./Baysor/bin/segment_data.jl -i 500 [-s SCALE -x X_COL -y Y_COL --gene GENE_COL] -c config.toml MOLECULES_CSV [CENTERS_CSV]
```

Here:

- MOLECULES_CSV must contain information about *x* and *y* positions and gene assignment for each transcript
- CENTERS_CSV must have the same column name for *x* and *y* coordinates as MOLECULES_CSV
- Paremters X_COL, Y_COL and GENE_COL must specify corresponding column names. Default values are "x", "y" and "gene" correspondingly.

To see full list of command-line options run

```bash
./Baysor/bin/segment_data.jl --help
```

For more info see [examples](https://github.com/hms-dbmi/Baysor/tree/master/examples).  
For the description of all config parameters, see [example_config.toml](https://github.com/hms-dbmi/Baysor/blob/master/configs/example_config.toml).


### Choise of parameters

Most important parameters:

- `scale` is the most sensitive parameter, which specifies expected radius of a cell. It doesn't have to be precise, but wrong set-up can lead to over- or under-segmentation. This parameter is inferred automatically if cell senters are provided.
- `min-molecules-per-cell` is the number of molecules, required for a cell to be considered as real. It really depends on a protocol. For instance, for osm-FISH it's fine to set it to 5, while for MERFISH it can require hundreds of molecule.

Some other sensitive parameters (normally, shouldn't be changed):

- `new-component-weight` is proportional to probability of generating new cell for a molecule, instead of assigning it to one of existing cells. More precisely, probability to assign a molecule to a particular cell linearly depends on number of molecules, alreadu assigned to this cell. And this parameter is used as number of molecules for a cell, which is just generated for this new molecule. The algorithm is robust to small changes in this parameter. And normally values in range 0.1-0.9 should work fine. Smaller values would lead to slower convergence of the algorithm, while larger values force emergence of large number of small cells on each iteration, which can produce noise in the result. In general, default value should work well.
- `center-component-weight` has the same meaning as `new-component-weight`, but works for cell centers, estimated from DAPI staining (i.e. from CENTERS_CSV). In general, it should be larger then `new-component-weight` and it reflects your confidence in the estimated centers. For instance, if you want all centers, which have enough molecules around, to be used, you just need to set some huge value here (e.g. 100000).

Run paremeters:

- `num-cells-init` expected number of cells in data. This parameter influence only convergence speed of the algorithm. It's better to set larger values than smaller ones.
- `iters` number of iterations for the algorithm. **At the moment, no convergence criteria is implemented, so it will work exactly `iters` iterations**. Thus, to small values would lead to non-convergence of the algorithm, while larger ones would just increase working time. Optimal values can be estimated by the convergence plots, produced among the results.
- `n-frames` determines parallelization level. To allow milti-threaded processing of the data, all the space is splitted in `n_frames` parts, which contain approximately equal number of molecules. Results can be a bit inaccurate on the borders of the frames, so it's better to avoid very large values.
