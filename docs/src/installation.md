# Installation

There are several ways to install the package.

```@contents
Pages = ["installation.md"]
```

## Linux binary download

The easiest way to install Baysor on Linux is to download a binary from the [release section](https://github.com/kharchenkolab/Baysor/releases) (see *Assets*). There, you can use *bin/baysor* executable. For other platforms, "Install as a Julia package" is a recommended way.
*If you know how to reliably compile binaries for MacOS or Windows, please, let me know in issues or over email!*

## Install as a Julia package

If you need to install julia, [`juliaup`](https://github.com/julialang/juliaup#installation) is a recommended way. It's cross-platform and doesn't require admin privileges. *TL;DR: `curl -fsSL https://install.julialang.org | sh`* .

**The current version of Baysor isn't compatible with Julia 1.11.**

To install Julia 1.10, use the following command:

```bash
juliaup add 1.10
juliaup default 1.10
```

To install Baysor as a Julia package run the following command from your CLI *(it requires `gcc` or `clang` installed)*:

```bash
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")); Pkg.build()'
```

It will install all dependencies, compile the package and create an executable in `~/.julia/bin/baysor`. *This executable can be moved to any other place if you need it.*

The same command can be used to update the package to the latest version.

### Installing other versions

To install development version, use `Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="develop"))`:

```bash
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="develop")); Pkg.build()'
```

Other versions may be similarly installed by passing the version to the `rev` parameter (e.g., `rev="v0.5.2"`).


## Docker

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
