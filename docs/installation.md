# Installation

## Release binaries

**Linux x86-64:**

```bash
curl -fLO https://github.com/kharchenkolab/Baysor/releases/download/cpp-0.9.0/baysor-0.9.0-linux-x86_64.tar.gz
tar -xzf baysor-0.9.0-linux-x86_64.tar.gz
./baysor-0.9.0-linux-x86_64/bin/baysor run -m 30 -s 8 molecules.csv
```

`30` and `8` are example parameters; see [Cell segmentation](run.md) before
choosing them for your data.

Download other platforms and `SHA256SUMS` from the
[cpp-0.9.0 release](https://github.com/kharchenkolab/Baysor/releases/tag/cpp-0.9.0):

| Archive | Requirements |
| --- | --- |
| `baysor-0.9.0-linux-x86_64.tar.gz` | Linux x86-64, glibc 2.28 or newer <!-- GLIBC_FLOOR --> (e.g. Debian 10+, RHEL 8+) |
| `baysor-0.9.0-macos-arm64.tar.gz` | macOS 12 or newer, Apple Silicon |
| `baysor-0.9.0-windows-x86_64.zip` | 64-bit Windows 10 or newer |

Extract the archive and use `bin/baysor` (`bin/baysor.exe` on Windows)
inside the extracted directory. Keep the Windows DLLs beside the executable.
The binaries use generic CPU baselines; AVX2 is not required.

Optionally verify the Linux archive before extracting it:

```bash
curl -fL -o SHA256SUMS \
  https://github.com/kharchenkolab/Baysor/releases/download/cpp-0.9.0/SHA256SUMS
sha256sum --check --ignore-missing SHA256SUMS
```

## Docker

From the directory containing `molecules.csv`:

```bash
docker run --rm -v "$PWD:/data" ghcr.io/kharchenkolab/baysor:0.9.0 run -m 30 -s 8 molecules.csv
```

The release image is Linux x86-64 and uses `baysor` as its entrypoint: pass
`run`, `preview` or `segfree` directly after the image name. `/data` is its
working directory, so paths are relative to the mounted directory. The image
runs as UID/GID 1000; if your user has another ID, add
`--user "$(id -u):$(id -g)"` so the results are writable.

Version tags (such as `0.9.0`) pin a release; `latest` tracks the newest stable
release. The older Docker Hub images `vpetukhov/baysor` (`v0.4`–`v0.7.1`)
contain the Julia implementation.

## Building from source

Use this only if a release binary does not suit your platform or you need to
modify Baysor. On Linux / macOS, the simplest build uses a private conda-forge
environment bootstrapped by `configure.sh`; no root access is needed:

```bash
git clone -b cpp-0.9.0 https://github.com/kharchenkolab/Baysor.git
cd Baysor
./configure.sh --deps=conda --install
./install/bin/baysor --version
```

The default dependency directory is `.deps`; the binary is installed in
`install/bin`. Other dependency modes are `system` (packages already
installed), `vcpkg` (requires `VCPKG_ROOT`) and `auto` (default: system if the
build tools and Arrow/Parquet are available, otherwise conda).

??? note "Using system packages"

    A source build needs a C++17 compiler, CMake ≥ 3.20, a build tool, Eigen ≥
    3.3, spdlog, CGAL, Arrow/Parquet, HDF5, nlohmann_json, libtiff and zlib.
    CMake fetches the header-only UMAP dependencies automatically.

    On Ubuntu 24.04, first enable the Apache Arrow apt repository:

    ```bash
    sudo apt-get update
    sudo apt-get install -y --no-install-recommends ca-certificates lsb-release wget
    wget https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
    sudo apt-get install -y --no-install-recommends ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
    rm ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
    sudo apt-get update
    sudo apt-get install -y --no-install-recommends \
      build-essential cmake ninja-build pkg-config git \
      libeigen3-dev libspdlog-dev libcgal-dev libarrow-dev libparquet-dev \
      libhdf5-dev nlohmann-json3-dev libtiff-dev zlib1g-dev
    ./configure.sh --deps=system --install
    ```

    On macOS with Homebrew:

    ```bash
    brew install cmake ninja pkg-config eigen spdlog cgal apache-arrow hdf5 nlohmann-json libtiff
    ./configure.sh --deps=system --install
    ```

    For dependencies in a non-standard location, pass their prefix to CMake:

    ```bash
    ./configure.sh --deps=system -- -DCMAKE_PREFIX_PATH=/path/to/prefix
    ```

The repository also provides a source-build Dockerfile:
`docker build -t baysor .` from the repository root.
See `./configure.sh --help` and [Development](development.md) for build and
test options. If configuration fails, the error lists the missing dependency;
include the configure log when [opening an issue](https://github.com/kharchenkolab/Baysor/issues).

## Legacy Julia implementation

For Baysor.jl v0.7.1, follow the [migration page](migrating.md#installation)
and [archived installation guide](https://kharchenkolab.github.io/Baysor/0.7.1/installation/).
Do not use an unpinned Julia package install: the default branch is now C++.
