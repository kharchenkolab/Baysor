# Installation

There are three ways to get the `baysor` binary: download a release archive,
pull a published Docker image, or build from source (the repository also
ships a Dockerfile for a self-built image).

## Legacy Julia implementation

The recommended path is the C++ release binary below. If you still need the
last Julia implementation (v0.7.1), pin the package revision explicitly:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/kharchenkolab/Baysor.git", rev="v0.7.1"))
Pkg.build("Baysor")
```

The repository's default branch is now the C++ line, so omitting `rev="v0.7.1"`
will not install the Julia package. See the [archived Julia documentation](https://kharchenkolab.github.io/Baysor/0.7.1/).

## Release binaries

Every GitHub release publishes prebuilt archives on the
[releases page](https://github.com/kharchenkolab/Baysor/releases):

| Asset | Platform |
| --- | --- |
| `baysor-<version>-linux-x86_64.tar.gz` | Linux x86-64 |
| `baysor-<version>-macos-arm64.tar.gz` | macOS Apple Silicon (arm64) |
| `baysor-<version>-windows-x86_64.zip` | Windows x86-64 |
| `SHA256SUMS` | checksums for all archives |

Each archive contains a single top-level directory
`baysor-<version>-<platform>/` with:

```text
baysor-<version>-<platform>/
  bin/baysor        # bin/baysor.exe (+ bundled DLLs) on Windows
  LICENSE
  README.md
```

The binaries are built against generic CPU baselines (plain x86-64 / arm64),
so they run on any CPU of their architecture — no AVX2-class CPU is required.
They need no additional system libraries: the Linux binary only needs glibc
2.28 or newer <!-- GLIBC_FLOOR --> (e.g. RHEL/Alma/Rocky 8+, Debian 10+,
Ubuntu 18.10+), the macOS binary needs macOS 12 or newer on Apple silicon,
and the Windows archive carries the DLLs it needs next to `baysor.exe`
(64-bit Windows 10 or newer). `baysor --version` prints the release version.

Example on Linux:

```bash
tar xzf baysor-0.8.3-linux-x86_64.tar.gz
./baysor-0.8.3-linux-x86_64/bin/baysor --version
./baysor-0.8.3-linux-x86_64/bin/baysor run --help
```

Optionally verify the download against `SHA256SUMS`:

```bash
sha256sum --check --ignore-missing SHA256SUMS
```

## Docker

Every release publishes a Docker image of the portable release binary,
built from the same archive as the other platforms:

- Docker Hub: [`vpetukhov/baysor`](https://hub.docker.com/r/vpetukhov/baysor)
- GitHub Container Registry: `ghcr.io/kharchenkolab/baysor`

Images exist from the C++ release that introduced them (the first one is
tagged with its version like all later ones); older tags on Docker Hub
(`v0.4`–`v0.7.1`) are the Julia-era images. Tags are the release version
(e.g. `0.8.3`) plus `latest`, which always points at the newest stable
release.

```bash
docker pull vpetukhov/baysor:latest
docker run --rm vpetukhov/baysor:latest --version

# segment a dataset from the current directory, mounted at /data
# (the image's working directory):
docker run --rm -v "$PWD:/data" vpetukhov/baysor:latest \
  run -m 30 --scale 8 -o /data/out /data/molecules.csv
```

The same images are on the GitHub Container Registry (handy when Docker Hub
is unreachable):

```bash
docker pull ghcr.io/kharchenkolab/baysor:latest
```

The image contains only the release binary, its license and a small Debian
base — no compilers. It runs as the non-root user `baysor` (UID/GID 1000)
with `/data` as the working directory, so files written to a bind-mounted
data directory are owned by UID 1000 on the host; pass
`--user "$(id -u):$(id -g)"` to write them as yourself instead
(pre-create the host directory, e.g. `mkdir -p data`, so Docker does not
create it root-owned).

### Building the image yourself

The repository ships a Dockerfile that builds `baysor` on Ubuntu 24.04 and
installs it as the image entrypoint (for the prebuilt-binary image see
`packaging/docker/Dockerfile` and `RELEASING.md`):

```bash
git clone https://github.com/kharchenkolab/Baysor.git
cd Baysor
docker build -t baysor .
docker run --rm -v "$PWD/data:/data" baysor run -m 30 --scale 8 -o /data/out /data/molecules.csv
```

The image entrypoint is `/usr/local/bin/baysor`, so arguments passed to
`docker run` are forwarded to the CLI.

## Building from source

`./configure.sh` is the supported front-end to the CMake build. It can take
its C++ dependencies from three places (`--deps=<mode>`):

- `conda` — a private conda-forge environment bootstrapped with micromamba
  into `--deps-dir` (default: `.deps`). No root and no system packages
  needed; the most reproducible option.
- `system` — packages already installed with your platform package manager
  (apt, Homebrew, …).
- `vcpkg` — vcpkg manifest mode (requires `VCPKG_ROOT`).
- `auto` (default) — `system` if CMake, a build tool, and Arrow/Parquet are
  already available, otherwise `conda`.

Minimal build + install with conda dependencies:

```bash
./configure.sh --deps=conda --install
./install/bin/baysor --help
```

With system packages on Ubuntu 24.04 (the Apache Arrow apt source is needed
for `libarrow-dev` / `libparquet-dev`):

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  ca-certificates lsb-release wget

wget https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
sudo apt-get install -y --no-install-recommends ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
rm ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb

sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential cmake ninja-build pkg-config git \
  libeigen3-dev libspdlog-dev libcgal-dev \
  libarrow-dev libparquet-dev libhdf5-dev nlohmann-json3-dev libtiff-dev zlib1g-dev

./configure.sh --deps=system --install
```

With Homebrew on macOS:

```bash
brew install cmake ninja pkg-config eigen spdlog cgal apache-arrow hdf5 nlohmann-json libtiff
./configure.sh --deps=system --install
```

With vcpkg on Windows (Visual Studio 2022 or newer with the C++ workload):

```powershell
git clone https://github.com/microsoft/vcpkg "$env:USERPROFILE\vcpkg"
& "$env:USERPROFILE\vcpkg\bootstrap-vcpkg.bat" -disableMetrics
$env:VCPKG_ROOT = "$env:USERPROFILE\vcpkg"

./configure.sh --deps=vcpkg --install
```

The binary lands in `install/bin/baysor` (`install/bin/baysor.exe` on
Windows). `./configure.sh --help` lists all build options (build types, tests,
coverage, CUDA, custom prefixes and compilers); see
[Development](development.md) for the test and coverage workflows.

### Dependencies

Baysor needs a C++17 toolchain, CMake `>= 3.20`, Ninja, and the following
libraries. Only CMake, the C++ standard, and Eigen have explicit minimums;
the rest are intentionally not pinned so that system package managers,
Homebrew, and vcpkg can provide compatible versions.

| Dependency | Required / known-working version |
| --- | --- |
| CMake | `>= 3.20` |
| C++ compiler | C++17 compiler; GCC 9.4.0 and Visual Studio 2022 are known to work |
| Ninja | Recent Ninja; 1.10.0 is known to work |
| Eigen3 | `>= 3.3` (>= 3.4.90 enables Eigen's own threaded GEMM; older versions run dense products single-threaded) |
| Threads | A C++17 `std::thread` implementation (pthreads on Linux/macOS; Baysor runs its own thread pool) |
| spdlog | Not pinned; 1.5.0 is known to work |
| CGAL | Not pinned; 5.0.2 is known to work |
| Arrow / Parquet | Not pinned; 19.0.1 is known to work; Arrow must include compute, CSV, and Parquet support |
| HDF5 | Not pinned; 1.10.x is known to work |
| nlohmann_json | Not pinned; 3.7.3 is known to work |
| libtiff | Not pinned; 4.1.0 is known to work |
| zlib | Not pinned; 1.3.2 is known to work (PNG images of the HTML reports; usually already present as a dependency of HDF5 and libtiff) |

Several header-only UMAP dependencies are fetched automatically by CMake with
pinned source tags: `aarand` `v1.0.2`, `CppKmeans` `v3.1.1`, `subpar` `v0.3.1`,
`knncolle` `v2.3.0`, `CppIrlba` `v2.0.2`, `umappp` `v2.0.1`.

If a dependency is installed in a non-standard location, point CMake at it:

```bash
./configure.sh --deps=system -- -DCMAKE_PREFIX_PATH=/path/to/prefix
```

or at the package-specific config directory:

```bash
./configure.sh --deps=system -- -DArrow_DIR=/path/to/lib/cmake/arrow
```

## Troubleshooting

The CMake configure step checks each required dependency and prints the
package-manager command to install it when it is missing. For anything else,
please [open an issue](https://github.com/kharchenkolab/Baysor/issues) with the
configure log.
