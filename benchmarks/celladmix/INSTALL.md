# Installing the cellAdmix Python bindings for the benchmark suite

`install.sh` performs exactly the steps below; this file records the commands
as they were run on the development machine so a fresh machine can reproduce
the build.

## Pinned source

* repository: <https://github.com/kharchenkolab/cellAdmix-core>
* commit: `7d3fe7ae70c61d2b9e57469d38d9a88fcf6ac14d`
* clone location: `$BAYSOR_BENCH_DATA/cache/celladmix/src`

The clone is patched with three committed patches (applied by `install.sh`,
idempotently):

| patch | why |
|---|---|
| `patches/0001-bind-tabular-store.patch` | The stock Python bindings only expose `build_xenium_store` (`celladmix.dataset.CellAdmix.ensure_store` raises `NotImplementedError` for other formats), while the C++ core has a full tabular store builder (`build_tabular_input_store`). The patch adds `_core.build_tabular_store(...)` — ~50 lines mirroring the Xenium binding — so contract-format molecule tables work through the normal `CellAdmix` lifecycle. |
| `patches/0002-keepalive-tabular-parquet-reader.patch` | **Upstream bug**: `make_tabular_stream_source()` (src/input_store.cpp) creates the parquet reader as a local `std::unique_ptr` and returns a record-batch reader derived from it; the `FileReader` is destroyed on return while Arrow's async record-batch pipeline still references it. Every parquet tabular store build then segfaults in `parquet::ParquetFileReader::metadata()` on an Arrow IO thread (reproduced with a pure-C++ repro outside Python; CSV inputs and the Xenium path are unaffected). The patch keeps the reader alive alongside the record-batch reader. Worth reporting upstream. |
| `patches/0003-tests-include-algorithm.patch` | `tests/test_bridge.cpp` misses `#include <algorithm>`; gcc 15 rejects it. Only needed to build the C++ test suite. |

## System dependencies

Everything except OpenJPEG already exists in the C++ toolchain env
`.deps/env` (conda-forge): `libarrow`/`libparquet` 25.0.0 (with
`ArrowConfig.cmake`/`ParquetConfig.cmake`), `eigen`, `nlohmann_json`,
`libtiff`, `zlib`, `pkg-config`, `cmake` 4.4.3, `ninja`, conda `gcc/g++` 15.3.

OpenJPEG is installed into a dedicated local prefix so the shared env is not
modified:

```bash
export MAMBA_ROOT_PREFIX=/home/vpetukhov/Projects/Baysor/.deps/mamba
micromamba create -y -p /home/vpetukhov/Projects/Baysor/.deps/celladmix-build \
  -c conda-forge openjpeg
```

Python build requirements go into the shared bench env (also listed in
`benchmarks/environment.yml`):

```bash
/home/vpetukhov/Projects/Baysor/.deps/bench/bin/pip install scikit-build-core pybind11
```

## Clone, patch, build

```bash
git clone https://github.com/kharchenkolab/cellAdmix-core \
  "$BAYSOR_BENCH_DATA/cache/celladmix/src"
git -C "$BAYSOR_BENCH_DATA/cache/celladmix/src" checkout 7d3fe7ae70c61d2b9e57469d38d9a88fcf6ac14d
for p in patches/*.patch; do git -C "$BAYSOR_BENCH_DATA/cache/celladmix/src" apply "$p"; done

DEPS=/home/vpetukhov/Projects/Baysor/.deps
SRC=$BAYSOR_BENCH_DATA/cache/celladmix/src
export PATH=$DEPS/env/bin:$PATH \
       CC=$DEPS/env/bin/x86_64-conda-linux-gnu-cc \
       CXX=$DEPS/env/bin/x86_64-conda-linux-gnu-c++ \
       PKG_CONFIG_PATH=$DEPS/celladmix-build/lib/pkgconfig
CMAKE_ARGS="-DCMAKE_PREFIX_PATH=$DEPS/env;$DEPS/celladmix-build \
            -Dpybind11_DIR=$($DEPS/bench/bin/python -m pybind11 --cmakedir)" \
CMAKE_BUILD_PARALLEL_LEVEL=6 \
  $DEPS/bench/bin/pip install --no-build-isolation "$SRC/python"
```

`--no-build-isolation` keeps `scikit-build-core`/`pybind11` from the bench
env and lets CMake/Ninja/compilers resolve to `.deps/env` (the default
isolated build would download its own CMake/Ninja wheels, which the suite
rules disallow).

## Verification

```bash
$DEPS/bench/bin/python -c \
  "import celladmix as ca; from celladmix import _core; \
   print(ca.__version__, hasattr(_core, 'build_tabular_store'))"
# 0.0.1 True
```

The C++ test suite (optional; needs patch 0003):

```bash
cd "$SRC"
cmake --preset tests
cmake --build --preset tests -j6
./build/tests/celladmix_core_tests    # 46/46 pass at the pinned commit
```

## Why not the R bindings?

The R bindings do support tabular input (`cellAdmix(format = "tabular", ...)`)
and would not need patch 0001, but they would need a separate R micromamba
env plus `R CMD INSTALL` of the whole C++ core, and the whole benchmark suite
is Python-based. The Python path builds and works (the only blocker was the
missing binding plus the upstream parquet reader lifetime bug, both patched
above), so the R fallback described in the task was not needed.

## Installed artifact

* `celladmix-0.0.1` (compiled `_core` extension) in
  `.deps/bench/lib/python3.11/site-packages/celladmix/`
* runtime dependencies (numpy, pandas, pyarrow, scipy) were already present
  in the bench env
* wheel size ≈ 1.4 MB; build takes ≈ 3–4 minutes at `-j6`
