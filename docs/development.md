# Development

## Building from source

`./configure.sh` is the front-end to the CMake build; `./configure.sh --help`
lists everything. The dependency modes (`--deps=conda|system|vcpkg|auto`) are
covered in [Installation](installation.md#building-from-source). The most
common build configurations:

```bash
# default developer build + tests + run them
./configure.sh --with-tests --test

# debug build with exported compile commands (clangd etc.)
./configure.sh --debug --compile-commands --build

# optimized install into ./install
./configure.sh --install

# optimise for this CPU
./configure.sh --native --build
```

Useful knobs:

- `--build-dir=DIR`, `--prefix=DIR`, `--build-type=Release|Debug|RelWithDebInfo|MinSizeRel`
- `--cc=PATH`, `--cxx=PATH`, `--generator=NAME`, `-DVAR=VALUE` (passed to
  CMake verbatim), and anything after `--` (forwarded to CMake verbatim)
- `--clean` reconfigures from scratch; `--build`, `--test`, `--install` run
  after configuring (in that order)
- environment variables `BAYSOR_DEPS`, `BAYSOR_DEPS_DIR`, `BAYSOR_BUILD_DIR`,
  `BAYSOR_PREFIX`, `BAYSOR_BUILD_TYPE`, `BAYSOR_GENERATOR`, `BAYSOR_JOBS`,
  `CC`, `CXX` provide defaults

After configuring, `build/baysor-env.sh` contains the toolchain environment
used for the build; source it to rebuild by hand with
`cmake --build build -j`.

The equivalent raw CMake flow for an end-user (no tests) build is
`cmake -P cmake/build_and_install.cmake`, which installs to `./install/bin`.

## Tests

The test suite is GTest-based and lives in `tests/`. Build and run it with:

```bash
./configure.sh --with-tests --test
```

`--test` runs `ctest` in the build directory. Please run the suite before
submitting changes; do not weaken tests to make them pass.

## Coverage

```bash
./configure.sh --coverage --test
```

`--coverage` instruments the build for gcov (implies `--with-tests` and adds
`gcovr` to conda deps). With `--test`, the `coverage` target writes HTML and
text reports to `<build-dir>/coverage/`.

## Benchmarks

The regression / accuracy suite lives in
[baysor-benchmarks](https://github.com/VPetukhov/baysor-benchmarks). It covers
real crops, simulations with ground truth and a cellAdmix audit. Clone that
repository, create its environment from `environment.yml` and point the
harness at your build:
`harness/bench.sh --baysor /path/to/baysor --preset regular --run-id r1`.
Use `--preset release` before a release; see its README for datasets,
baselines and comparison rules.

## Releasing

Release binaries are built and published automatically for every GitHub
release. The release procedure is documented in
[RELEASING.md](https://github.com/kharchenkolab/Baysor/blob/HEAD/RELEASING.md).

## Documentation site

The site uses MkDocs + Material and is versioned with mike. In an environment
with `docs/requirements.txt` installed, run from the repository root:

```bash
mkdocs build --strict
python3 docs/tools/check_cli_docs.py
```

See [docs/README.md](https://github.com/kharchenkolab/Baysor/blob/HEAD/docs/README.md)
for setup and editing notes. Publishing is handled by
[the docs workflow](https://github.com/kharchenkolab/Baysor/blob/HEAD/.github/workflows/docs.yml);
`latest` follows the newest stable release.
