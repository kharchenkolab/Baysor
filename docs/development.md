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

The repository carries a regression/quality benchmark suite under
[benchmarks/](https://github.com/kharchenkolab/Baysor/tree/HEAD/benchmarks) —
cropped real datasets and simulated datasets with known ground truth, with
committed baselines and a cellAdmix admixture audit. See
[benchmarks/README.md](https://github.com/kharchenkolab/Baysor/blob/HEAD/benchmarks/README.md)
for the layout, dataset contract, and the `bench.sh` runner.

## Releasing

Release binaries are built and published automatically for every GitHub
release. The release procedure is documented in
[RELEASING.md](https://github.com/kharchenkolab/Baysor/blob/HEAD/RELEASING.md).

## Documentation site

The docs site (this site) is built with MkDocs + Material and versioned with
mike:

```bash
python -m pip install -r docs/requirements.txt
mkdocs build --strict     # build into site/
mkdocs serve              # local preview
python docs/tools/check_cli_docs.py   # docs <-> CLI/config consistency check
```

- `.github/workflows/docs.yml` builds the docs strictly on every docs-related
  push/PR and deploys a new version to the `gh-pages` branch on every GitHub
  release (`mike deploy --update-aliases <version> latest` +
  `mike set-default latest`).
- `docs/tools/check_cli_docs.py` fails if the docs mention a CLI option or
  config key that does not exist in the sources, and warns about options that
  are not documented. It parses `src/cli/main.cpp`, `src/utils/options.cpp`,
  and `configure.sh`, so it runs without building Baysor.
- `docs/tools/migrate_gh_pages.py` is the one-time migration that moved the
  archived Julia site from `dev/` to `0.7.1/` on `gh-pages` (see the
  `Migrating from Baysor.jl` page and the script header for usage).
