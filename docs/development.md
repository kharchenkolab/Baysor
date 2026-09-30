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

The regression/quality benchmark suite lives in its own repository,
[baysor-benchmarks](https://github.com/VPetukhov/baysor-benchmarks) (this
repo's former `benchmarks/` directory, history preserved): cropped real
datasets and simulated datasets with known ground truth, a runner with
metrics and baseline comparison, and a cellAdmix admixture audit. It checks
that a change keeps simulated-data metrics within the noise floor of a
stored baseline and real-data segmentations essentially unchanged
(`--expect identical|same`), or improves accuracy without regressions
(`--expect improved`); datasets and baselines live under a local data dir
(`.bench-data`, never in git). To run it against a Baysor build, clone that
repository, create its Python env from `environment.yml`, and point the
harness at your binary explicitly, e.g.
`harness/bench.sh --baysor /path/to/baysor --preset regular --run-id r1`
(≈20 min; `--preset release` before a release, `--dry-run` to resolve the
plan without running Baysor). See its README for the dataset contract,
suites and baselines.

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
  push/PR and deploys one version per GitHub release to the `gh-pages` branch
  (`mike deploy --push --update-aliases <version> latest` +
  `mike set-default --push latest`). The `latest` alias only moves to the
  newest stable release: pre-releases, `workflow_dispatch` redeploys of older
  tags, and backport patch releases are deployed without touching `latest`.
  Promoting a pre-release to a full release redeploys it and moves `latest`
  if it is the newest stable version.
- `docs/tools/check_cli_docs.py` fails if the docs mention a CLI option or
  config key that does not exist in the sources, and warns about options that
  are not documented. It parses `src/cli/main.cpp`, `src/utils/options.cpp`,
  and `configure.sh`, so it runs without building Baysor.
- `docs/tools/migrate_gh_pages.py` is the one-time migration that moved the
  archived Julia site from `dev/` to `0.7.1/` on `gh-pages` (see the
  `Migrating from Baysor.jl` page and the script header for usage).
