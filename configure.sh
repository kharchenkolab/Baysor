#!/usr/bin/env bash
# configure.sh — one-stop configuration helper for the Baysor C++ build.
#
# Typical usage:
#   ./configure.sh                     # detect deps, configure ./build (Release)
#   ./configure.sh --build --install   # ...and compile, install to ./install
#   ./configure.sh --deps=conda --with-tests --build --test
#   ./configure.sh --coverage --debug --test --build-dir=build-cov  # coverage report
#
# Run `./configure.sh --help` for all options.

set -euo pipefail

SRC_DIR="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIGURE_ARGS=""
[[ $# -gt 0 ]] && CONFIGURE_ARGS="$(printf ' %q' "$@")"

# ----------------------------------------------------------------------------
# Defaults (every one of them can also be set through the environment)
# ----------------------------------------------------------------------------
DEPS_MODE="${BAYSOR_DEPS:-auto}"             # auto | system | conda | vcpkg
DEPS_DIR="${BAYSOR_DEPS_DIR:-$SRC_DIR/.deps}" # used by --deps=conda
BUILD_DIR="${BAYSOR_BUILD_DIR:-$SRC_DIR/build}"
PREFIX="${BAYSOR_PREFIX:-$SRC_DIR/install}"
BUILD_TYPE="${BAYSOR_BUILD_TYPE:-Release}"
GENERATOR="${BAYSOR_GENERATOR:-}"
JOBS="${BAYSOR_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
WITH_TESTS=OFF
WITH_CUDA=OFF
WITH_REPORTING=ON
COMPILE_COMMANDS=OFF
WITH_COVERAGE=OFF
NATIVE=OFF
CC_OVERRIDE="${CC:-}"
CXX_OVERRIDE="${CXX:-}"
CC_FROM_ENV="${CC:+1}"     # non-empty while CC/CXX come from the environment,
CXX_FROM_ENV="${CXX:+1}"   # not from --cc/--cxx
DO_BUILD=0
DO_INSTALL=0
DO_TEST=0
DO_CLEAN=0
EXTRA_CMAKE_ARGS=()

CONDA_CHANNEL="conda-forge"
CONDA_PACKAGES=(
    cmake ninja make pkg-config git
    c-compiler cxx-compiler
    eigen spdlog cgal-cpp
    libarrow libarrow-compute libparquet
    hdf5 nlohmann_json libtiff zlib
)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
if [[ -t 1 ]]; then
    C_B=$'\033[1m'; C_G=$'\033[32m'; C_Y=$'\033[33m'; C_R=$'\033[31m'; C_0=$'\033[0m'
else
    C_B=; C_G=; C_Y=; C_R=; C_0=
fi
info() { printf '%s==>%s %s\n' "$C_G$C_B" "$C_0" "$*"; }
warn() { printf '%swarning:%s %s\n' "$C_Y$C_B" "$C_0" "$*" >&2; }
die()  { printf '%serror:%s %s\n' "$C_R$C_B" "$C_0" "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }
abspath() { case "$1" in /*) printf '%s' "$1" ;; *) printf '%s' "$PWD/$1" ;; esac; }
cache_var() { sed -n "s/^$1:[A-Z]*=//p" "$BUILD_DIR/CMakeCache.txt" 2>/dev/null; }

usage() {
    cat <<EOF
Usage: ./configure.sh [options] [-- extra cmake args]

Dependencies:
  --deps=MODE           Where to take C++ dependencies from (default: auto)
                          system  packages already installed (apt, brew, ...)
                          conda   private conda-forge env in --deps-dir,
                                  bootstrapped with micromamba; no root needed
                          vcpkg   vcpkg manifest mode (requires VCPKG_ROOT)
                          auto    system if cmake (>= 3.20), ninja or make,
                                  and Arrow/Parquet are found on this machine,
                                  otherwise conda
  --deps-dir=DIR        Location of the conda env / micromamba (default: .deps)
  --update-deps         Re-solve and update the conda env even if it exists
                          (conda mode only; in other modes it is ignored with
                          a warning). Pinned conda specs must use the name=ver
                          form, e.g. sysroot_linux-64=2.28

Build configuration:
  --build-dir=DIR       CMake binary directory (default: build)
  --prefix=DIR          Install prefix (default: install)
  --build-type=TYPE     Release | Debug | RelWithDebInfo | MinSizeRel (default: Release)
  --debug               Shortcut for --build-type=Debug
  --with-tests          Build the GTest test suite (adds gtest to conda deps)
  --without-tests       Do not build the test suite (default)
  --with-cuda           Enable CUDA support (BAYSOR_WITH_CUDA=ON)
  --without-cuda        Disable CUDA support (default)
  --with-reporting      Enable reporting/visualisation (default)
  --without-reporting   Disable reporting/visualisation (BAYSOR_WITH_REPORTING=OFF)
  --coverage            Instrument for gcov coverage (implies --with-tests; adds
                        gcovr to conda deps). With --test, runs the 'coverage'
                        target, which writes <build-dir>/coverage/ reports
  --native              Optimise for this CPU (-march=native)
  --compile-commands    Export compile_commands.json (for clangd etc.)
  --cc=PATH, --cxx=PATH Override C / C++ compiler
  --generator=NAME      CMake generator (default: Ninja if available)
  -DVAR=VALUE           Passed to cmake verbatim (may be repeated)

Actions (--clean runs before configuring; the other actions run after
configuring, in the order build -> test -> install):
  --clean               Remove the build directory before configuring
  --build               Compile the baysor binary (and tests with --with-tests)
  --install             Install into --prefix (implies --build)
  --test                Run ctest (implies --with-tests and --build)
  -j N, --jobs=N        Parallel build jobs (default: $JOBS)

  -h, --help            Show this help

Environment variables BAYSOR_DEPS, BAYSOR_DEPS_DIR, BAYSOR_BUILD_DIR,
BAYSOR_PREFIX, BAYSOR_BUILD_TYPE, BAYSOR_GENERATOR, BAYSOR_JOBS, CC, CXX
provide defaults for the corresponding options.

After configuring, an env file is written to <build-dir>/baysor-env.sh; source
it to get the same toolchain on PATH, e.g. to rebuild by hand with
  cmake --build <build-dir> -j
EOF
}

# ----------------------------------------------------------------------------
# Argument parsing
# ----------------------------------------------------------------------------
UPDATE_DEPS=0
while [[ $# -gt 0 ]]; do
    arg="$1"; shift
    case "$arg" in
        --deps=*)            DEPS_MODE="${arg#*=}" ;;
        --deps-dir=*)        DEPS_DIR="${arg#*=}" ;;
        --update-deps)       UPDATE_DEPS=1 ;;
        --build-dir=*)       BUILD_DIR="${arg#*=}" ;;
        --prefix=*)          PREFIX="${arg#*=}" ;;
        --build-type=*)      BUILD_TYPE="${arg#*=}" ;;
        --debug)             BUILD_TYPE=Debug ;;
        --with-tests)        WITH_TESTS=ON ;;
        --without-tests)     WITH_TESTS=OFF ;;
        --with-cuda)         WITH_CUDA=ON ;;
        --without-cuda)      WITH_CUDA=OFF ;;
        --with-reporting)    WITH_REPORTING=ON ;;
        --without-reporting) WITH_REPORTING=OFF ;;
        --native)            NATIVE=ON ;;
        --compile-commands)  COMPILE_COMMANDS=ON ;;
        --coverage)          WITH_COVERAGE=ON; WITH_TESTS=ON ;;
        --cc=*)              CC_OVERRIDE="${arg#*=}"; CC_FROM_ENV= ;;
        --cxx=*)             CXX_OVERRIDE="${arg#*=}"; CXX_FROM_ENV= ;;
        --generator=*)       GENERATOR="${arg#*=}" ;;
        --clean)             DO_CLEAN=1 ;;
        --build)             DO_BUILD=1 ;;
        --install)           DO_BUILD=1; DO_INSTALL=1 ;;
        --test)              DO_BUILD=1; DO_TEST=1; WITH_TESTS=ON ;;
        -j)                  [[ $# -gt 0 ]] || die "-j needs a value"; JOBS="$1"; shift ;;
        -j*)                 JOBS="${arg#-j}" ;;
        --jobs=*)            JOBS="${arg#*=}" ;;
        -D*)                 EXTRA_CMAKE_ARGS+=("$arg") ;;
        --)                  EXTRA_CMAKE_ARGS+=("$@"); break ;;
        -h|--help)           usage; exit 0 ;;
        *)                   die "unknown option '$arg' (see --help)" ;;
    esac
done

case "$DEPS_MODE" in auto|system|conda|vcpkg) ;; *) die "--deps must be auto, system, conda or vcpkg" ;; esac
case "$BUILD_TYPE" in Release|Debug|RelWithDebInfo|MinSizeRel) ;; *) die "unsupported build type '$BUILD_TYPE'" ;; esac
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || die "jobs must be a positive integer, got '$JOBS'"
[[ -n "$DEPS_DIR"  ]] || die "--deps-dir must not be empty"
[[ -n "$BUILD_DIR" ]] || die "--build-dir must not be empty"
[[ -n "$PREFIX"    ]] || die "--prefix must not be empty"

DEPS_DIR="$(abspath "$DEPS_DIR")"
BUILD_DIR="$(abspath "$BUILD_DIR")"
PREFIX="$(abspath "$PREFIX")"
[[ "$WITH_TESTS" == ON ]] && CONDA_PACKAGES+=(gtest)
[[ "$WITH_COVERAGE" == ON ]] && CONDA_PACKAGES+=(gcovr)
# conda-forge's newest Linux sysroot (glibc 2.39) ships a crt1.o tagged as
# requiring x86-64-v3 (AVX2), so binaries refuse to start on older CPUs
# ("CPU ISA level is lower than required"). The 2.28 sysroot is conda-forge's
# default baseline and also keeps the binary portable to older distros.
case "$(uname -s)/$(uname -m)" in
    Linux/x86_64)  CONDA_PACKAGES+=("sysroot_linux-64=2.28") ;;
    Linux/aarch64) CONDA_PACKAGES+=("sysroot_linux-aarch64=2.28") ;;
esac

# ----------------------------------------------------------------------------
# Dependency detection
# ----------------------------------------------------------------------------
# Portable numeric cmake >= 3.20 check (no `sort -V`, which old sort lacks).
cmake_version_ok() {
    local v major minor
    v="$(cmake --version 2>/dev/null | head -n1 | awk '{print $3}')" || return 1
    major="${v%%.*}"; minor="${v#*.}"; minor="${minor%%.*}"
    [[ "$major" =~ ^[0-9]+$ && "$minor" =~ ^[0-9]+$ ]] || return 1
    (( 10#$major > 3 || (10#$major == 3 && 10#$minor >= 20) ))
}

system_deps_look_ok() {
    have cmake || return 1
    cmake_version_ok || return 1
    have ninja || have ninja-build || have make || return 1
    have "${CXX_OVERRIDE:-c++}" || return 1
    # Arrow/Parquet are the hardest dependencies to get; if they are present
    # the remaining ones (installed alongside in every guide) usually are too.
    if have pkg-config; then
        pkg-config --exists arrow parquet 2>/dev/null && return 0
    fi
    local d
    for d in /usr /usr/local /opt/homebrew; do
        [[ -e "$d/include/parquet/api/reader.h" ]] && return 0
    done
    return 1
}

if [[ "$DEPS_MODE" == auto ]]; then
    if [[ -x "$DEPS_DIR/env/bin/cmake" ]]; then
        DEPS_MODE=conda   # reuse a previously bootstrapped env
    elif system_deps_look_ok; then
        DEPS_MODE=system
    else
        DEPS_MODE=conda
    fi
    info "Dependency mode (auto-detected): $DEPS_MODE"
else
    info "Dependency mode: $DEPS_MODE"
fi

if [[ "$UPDATE_DEPS" == 1 && "$DEPS_MODE" != conda ]]; then
    warn "--update-deps only applies to conda mode; ignoring it in '$DEPS_MODE' mode"
fi

ENV_EXPORTS=()   # lines written to baysor-env.sh
CMAKE_ARGS=()

# ----------------------------------------------------------------------------
# conda: bootstrap micromamba + a private env with every build dependency
# ----------------------------------------------------------------------------
setup_conda() {
    local mm env_dir="$DEPS_DIR/env"
    mm="$(command -v micromamba 2>/dev/null || true)"
    if [[ -z "$mm" ]]; then
        mm="$DEPS_DIR/bin/micromamba"
        if [[ ! -x "$mm" ]]; then
            local os arch plat
            os="$(uname -s)"; arch="$(uname -m)"
            case "$os/$arch" in
                Linux/x86_64)            plat=linux-64 ;;
                Linux/aarch64)           plat=linux-aarch64 ;;
                Linux/ppc64le)           plat=linux-ppc64le ;;
                Darwin/x86_64)           plat=osx-64 ;;
                Darwin/arm64)            plat=osx-arm64 ;;
                *) die "no micromamba build for $os/$arch; install deps manually and use --deps=system" ;;
            esac
            info "Downloading micromamba ($plat) into $DEPS_DIR/bin"
            mkdir -p "$DEPS_DIR"
            have curl || die "curl is required to download micromamba"
            curl -fsSL "https://micro.mamba.pm/api/micromamba/$plat/latest" | tar -xj -C "$DEPS_DIR" bin/micromamba
        fi
    fi
    export MAMBA_ROOT_PREFIX="$DEPS_DIR/mamba"

    if [[ ! -x "$env_dir/bin/cmake" ]]; then
        info "Creating conda env $env_dir (this downloads ~1 GB once)"
        "$mm" create -y -q -p "$env_dir" -c "$CONDA_CHANNEL" --override-channels "${CONDA_PACKAGES[@]}"
    fi

    # Write the version pins (name=ver specs only) after env creation and
    # before any install/update, so `micromamba update --all` cannot silently
    # upgrade them (e.g. the sysroot, whose newest glibc ships an AVX2 crt1.o).
    # Pins are appended, so ones the user added by hand are kept.
    local pin pinned="$env_dir/conda-meta/pinned"
    for pin in "${CONDA_PACKAGES[@]}"; do
        [[ "$pin" == *=* ]] || continue
        grep -qxF "$pin" "$pinned" 2>/dev/null || printf '%s\n' "$pin" >> "$pinned"
    done

    if [[ "$UPDATE_DEPS" == 1 ]]; then
        info "Updating conda env $env_dir"
        "$mm" install -y -q -p "$env_dir" -c "$CONDA_CHANNEL" --override-channels "${CONDA_PACKAGES[@]}"
        "$mm" update -y -q -p "$env_dir" -c "$CONDA_CHANNEL" --override-channels --all
    else
        # Make sure optional packages (e.g. gtest) requested now are present.
        local missing=() p
        for p in "${CONDA_PACKAGES[@]}"; do
            local name="${p%%[=<>]*}" want="${p#*=}" hit
            hit="$(ls "$env_dir/conda-meta/$name-"[0-9]*.json 2>/dev/null | head -n1 || true)"
            if [[ -z "$hit" ]] || [[ "$p" == *=* && "$(basename "$hit")" != "$name-$want"* ]]; then
                missing+=("$p")
            fi
        done
        if [[ ${#missing[@]} -gt 0 ]]; then
            info "Adding to conda env: ${missing[*]}"
            "$mm" install -y -q -p "$env_dir" -c "$CONDA_CHANNEL" --override-channels "${missing[@]}"
        else
            info "Reusing conda env $env_dir"
        fi
    fi

    export PATH="$env_dir/bin:$PATH"
    export PKG_CONFIG_PATH="$env_dir/lib/pkgconfig:$env_dir/share/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
    ENV_EXPORTS+=("$(printf 'export PATH=%q:"$PATH"' "$env_dir/bin")")
    ENV_EXPORTS+=("$(printf 'export PKG_CONFIG_PATH=%q' "$PKG_CONFIG_PATH")")

    # Use the env's compilers so that libstdc++ matches the prebuilt libs.
    # Linux conda-forge names them *-conda-*-gcc/g++, macOS *-apple-darwin*-clang(++).
    TC_CC="$(ls "$env_dir"/bin/*-conda-*-gcc "$env_dir"/bin/*-apple-darwin*-clang 2>/dev/null | head -n1 || true)"
    TC_CXX="$(ls "$env_dir"/bin/*-conda-*-g++ "$env_dir"/bin/*-apple-darwin*-clang++ 2>/dev/null | head -n1 || true)"
    [[ -z "$CC_OVERRIDE"  && -n "$TC_CC"  ]] && CC_OVERRIDE="$TC_CC"
    [[ -z "$CXX_OVERRIDE" && -n "$TC_CXX" ]] && CXX_OVERRIDE="$TC_CXX"

    CMAKE_ARGS+=(
        "-DCMAKE_PREFIX_PATH=$env_dir"
        # Keep the link-time rpath so the installed binary finds the env's libs.
        "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON"
        "-DCMAKE_BUILD_RPATH=$env_dir/lib"
        "-DCMAKE_INSTALL_RPATH=$env_dir/lib"
    )
}

setup_vcpkg() {
    [[ -n "${VCPKG_ROOT:-}" ]] || die "--deps=vcpkg needs VCPKG_ROOT pointing to a bootstrapped vcpkg checkout"
    [[ -f "$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" ]] || die "VCPKG_ROOT=$VCPKG_ROOT does not look like a vcpkg checkout"
    CMAKE_ARGS+=(
        "-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
        "-DVCPKG_INSTALLED_DIR=$SRC_DIR/vcpkg_installed"
    )
    [[ "$WITH_TESTS" == ON ]] && CMAKE_ARGS+=("-DVCPKG_MANIFEST_FEATURES=tests")
    ENV_EXPORTS+=("$(printf 'export VCPKG_ROOT=%q' "$VCPKG_ROOT")")
}

setup_system() {
    have cmake || die "cmake not found. Install it (e.g. 'sudo apt-get install cmake ninja-build') or use --deps=conda"
    if have pkg-config && ! pkg-config --exists arrow parquet 2>/dev/null; then
        warn "Arrow/Parquet not visible to pkg-config; CMake may fail. See docs/installation.md or use --deps=conda"
    fi
}

case "$DEPS_MODE" in
    conda)  setup_conda ;;
    vcpkg)  setup_vcpkg ;;
    system) setup_system ;;
esac

# Resolve bare compiler names (e.g. CXX=g++) to full paths, so they compare
# equal to what CMake records in CMakeCache.txt and do not force a cache reset.
full_path() { case "$1" in ''|/*) printf '%s' "$1" ;; *) command -v "$1" 2>/dev/null || printf '%s' "$1" ;; esac; }
CC_OVERRIDE="$(full_path "$CC_OVERRIDE")"
CXX_OVERRIDE="$(full_path "$CXX_OVERRIDE")"

# In conda mode an exported CC/CXX would silently replace the env's toolchain
# and can produce GLIBCXX mismatches against the prebuilt conda libraries.
if [[ "$DEPS_MODE" == conda ]]; then
    if [[ -n "$CC_FROM_ENV" && -n "$TC_CC" && "$CC_OVERRIDE" != "$TC_CC" ]]; then
        warn "CC='$CC_OVERRIDE' comes from the environment and overrides the conda toolchain compiler '$TC_CC'"
    fi
    if [[ -n "$CXX_FROM_ENV" && -n "$TC_CXX" && "$CXX_OVERRIDE" != "$TC_CXX" ]]; then
        warn "CXX='$CXX_OVERRIDE' comes from the environment and overrides the conda toolchain compiler '$TC_CXX'"
    fi
fi

# ----------------------------------------------------------------------------
# Generic cmake arguments
# ----------------------------------------------------------------------------
have cmake || die "cmake not found"
CMAKE_VERSION="$(cmake --version | head -n1 | awk '{print $3}')"
cmake_version_ok || die "cmake >= 3.20 required, found $CMAKE_VERSION (use --deps=conda for a recent one)"

if [[ -z "$GENERATOR" ]]; then
    if have ninja || have ninja-build; then GENERATOR=Ninja; else GENERATOR="Unix Makefiles"; fi
fi

CMAKE_ARGS+=(
    "-G" "$GENERATOR"
    "-DCMAKE_BUILD_TYPE=$BUILD_TYPE"
    "-DCMAKE_INSTALL_PREFIX=$PREFIX"
    "-DBAYSOR_WITH_TESTS=$WITH_TESTS"
    "-DBAYSOR_WITH_CUDA=$WITH_CUDA"
    "-DBAYSOR_WITH_REPORTING=$WITH_REPORTING"
    "-DBAYSOR_EXPORT_COMPILE_COMMANDS=$COMPILE_COMMANDS"
    "-DBAYSOR_WITH_COVERAGE=$WITH_COVERAGE"
)
[[ -n "$CC_OVERRIDE"  ]] && CMAKE_ARGS+=("-DCMAKE_C_COMPILER=$CC_OVERRIDE")
[[ -n "$CXX_OVERRIDE" ]] && CMAKE_ARGS+=("-DCMAKE_CXX_COMPILER=$CXX_OVERRIDE")
# Always pass -DCMAKE_CXX_FLAGS (built from $CXXFLAGS) so that dropping
# --native takes effect on the next configure instead of sticking in the cache.
CXX_FLAGS="${CXXFLAGS:-}"
[[ "$NATIVE" == ON ]] && CXX_FLAGS+=" -march=native"
CMAKE_ARGS+=("-DCMAKE_CXX_FLAGS=$CXX_FLAGS")
CMAKE_ARGS+=("${EXTRA_CMAKE_ARGS[@]+"${EXTRA_CMAKE_ARGS[@]}"}")

# --clean runs before configuring and only removes a CMake build dir or an empty
# dir, never the source dir, one of its ancestors (e.g. /) or $HOME. Physical
# paths, so that symlinks and '..' cannot sneak past the checks.
if [[ "$DO_CLEAN" == 1 && -e "$BUILD_DIR" ]]; then
    real_dir() { CDPATH='' cd -- "$1" 2>/dev/null && pwd -P; }
    clean_dir="$(real_dir "$BUILD_DIR")" || die "--clean: '$BUILD_DIR' is not a directory"
    case "$(real_dir "$SRC_DIR")/" in
        "${clean_dir%/}/"*) die "--clean refuses to remove '$BUILD_DIR': it is or contains the source dir $SRC_DIR" ;;
    esac
    [[ "$clean_dir" != "$(real_dir "${HOME:-/}")" ]] || die "--clean refuses to remove '$BUILD_DIR' (\$HOME)"
    [[ -f "$BUILD_DIR/CMakeCache.txt" || -z "$(ls -A "$BUILD_DIR")" ]] \
        || die "--clean: '$BUILD_DIR' has no CMakeCache.txt and is not empty; refusing to delete it"
    info "Removing $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

# A cached generator/compiler/deps mode cannot be changed in place; start fresh.
if [[ "$DO_CLEAN" == 0 && -f "$BUILD_DIR/CMakeCache.txt" ]]; then
    cached_gen="$(cache_var CMAKE_GENERATOR)"
    cached_cxx="$(cache_var CMAKE_CXX_COMPILER)"
    cached_deps=""
    [[ -f "$BUILD_DIR/.baysor-deps" ]] && cached_deps="$(cat "$BUILD_DIR/.baysor-deps")"
    if [[ "$cached_gen" != "$GENERATOR" || "$cached_deps" != "$DEPS_MODE" ||
          ( -n "$CXX_OVERRIDE" && "$cached_cxx" != "$CXX_OVERRIDE" ) ]]; then
        warn "generator/compiler/deps mode changed since last configure; resetting CMake cache"
        rm -rf "$BUILD_DIR/CMakeCache.txt" "$BUILD_DIR/CMakeFiles" "$BUILD_DIR"/_deps/*-subbuild
    fi
fi

# ----------------------------------------------------------------------------
# Configure
# ----------------------------------------------------------------------------
info "Configuring in $BUILD_DIR"
printf '    cmake -S %q -B %q' "$SRC_DIR" "$BUILD_DIR"; printf ' %q' "${CMAKE_ARGS[@]}"; echo
cmake -S "$SRC_DIR" -B "$BUILD_DIR" "${CMAKE_ARGS[@]}"

# Stamp the deps mode this cache was configured for (used to detect changes).
printf '%s\n' "$DEPS_MODE" > "$BUILD_DIR/.baysor-deps"
{
    echo "# Generated by configure.sh on $(date -u '+%Y-%m-%d %H:%M UTC'). Source this file"
    echo "# to reproduce the configure-time environment."
    echo "# Invocation: ./configure.sh${CONFIGURE_ARGS}"
    [[ ${#ENV_EXPORTS[@]} -eq 0 ]] || printf '%s\n' "${ENV_EXPORTS[@]}"
    printf 'export BAYSOR_BUILD_DIR=%q\n' "$BUILD_DIR"
    printf 'export BAYSOR_PREFIX=%q\n' "$PREFIX"
} > "$BUILD_DIR/baysor-env.sh"

# ----------------------------------------------------------------------------
# Optional actions
# ----------------------------------------------------------------------------
TARGETS=(baysor)
[[ "$WITH_TESTS" == ON ]] && TARGETS+=(baysor_tests)

if [[ "$DO_BUILD" == 1 ]]; then
    info "Building ${TARGETS[*]} with $JOBS jobs"
    cmake --build "$BUILD_DIR" --config "$BUILD_TYPE" --parallel "$JOBS" --target "${TARGETS[@]}"
fi
if [[ "$DO_TEST" == 1 && "$WITH_COVERAGE" == ON ]] && have gcovr; then
    info "Running tests with coverage (reports in $BUILD_DIR/coverage)"
    cmake --build "$BUILD_DIR" --config "$BUILD_TYPE" --target coverage
elif [[ "$DO_TEST" == 1 ]]; then
    [[ "$WITH_COVERAGE" == ON ]] && warn "gcovr not found; running tests without a coverage report"
    info "Running tests"
    ctest --test-dir "$BUILD_DIR" --output-on-failure -j "$JOBS" -C "$BUILD_TYPE"
fi
if [[ "$DO_INSTALL" == 1 ]]; then
    info "Installing to $PREFIX"
    cmake --install "$BUILD_DIR" --config "$BUILD_TYPE"
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
summary_cxx="$(cache_var CMAKE_CXX_COMPILER)"
cat <<EOF

${C_B}Baysor configured${C_0}
  deps mode   : $DEPS_MODE
  build dir   : $BUILD_DIR
  build type  : $BUILD_TYPE
  generator   : $GENERATOR
  compiler    : ${summary_cxx:-(cmake default)}
  tests       : $WITH_TESTS   coverage: $WITH_COVERAGE   cuda: $WITH_CUDA   reporting: $WITH_REPORTING
  prefix      : $PREFIX

Next steps:
EOF
printf '  source %q\n' "$BUILD_DIR/baysor-env.sh"
printf '  cmake --build %q -j%s        # compile\n' "$BUILD_DIR" "$JOBS"
printf '  cmake --install %q              # install to %s/bin/baysor\n' "$BUILD_DIR" "$PREFIX"
[[ -x "$BUILD_DIR/baysor" ]] && printf '  %q --help           # binary is already built\n' "$BUILD_DIR/baysor"
exit 0
