#!/usr/bin/env bash
# Run the Baysor test suite and generate gcov coverage reports.
#
# Invoked by the CMake custom target 'coverage' (BAYSOR_WITH_COVERAGE=ON):
#   cmake --build <build-dir> --target coverage
#
# Usage: run_coverage.sh <build-dir> <source-dir> <ctest> <gcovr> <gcov>
#
# Outputs (under <build-dir>/coverage/):
#   summary.txt   - text report (also printed to stdout)
#   coverage.json - pretty-printed JSON summary
#   index.html    - HTML report with per-file details
set -euo pipefail

BUILD_DIR=$1
SOURCE_DIR=$2
CTEST=$3
GCOVR=$4
GCOV=$5

# 1. Drop stale counters so the report only reflects this test run.
find "$BUILD_DIR" -name '*.gcda' -delete

# 2. Run the tests (from the build dir so ctest picks up CTestTestfile.cmake).
"$CTEST" --test-dir "$BUILD_DIR" --output-on-failure

# 3. Generate all three reports in one gcovr run. Lines that are not code
#    (e.g. a lone '}' that GCC tags with exception-cleanup blocks at -O0) and
#    compiler-generated exception/unreachable branches are left out, so the
#    numbers reflect the source that tests can actually exercise. Run from
#    the source dir so the relative filters src/ and include/baysor/ match;
#    the build dir is passed explicitly as the search path for .gcda files.
#    include/third_party/ is excluded defensively (the filters already drop it).
mkdir -p "$BUILD_DIR/coverage"
cd "$SOURCE_DIR"
"$GCOVR" \
    --root "$SOURCE_DIR" \
    --gcov-executable "$GCOV" \
    --filter 'src/' \
    --filter 'include/baysor/' \
    --exclude 'include/third_party/' \
    --exclude-noncode-lines \
    --exclude-throw-branches \
    --exclude-unreachable-branches \
    --txt "$BUILD_DIR/coverage/summary.txt" \
    --json-summary "$BUILD_DIR/coverage/coverage.json" \
    --json-summary-pretty \
    --html-details "$BUILD_DIR/coverage/index.html" \
    "$BUILD_DIR"

# 4. Print the text summary.
cat "$BUILD_DIR/coverage/summary.txt"
