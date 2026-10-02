#!/usr/bin/env bash
# Smoke-test a Baysor binary from a release archive.
#
#   packaging/smoke_test.sh <path/to/bin/baysor> [expected-version]
#
# Runs `--version` (compared with expected-version if given), `--help`, and a
# full `baysor run` on a small synthetic dataset (a grid of round cells,
# generated here with awk), then checks the number of segmented cells.
#
# Needs only bash and POSIX tools (awk, mktemp), so it works in bare distro
# containers and in Git Bash on Windows runners. BAYSOR_SMOKE_WRAPPER is an
# optional command prefix, e.g. "qemu-x86_64 -cpu qemu64".
set -euo pipefail

die() { printf 'smoke_test: FAIL: %s\n' "$*" >&2; exit 1; }
say() { printf 'smoke_test: %s\n' "$*"; }

[[ $# -ge 1 && $# -le 2 ]] || die "usage: $0 <baysor-binary> [expected-version]"
BIN="$1"
EXPECTED_VERSION="${2:-}"
[[ -f "$BIN" ]] || die "binary not found: $BIN"
case "$BIN" in /*|[A-Za-z]:*) ;; *) BIN="$PWD/$BIN" ;; esac
read -r -a WRAPPER <<< "${BAYSOR_SMOKE_WRAPPER:-}"

baysor() { "${WRAPPER[@]+"${WRAPPER[@]}"}" "$BIN" "$@"; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/baysor-smoke.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

# --- --version / --help -------------------------------------------------------
version_out="$(baysor --version)" || die "--version exited with $?"
version_out="${version_out%$'\r'}"
say "--version: $version_out"
if [[ -n "$EXPECTED_VERSION" && "$version_out" != "baysor $EXPECTED_VERSION" ]]; then
    die "expected 'baysor $EXPECTED_VERSION', got '$version_out'"
fi

baysor --help > help.txt || die "--help exited with $?"
for sub in run preview segfree; do
    grep -q "$sub" help.txt || die "--help does not list the '$sub' subcommand"
done
say "--help: OK"

# --- baysor run ---------------------------------------------------------------
# 12 x 12 grid of round cells (radius 6, spacing 16) with ~120 molecules each;
# 4 cell types with 5 marker genes each over a 20-gene panel, plus 1% uniform
# background.
awk -v seed=42 'BEGIN {
    srand(seed); OFS = ","; pi = 3.141592653589793
    n_grid = 12; spacing = 16; radius = 6; n_genes = 20
    print "x", "y", "gene"
    for (i = 0; i < n_grid; i++) for (j = 0; j < n_grid; j++) {
        cx = (i + 0.5) * spacing; cy = (j + 0.5) * spacing
        type = (i * n_grid + j) % 4
        n = 90 + int(rand() * 60)
        for (k = 0; k < n; k++) {
            r = radius * sqrt(rand()); a = 2 * pi * rand()
            if (rand() < 0.6) g = type * 5 + int(rand() * 5)
            else g = int(rand() * n_genes)
            printf "%.3f,%.3f,g%d\n", cx + r * cos(a), cy + r * sin(a), g
        }
    }
    extent = n_grid * spacing
    for (k = 0; k < 170; k++)
        printf "%.3f,%.3f,g%d\n", rand() * extent, rand() * extent, int(rand() * n_genes)
}' > molecules.csv
run_args=(molecules.csv -x x -y y -g gene -s 6 -m 20)
min_cells=100; max_cells=190

say "running: baysor run ${run_args[*]} -o seg"
start=$(date +%s)
if ! baysor run "${run_args[@]}" -o seg > run.log 2>&1; then
    tail -n 30 run.log >&2
    die "baysor run failed"
fi
say "baysor run finished in $(( $(date +%s) - start )) s"

for f in segmentation.csv segmentation_cell_stats.csv segmentation_counts.loom; do
    [[ -s "seg/$f" ]] || die "missing or empty output seg/$f"
done
n_cells=$(awk 'NR > 1' seg/segmentation_cell_stats.csv | wc -l | tr -d ' ')
say "segmented cells: $n_cells (expected $min_cells..$max_cells)"
if (( n_cells < min_cells || n_cells > max_cells )); then
    die "unexpected number of cells"
fi
say "OK"
