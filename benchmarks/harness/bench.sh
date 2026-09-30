#!/usr/bin/env bash
# One-command benchmark: run (a suite of) datasets, compare against a
# baseline and print the report (compare.py prints it; bench.sh does not
# repeat it).
#
#   bench.sh --baysor PATH [--baseline NAME] [--expect same|improved|identical]
#            [--preset regular|release|refactor|algorithm]
#            [--datasets quick|full|all|<ids/globs>] [--kind sim|real]
#            [--replicates N] [--threads N] [--run-id ID] [--timeout S]
#            [--data-root PATH] [--celltypes-from NAME]
#            [--create-baseline NAME]
#
# Presets:
#   regular    suite `regular` from datasets/suites.yaml: runs the 1-thread
#              bitwise step (expect identical, bare run-id) and the 6-thread
#              coverage step with the cellAdmix audit (expect same, run-id
#              <id>-noise). ~20 min wall. compare.py --suite compares both
#              groups; a `same` group that fails ONLY on the binary-sha256
#              provenance rows (rebuilt binary) is downgraded to a pass with
#              a note. --expect overrides the non-bitwise group only, so
#              `--expect improved` judges an algorithm change (the bitwise
#              group is then skipped) and `--expect same` is the default.
#   release    suite `release`: quick + full at 6 threads x 3 replicates with
#              the audit (expect same, run-id <id>-noise) plus the quick tier
#              at 1 thread x 1 replicate (expect identical, bare run-id).
#              ~8 h wall. Same --expect/--baseline override semantics.
#   refactor   legacy single step: threads 1, 1 replicate, expect identical
#              (exact refactor gate on --datasets, default quick).
#   algorithm  legacy single step: threads 6, 3 replicates, expect improved
#              (algorithm gate on --datasets, default quick).
#
# For the suite presets the per-step --datasets/--threads/--replicates/
# --timeout come from the manifest and must not be given on the command
# line; --baseline and --expect override every compared group.
#
# Environment overrides:
#   BAYSOR_BIN   path to the Baysor binary (alternative to --baysor)
#   BENCH_PY     Python interpreter (default: $REPO/.deps/bench/bin/python,
#                then $BAYSOR_BENCH_DATA/../.deps/bench/bin/python, else python3)
#
# Exit code: 0 = comparison passed, 1 = comparison failed, 2 = usage/setup error.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

BAYSOR="${BAYSOR_BIN:-}"
BASELINE=""
EXPECT="same"
EXPECT_GIVEN=0
PRESET=""
SUITE=""
DATASETS="quick"
KIND=""
THREADS=6
REPLICATES=3
RUN_ID=""
TIMEOUT=""
DATA_ROOT=""
CREATE_BASELINE=""
CELLTYPES_FROM=""
PY="${BENCH_PY:-}"

# Preset pre-scan: applied as defaults, explicit flags below still win.
# EXPL_CONFIG tracks flags that a suite manifest owns per step.
prev=""
EXPL_CONFIG=0
for arg in "$@"; do
  case "$prev" in
    --preset)     PRESET="$arg" ;;
    --expect)     EXPECT_GIVEN=1 ;;
    --datasets|--threads|--replicates|--timeout) EXPL_CONFIG=1 ;;
  esac
  prev="$arg"
done
case "$PRESET" in
  "") ;;
  regular)     SUITE="regular" ;;
  release)     SUITE="release" ;;
  refactor)    EXPECT="identical"; THREADS=1; REPLICATES=1 ;;
  algorithm)   EXPECT="improved";  THREADS=6; REPLICATES=3 ;;
  *) echo "unknown preset: $PRESET (expected regular|release|refactor|algorithm)" >&2; exit 2 ;;
esac

while [[ $# -gt 0 ]]; do
  case "$1" in
    --baysor)           BAYSOR="$2"; shift 2 ;;
    --baseline)         BASELINE="$2"; shift 2 ;;
    --expect)           EXPECT="$2"; shift 2 ;;
    --preset)           shift 2 ;;          # handled in the pre-scan
    --datasets)         DATASETS="$2"; shift 2 ;;
    --kind)             KIND="$2"; shift 2 ;;
    --replicates)       REPLICATES="$2"; shift 2 ;;
    --threads)          THREADS="$2"; shift 2 ;;
    --run-id)           RUN_ID="$2"; shift 2 ;;
    --timeout)          TIMEOUT="$2"; shift 2 ;;
    --data-root)        DATA_ROOT="$2"; shift 2 ;;
    --celltypes-from)   CELLTYPES_FROM="$2"; shift 2 ;;
    --create-baseline)  CREATE_BASELINE="$2"; shift 2 ;;
    -h|--help)          grep '^# ' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

if [[ -n "$SUITE" && $EXPL_CONFIG -eq 1 ]]; then
  echo "error: --datasets/--threads/--replicates/--timeout are defined per" >&2
  echo "       step by the suite manifest (benchmarks/datasets/suites.yaml)" >&2
  exit 2
fi
if [[ -n "$SUITE" && -n "$CREATE_BASELINE" ]]; then
  echo "error: --create-baseline needs a single run; with a suite run" >&2
  echo "       baseline.py create --run-id <id> yourself afterwards" >&2
  exit 2
fi

[[ -n "$BAYSOR" ]] || { echo "error: give --baysor PATH or set BAYSOR_BIN" >&2; exit 2; }
[[ -x "$BAYSOR" ]] || { echo "error: not executable: $BAYSOR" >&2; exit 2; }
if [[ -z "$PY" ]]; then
  if [[ -x "$REPO/.deps/bench/bin/python" ]]; then
    PY="$REPO/.deps/bench/bin/python"
  elif [[ -n "${BAYSOR_BENCH_DATA:-}" && -x "$BAYSOR_BENCH_DATA/../.deps/bench/bin/python" ]]; then
    PY="$(cd "$BAYSOR_BENCH_DATA/.." && pwd)/.deps/bench/bin/python"
  else
    PY=python3
  fi
fi
# Keep the default run id short (13 chars): at 1 thread Baysor output currently
# depends on the output-path length (>= 18-char run ids flip some datasets),
# see "Determinism" in harness/README.md. The suite flow keeps the bare id for
# the 1-thread (identical) group and suffixes the other groups.
[[ -z "$RUN_ID" ]] && RUN_ID="b$(date +%y%m%d%H%M%S)"

if [[ -n "$SUITE" ]]; then
  ARGS=(--baysor "$BAYSOR" --suite "$SUITE" --run-id "$RUN_ID")
  [[ -n "$KIND" ]]             && ARGS+=(--kind "$KIND")
  [[ -n "$DATA_ROOT" ]]        && ARGS+=(--data-root "$DATA_ROOT")
  [[ -n "$CELLTYPES_FROM" ]]   && ARGS+=(--celltypes-from "$CELLTYPES_FROM")

  echo "== run: suite $SUITE (run-id base $RUN_ID) =="
  "$PY" "$HERE/run.py" "${ARGS[@]}"
  RUN_RC=$?

  CARGS=(--run-id "$RUN_ID" --suite "$SUITE")
  [[ -n "$BASELINE" ]]       && CARGS+=(--baseline "$BASELINE")
  [[ $EXPECT_GIVEN -eq 1 ]]  && CARGS+=(--expect "$EXPECT")
  [[ -n "$DATA_ROOT" ]]      && CARGS+=(--data-root "$DATA_ROOT")
  echo "== compare: suite $SUITE (run-id base $RUN_ID) =="
  "$PY" "$HERE/compare.py" "${CARGS[@]}"
  CMP_RC=$?

  if [[ $RUN_RC -ne 0 && $CMP_RC -eq 0 ]]; then
    exit 1
  fi
  exit $CMP_RC
fi

ARGS=(--baysor "$BAYSOR" --datasets "$DATASETS" --run-id "$RUN_ID"
      --replicates "$REPLICATES" --threads "$THREADS")
[[ -n "$KIND" ]]             && ARGS+=(--kind "$KIND")
[[ -n "$TIMEOUT" ]]          && ARGS+=(--timeout "$TIMEOUT")
[[ -n "$DATA_ROOT" ]]        && ARGS+=(--data-root "$DATA_ROOT")
[[ -n "$CELLTYPES_FROM" ]]   && ARGS+=(--celltypes-from "$CELLTYPES_FROM")

echo "== run: $RUN_ID (datasets=$DATASETS replicates=$REPLICATES threads=$THREADS expect=$EXPECT) =="
"$PY" "$HERE/run.py" "${ARGS[@]}"
RUN_RC=$?

if [[ -n "$CREATE_BASELINE" ]]; then
  BARGS=(create --run-id "$RUN_ID" --name "$CREATE_BASELINE")
  [[ -n "$DATA_ROOT" ]] && BARGS+=(--data-root "$DATA_ROOT")
  # single-replicate / identical runs legitimately have no noise floor
  [[ "$REPLICATES" -lt 3 ]] && BARGS+=(--allow-incomplete)
  [[ "$EXPECT" == "identical" ]] && BARGS+=(--identical)
  "$PY" "$HERE/baseline.py" "${BARGS[@]}"
  exit $?
fi

if [[ -z "$BASELINE" ]]; then
  echo "error: give --baseline NAME or --create-baseline NAME" >&2
  exit 2
fi

CARGS=(--run-id "$RUN_ID" --baseline "$BASELINE" --expect "$EXPECT")
[[ -n "$DATA_ROOT" ]] && CARGS+=(--data-root "$DATA_ROOT")
echo "== compare: run $RUN_ID vs baseline $BASELINE (expect $EXPECT) =="
"$PY" "$HERE/compare.py" "${CARGS[@]}"
CMP_RC=$?

if [[ $RUN_RC -ne 0 && $CMP_RC -eq 0 ]]; then
  # runs failed but comparisons somehow passed: still fail
  exit 1
fi
exit $CMP_RC
