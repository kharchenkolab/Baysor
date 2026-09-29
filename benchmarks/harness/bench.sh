#!/usr/bin/env bash
# One-command benchmark: run the quick tier, compare against a baseline and
# print the report (compare.py prints it; bench.sh does not repeat it).
#
#   bench.sh --baysor PATH --baseline NAME [--expect same|improved|identical]
#            [--preset refactor|algorithm]
#            [--datasets quick|full|all|<ids/globs>] [--kind sim|real]
#            [--replicates N] [--threads N] [--run-id ID] [--timeout S]
#            [--data-root PATH] [--celltypes-from NAME]
#            [--create-baseline NAME]
#
# Presets (set defaults; any explicit flag below overrides them):
#   refactor    threads 1, 1 replicate,  expect identical  (exact refactor gate)
#   algorithm   threads 6, 3 replicates, expect improved   (algorithm gate)
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
PRESET=""
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
prev=""
for arg in "$@"; do
  [[ "$prev" == "--preset" ]] && PRESET="$arg"
  prev="$arg"
done
case "$PRESET" in
  "") ;;
  refactor)   EXPECT="identical"; THREADS=1; REPLICATES=1 ;;
  algorithm)  EXPECT="improved";  THREADS=6; REPLICATES=3 ;;
  *) echo "unknown preset: $PRESET (expected refactor|algorithm)" >&2; exit 2 ;;
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
[[ -z "$RUN_ID" ]] && RUN_ID="bench-$(date +%Y%m%d-%H%M%S)"

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
