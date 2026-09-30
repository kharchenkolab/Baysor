#!/usr/bin/env bash
# Validate that the cellAdmix admixture audit discriminates segmentation quality
# on a real Xenium pancreas crop: vendor baseline vs deliberately degraded
# versions of the vendor segmentation, seed stochasticity, and runtime.
#
# Prerequisites (see INSTALL.md / install.sh):
#   * celladmix Python bindings installed into the bench environment
#   * the crop built by fetch_pancreas.py
#
# Environment overrides:
#   BENCH_PYTHON  python of the bench env   (default /home/vpetukhov/Projects/Baysor/.deps/bench/bin/python)
#   BAYSOR_BENCH_DATA                      (default /home/vpetukhov/Projects/Baysor/.bench-data)
#   THREADS       threads per audit run     (default 6)
#   SEED          baseline seed             (default 1)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAYSOR_MAIN=${BAYSOR_MAIN:-/home/vpetukhov/Projects/Baysor}
DATA=${BAYSOR_BENCH_DATA:-$BAYSOR_MAIN/.bench-data}
PY=${BENCH_PYTHON:-$BAYSOR_MAIN/.deps/bench/bin/python}
DS=$DATA/cache/celladmix/datasets/pancreas_crop_quick
WORK=${WORK:-$DATA/cache/celladmix/work/validation}
THREADS=${THREADS:-6}
SEED=${SEED:-1}
SRC=$DATA/cache/celladmix/src

if [[ ! -f "$DS/molecules.parquet" ]]; then
  echo "missing crop; run: python fetch_pancreas.py --crop-id pancreas_crop_quick" >&2
  exit 1
fi
mkdir -p "$WORK"

audit() { # audit <name> <extra audit.py args...>
  local name=$1; shift
  echo "=== audit: $name"
  "$PY" -u "$ROOT/audit.py" --molecules "$DS/molecules.parquet" \
    --out "$WORK/$name.json" --work-dir "$WORK/$name.work" \
    --celladmix-source "$SRC" --threads "$THREADS" "$@"
}

# --- 1. vendor baseline; typing = cellAdmix quick clustering (fixed seed) ---
audit vendor --cell-column cell_vendor --seed "$SEED" \
  --save-celltypes "$WORK/vendor_celltypes.parquet"

# --- 2/3. stochasticity: different seed re-draws the quick clustering; ---
#          same seed must reproduce the result exactly.
audit vendor_seed2 --cell-column cell_vendor --seed "$((SEED + 1))"
audit vendor_repeat --cell-column cell_vendor --seed "$SEED"

# --- 4/5. border reassignments at 10% / 30%, typed by transferring the ---
#          baseline labels (same cell ids => direct reuse).
for spec in "10 0.1" "30 0.3"; do
  pct=${spec%% *}; frac=${spec##* }
  echo "=== degrade: border$pct (${frac} of assigned molecules)"
  "$PY" "$ROOT/degrade.py" --molecules "$DS/molecules.parquet" --cell-column cell_vendor \
    --mode border --fraction "$frac" \
    --out "$WORK/border$pct.assignment.parquet" --report "$WORK/border$pct.degrade.json" \
    > "$WORK/border$pct.degrade.stdout"
  audit "border$pct" --assignment "$WORK/border$pct.assignment.parquet" \
    --celltypes "$WORK/vendor_celltypes.parquet" --seed "$SEED"
done

# --- 6. 2um dilation (takes molecules from neighbours and background) ---
echo "=== degrade: dilate 2um"
"$PY" "$ROOT/degrade.py" --molecules "$DS/molecules.parquet" --cell-column cell_vendor \
  --mode dilate --distance 2.0 \
  --out "$WORK/dilate2.assignment.parquet" --report "$WORK/dilate2.degrade.json" \
  > "$WORK/dilate2.degrade.stdout"
audit dilate2 --assignment "$WORK/dilate2.assignment.parquet" \
  --celltypes "$WORK/vendor_celltypes.parquet" --seed "$SEED"

# --- 7. comparability demo: the same degraded segmentation typed by ---
#        re-clustering instead of transferring the baseline labels.
audit border30_recluster --assignment "$WORK/border30.assignment.parquet" --seed "$SEED"

# --- 8. summary -----------------------------------------------------------
echo "=== summarize"
"$PY" "$ROOT/summarize.py" \
  --variant "vendor=$WORK/vendor.json" \
  --variant "vendor_seed2=$WORK/vendor_seed2.json" \
  --variant "vendor_repeat=$WORK/vendor_repeat.json" \
  --variant "border10=$WORK/border10.json" \
  --variant "border30=$WORK/border30.json" \
  --variant "dilate2=$WORK/dilate2.json" \
  --variant "border30_recluster=$WORK/border30_recluster.json" \
  --out "$WORK/summary.json" --markdown "$WORK/summary.md"

# Small results only: copy summary + audit JSONs into the data-dir results store.
RES="$DATA/results/celladmix"
mkdir -p "$RES"
cp "$WORK/summary.json" "$WORK/summary.md" "$RES/"
for name in vendor vendor_seed2 vendor_repeat border10 border30 dilate2 border30_recluster; do
  cp "$WORK/$name.json" "$RES/validation_$name.json"
done
echo "=== done; results in $RES"
