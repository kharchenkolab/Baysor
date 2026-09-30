#!/usr/bin/env bash
# Step 4 of the cellAdmix validation: audit a Baysor segmentation of the same
# crop and compare it with the vendor segmentation.
#
# Baysor parameters come from the dataset's meta.json (scale_um, scale_std,
# min_molecules_per_cell, prior=none, extra_args). The binary defaults to all
# cores; OMP_NUM_THREADS caps it per the suite rules.
#
# Environment overrides: BAYSOR_BIN, BENCH_PYTHON, BAYSOR_BENCH_DATA, THREADS.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAYSOR_MAIN=${BAYSOR_MAIN:-/home/vpetukhov/Projects/Baysor}
DATA=${BAYSOR_BENCH_DATA:-$BAYSOR_MAIN/.bench-data}
PY=${BENCH_PYTHON:-$BAYSOR_MAIN/.deps/bench/bin/python}
BAYSOR_BIN=${BAYSOR_BIN:-/home/vpetukhov/Projects/Baysor/.bench-data/binaries/baysor-bugfixes-35e8a7e}
DS=$DATA/cache/celladmix/datasets/pancreas_crop_quick
WORK=${WORK:-$DATA/cache/celladmix/work/validation}
THREADS=${THREADS:-6}
SEED=${SEED:-1}
SRC=$DATA/cache/celladmix/src

mkdir -p "$WORK"

# --- 1. Baysor segmentation (skip when present) ---------------------------
if [[ ! -f "$WORK/baysor_seg/molecules.parquet" ]]; then
  echo "=== baysor run"
  OMP_NUM_THREADS=$THREADS "$BAYSOR_BIN" run "$DS/molecules.parquet" \
    -x x -y y -g gene \
    -s "$(python3 -c "import json;print(json.load(open('$DS/meta.json'))['baysor']['scale_um'])")" \
    --scale-std "$(python3 -c "import json;print(json.load(open('$DS/meta.json'))['baysor']['scale_std'])")" \
    -m "$(python3 -c "import json;print(json.load(open('$DS/meta.json'))['baysor']['min_molecules_per_cell'])")" \
    --force-2d \
    -o "$WORK/baysor_seg" --output-style parquet
fi

# --- 2. row-align the per-molecule labels to the dataset ------------------
echo "=== align"
"$PY" "$ROOT/align.py" --molecules "$DS/molecules.parquet" \
  --segmentation "$WORK/baysor_seg/molecules.parquet" --cell-column cell \
  --out "$WORK/baysor.assignment.parquet" --report "$WORK/baysor.align.json"

# --- 3. transfer the baseline cell types onto Baysor cells by molecule overlap
echo "=== transfer typing"
"$PY" "$ROOT/transfer.py" --molecules "$DS/molecules.parquet" \
  --baseline-cell-column cell_vendor \
  --target-assignment "$WORK/baysor.assignment.parquet" \
  --baseline-celltypes "$WORK/vendor_celltypes.parquet" \
  --out "$WORK/baysor_celltypes.parquet" --report "$WORK/baysor.transfer.json"

# --- 4. audit --------------------------------------------------------------
echo "=== audit: baysor"
"$PY" -u "$ROOT/audit.py" --molecules "$DS/molecules.parquet" \
  --assignment "$WORK/baysor.assignment.parquet" \
  --celltypes "$WORK/baysor_celltypes.parquet" \
  --out "$WORK/baysor.json" --work-dir "$WORK/baysor.work" \
  --celladmix-source "$SRC" --threads "$THREADS" --seed "$SEED"

# --- 5. summary including the baysor variant ------------------------------
"$PY" "$ROOT/summarize.py" \
  --variant "vendor=$WORK/vendor.json" \
  --variant "vendor_seed2=$WORK/vendor_seed2.json" \
  --variant "vendor_repeat=$WORK/vendor_repeat.json" \
  --variant "border10=$WORK/border10.json" \
  --variant "border30=$WORK/border30.json" \
  --variant "dilate2=$WORK/dilate2.json" \
  --variant "border30_recluster=$WORK/border30_recluster.json" \
  --variant "baysor=$WORK/baysor.json" \
  --out "$WORK/summary.json" --markdown "$WORK/summary.md"

RES="$DATA/results/celladmix"
mkdir -p "$RES"
cp "$WORK/summary.json" "$WORK/summary.md" "$RES/"
for name in vendor vendor_seed2 vendor_repeat border10 border30 dilate2 border30_recluster baysor; do
  cp "$WORK/$name.json" "$RES/validation_$name.json"
done
echo "=== done; results in $RES"
