#!/usr/bin/env bash
# Run the Release Baysor binary on a built benchmark dataset and record
# wall time, peak RSS and output presence (BENCH-REALX sanity check).
#
# Usage: benchmarks/fetch/sanity_run.sh <dataset_id>
# Env:   BAYSOR_BIN, BAYSOR_BENCH_DATA, RUN_ID (default sanity_realx),
#        N_THREADS (default 6, the shared-machine limit for this suite).
set -euo pipefail

BAYSOR_BIN=${BAYSOR_BIN:-/home/vpetukhov/.bb/thread-storage/thr_cpwic2f6q3/baysor-bugfixes/build-rel/baysor}
BAYSOR_BENCH_DATA=${BAYSOR_BENCH_DATA:-/home/vpetukhov/Projects/Baysor/.bench-data}
BENCH_PY=${BENCH_PY:-/home/vpetukhov/Projects/Baysor/.deps/bench/bin/python}
RUN_ID=${RUN_ID:-sanity_realx}
N_THREADS=${N_THREADS:-6}
export N_THREADS

id=${1:?usage: sanity_run.sh <dataset_id>}
ds=$BAYSOR_BENCH_DATA/real/$id
out=$BAYSOR_BENCH_DATA/runs/$RUN_ID/$id
test -f "$ds/molecules.parquet" || { echo "missing $ds/molecules.parquet" >&2; exit 1; }
mkdir -p "$out"

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
cd "$repo_root"

# scale_um from the dataset meta, as the harness would pass it
scale=$("$BENCH_PY" -c "import json;print(json.load(open('$ds/meta.json'))['baysor']['scale_um'])")

echo "== baysor run: $id -> $out (OMP_NUM_THREADS=$N_THREADS, scale=$scale)"
OMP_NUM_THREADS=$N_THREADS /usr/bin/time -v -o "$out/time.txt" \
  "$BAYSOR_BIN" run -c configs/xenium.toml \
  -x x -y y -g gene --qv-column qv --unassigned-prior-label 0 --scale "$scale" \
  -o "$out" "$ds/molecules.parquet" :prior \
  >"$out/stdout.log" 2>"$out/stderr.log"

# Extract wall time and peak RSS into timing.json for the inventory report.
"$BENCH_PY" - "$out" <<'PYEOF'
import json, re, sys
from pathlib import Path

out = Path(sys.argv[1])
txt = (out / "time.txt").read_text()
m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", txt)
parts = [float(x) for x in m.group(1).split(":")]
wall_s = sum(v * 60**i for i, v in enumerate(reversed(parts)))
r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt)
timing = {
    "wall_s": round(wall_s, 2),
    "max_rss_kb": int(r.group(1)) if r else None,
    "threads": int(__import__("os").environ.get("N_THREADS", 6)),
}
(out / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
print(timing)
PYEOF

for f in segmentation.csv segmentation_polygons_2d.json segmentation_log.log; do
  test -s "$out/$f" || { echo "MISSING or empty output: $out/$f" >&2; exit 1; }
done
echo "OK: outputs present in $out"
