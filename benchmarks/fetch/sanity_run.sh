#!/usr/bin/env bash
# Run the Release Baysor binary on a built benchmark dataset and record
# wall time, peak RSS and output presence (BENCH-REALX sanity check).
#
# The command line is assembled from the dataset's meta.json the way the
# harness assembles it: config, extra_args (e.g. -z z for 3D crops), scale,
# prior (`column` -> :prior, `image:<path>` -> the label TIFF, `none` ->
# no positional), prior confidence and min-molecules-per-cell.
#
# Usage: benchmarks/fetch/sanity_run.sh <dataset_id>
# Env:   BAYSOR_BIN, BAYSOR_BENCH_DATA, RUN_ID (default sanity_realx),
#        N_THREADS (default 6, the shared-machine limit for this suite).
set -euo pipefail

BAYSOR_BIN=${BAYSOR_BIN:-/home/vpetukhov/Projects/Baysor/.bench-data/binaries/baysor-bugfixes-35e8a7e}
BAYSOR_BENCH_DATA=${BAYSOR_BENCH_DATA:-/home/vpetukhov/Projects/Baysor/.bench-data}
BENCH_PY=${BENCH_PY:-/home/vpetukhov/Projects/Baysor/.deps/bench/bin/python}
RUN_ID=${RUN_ID:-sanity_realx}
N_THREADS=${N_THREADS:-6}
export N_THREADS

id=${1:?usage: sanity_run.sh <dataset_id>}
ds=$BAYSOR_BENCH_DATA/real/$id
out=$BAYSOR_BENCH_DATA/runs/$RUN_ID/$id
test -f "$ds/molecules.parquet" || { echo "missing $ds/molecules.parquet" >&2; exit 1; }
test -f "$ds/meta.json" || { echo "missing $ds/meta.json" >&2; exit 1; }
mkdir -p "$out"

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
cd "$repo_root"

# options and prior positional, derived from meta.json (one per line)
mapfile -t baysor_opts < <("$BENCH_PY" - "$ds" <<'PYEOF'
import json, sys
from pathlib import Path

cfg = json.load(open(Path(sys.argv[1]) / "meta.json")).get("baysor", {})
opts = []
if cfg.get("config"):
    opts += ["-c", str(cfg["config"])]
if cfg.get("scale_um") is not None:
    opts += ["--scale", str(cfg["scale_um"])]
if cfg.get("scale_std") is not None:
    opts += ["--scale-std", str(cfg["scale_std"])]
if cfg.get("min_molecules_per_cell") is not None:
    opts += ["-m", str(cfg["min_molecules_per_cell"])]
prior = cfg.get("prior") or "none"
if prior != "none" and cfg.get("prior_confidence") is not None:
    opts += ["--prior-segmentation-confidence", str(cfg["prior_confidence"])]
opts += [str(a) for a in (cfg.get("extra_args") or [])]
print("\n".join(opts))
PYEOF
)
if [ ${#baysor_opts[@]} -eq 0 ]; then
  echo "failed to derive baysor options from $ds/meta.json" >&2
  exit 1
fi
prior_arg=$("$BENCH_PY" - "$ds" <<'PYEOF'
import json, sys
from pathlib import Path

ds = Path(sys.argv[1])
prior = json.load(open(ds / "meta.json")).get("baysor", {}).get("prior") or "none"
if prior == "column":
    print(":prior")
elif prior.startswith("image:"):
    img = ds / prior[len("image:"):]
    if not img.is_file():
        raise SystemExit(f"prior image not found: {img}")
    print(img)
elif prior in ("none", ""):
    print("")
else:
    raise SystemExit(f"unsupported baysor.prior: {prior!r}")
PYEOF
)

echo "== baysor run: $id -> $out (OMP_NUM_THREADS=$N_THREADS)"
cmd=("$BAYSOR_BIN" run "${baysor_opts[@]}" -o "$out" "$ds/molecules.parquet")
if [ -n "$prior_arg" ]; then
  cmd+=("$prior_arg")
fi
OMP_NUM_THREADS=$N_THREADS /usr/bin/time -v -o "$out/time.txt" \
  "${cmd[@]}" >"$out/stdout.log" 2>"$out/stderr.log"

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
