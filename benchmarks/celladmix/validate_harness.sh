#!/usr/bin/env bash
# Re-validate the cellAdmix audit on the *harness* datasets (not on the
# self-made cellAdmix crop), following the review's requirements:
#
#   1. choose n_pool on one crop (xenium_lung_cancer_quick) and validate it
#      on a held-out crop (xenium_pancreas_377_full, ~39k cells);
#   2. check vendor < border10 < border30 < dilate2 with the fixed typing
#      (vendor quick-clustering transferred to every variant) and a fixed
#      pair set taken from the vendor audit (--fixed-pairs);
#   3. measure the audit SD across 3 Baysor replicates with baseline-
#      transferred typing (harness run.py --celltypes-from) — FIX-A1's
#      admixture tolerance input — and contrast it with per-replicate
#      quick clustering.
#
# Phases (run any subset):
#   ./validate_harness.sh --select        # n_pool candidates on the lung crop
#   ./validate_harness.sh --pancreas      # held-out chain with the chosen n_pool
#   ./validate_harness.sh --baysor-sd     # harness runs + SD measurement
#   ./validate_harness.sh --all
#
# Environment:
#   BENCH_PYTHON       bench env python   (default $BAYSOR_MAIN/.deps/bench/bin/python)
#   BAYSOR_BENCH_DATA  data root          (default $BAYSOR_MAIN/.bench-data)
#   BAYSOR_BIN         Baysor binary for the --baysor-sd phase
#   THREADS=6 SEED=1 NPOOL_CANDIDATES="20 40 60 80"
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAYSOR_MAIN=${BAYSOR_MAIN:-/home/vpetukhov/Projects/Baysor}
PY=${BENCH_PYTHON:-$BAYSOR_MAIN/.deps/bench/bin/python}
DATA=${BAYSOR_BENCH_DATA:-$BAYSOR_MAIN/.bench-data}
BAYSOR_BIN=${BAYSOR_BIN:-/home/vpetukhov/Projects/Baysor/.bench-data/binaries/baysor-bugfixes-35e8a7e}
THREADS=${THREADS:-6}
SEED=${SEED:-1}
NPOOL_CANDIDATES=${NPOOL_CANDIDATES:-"20 40 60 80"}
NPOOL_PANCREAS=${NPOOL_PANCREAS:-}          # set explicitly to skip --select
WORKBASE=$DATA/cache/celladmix/work/harness-validation
HARNESS=$ROOT/../harness
LUNG=$DATA/real/xenium_lung_cancer_quick
PANC=$DATA/real/xenium_pancreas_377_full
LUNG_WORK=$WORKBASE/lung_quick
PANC_WORK=$WORKBASE/pancreas_full

for f in "$LUNG/molecules.parquet" "$PANC/molecules.parquet"; do
  [[ -f "$f" ]] || { echo "missing dataset: $f" >&2; exit 1; }
done
mkdir -p "$WORKBASE"

audit() { # audit <workdir> <out.json> <extra args...>
  local wd=$1 out=$2; shift 2
  "$PY" -u "$ROOT/audit.py" --out "$out" --work-dir "$wd" \
    --threads "$THREADS" --seed "$SEED" "$@"
}

degrade() { # degrade <dataset dir> <workdir>
  local ds=$1 wd=$2
  mkdir -p "$wd"
  echo "=== degrade $(basename "$ds")"
  "$PY" "$ROOT/degrade.py" --molecules "$ds/molecules.parquet" \
    --cell-column cell_vendor --mode border --fraction 0.1 \
    --out "$wd/border10.assignment.parquet" --report "$wd/border10.degrade.json" \
    > "$wd/border10.degrade.stdout"
  "$PY" "$ROOT/degrade.py" --molecules "$ds/molecules.parquet" \
    --cell-column cell_vendor --mode border --fraction 0.3 \
    --out "$wd/border30.assignment.parquet" --report "$wd/border30.degrade.json" \
    > "$wd/border30.degrade.stdout"
  "$PY" "$ROOT/degrade.py" --molecules "$ds/molecules.parquet" \
    --cell-column cell_vendor --mode dilate --distance 2.0 \
    --out "$wd/dilate2.assignment.parquet" --report "$wd/dilate2.degrade.json" \
    > "$wd/dilate2.degrade.stdout"
}

extract_pairs() { # extract_pairs <vendor audit json> <out fixed_pairs.json>
  "$PY" - "$1" "$2" <<'EOF'
import json, sys
src, dst = sys.argv[1], sys.argv[2]
audit = json.load(open(src))
pairs = [{"source": p["source"], "target": p["target"]}
         for p in audit.get("pairs_top", [])]
if not pairs:
    raise SystemExit(f"{src}: no pairs_top; cannot build a fixed pair set")
json.dump(pairs, open(dst, "w"), indent=2)
print(f"fixed pair set: {len(pairs)} pairs -> {dst}")
EOF
}

run_chain() { # run_chain <dataset dir> <workdir> <npool>
  local ds=$1 wd=$2 npool=$3
  mkdir -p "$wd"
  echo "=== chain: $(basename "$ds") n_pool=$npool"
  # vendor baseline: quick-cluster typing, saved once for all variants
  audit "$wd/vendor.work" "$wd/vendor_npool$npool.json" \
    --molecules "$ds/molecules.parquet" --cell-column cell_vendor \
    --n-pool "$npool" --save-celltypes "$wd/vendor_celltypes.parquet"
  extract_pairs "$wd/vendor_npool$npool.json" "$wd/fixed_pairs_npool$npool.json"
  # every variant scored on the same pair set with the same typing
  audit "$wd/vendor.work" "$wd/vendor_fixed_npool$npool.json" \
    --molecules "$ds/molecules.parquet" --cell-column cell_vendor \
    --n-pool "$npool" --celltypes "$wd/vendor_celltypes.parquet" \
    --fixed-pairs "$wd/fixed_pairs_npool$npool.json"
  for v in border10 border30 dilate2; do
    audit "$wd/$v.work" "$wd/${v}_npool$npool.json" \
      --molecules "$ds/molecules.parquet" \
      --assignment "$wd/$v.assignment.parquet" \
      --celltypes "$wd/vendor_celltypes.parquet" \
      --fixed-pairs "$wd/fixed_pairs_npool$npool.json" \
      --n-pool "$npool"
  done
  "$PY" "$ROOT/summarize.py" \
    --variant "vendor=$wd/vendor_fixed_npool$npool.json" \
    --variant "border10=$wd/border10_npool$npool.json" \
    --variant "border30=$wd/border30_npool$npool.json" \
    --variant "dilate2=$wd/dilate2_npool$npool.json" \
    --title "cellAdmix audit on $(basename "$ds") (n_pool=$npool, fixed typing, fixed pairs)" \
    --out "$wd/summary_npool$npool.json" --markdown "$wd/summary_npool$npool.md" \
    > /dev/null
}

chain_margins() { # chain_margins <summary json> -> prints "<ok> <min_margin>"
  "$PY" - "$1" <<'EOF'
import json, sys
s = json.load(open(sys.argv[1]))
r = s["results"]
chain = ["vendor", "border10", "border30", "dilate2"]
vals = [r[n]["total_admixture_rate"] for n in chain]
if any(v is None for v in vals):
    print("unavailable", "nan")
else:
    margins = [b - a for a, b in zip(vals, vals[1:])]
    ok = all(m > 0 for m in margins)
    print(("monotone" if ok else "VIOLATED"), f"{min(margins):.6f}")
EOF
}

select_npool() {
  degrade "$LUNG" "$LUNG_WORK"
  local results="$WORKBASE/npool_selection.json"
  : > "$WORKBASE/npool_selection.txt"
  for n in $NPOOL_CANDIDATES; do
    run_chain "$LUNG" "$LUNG_WORK" "$n"
    read -r status margin <<<"$(chain_margins "$LUNG_WORK/summary_npool$n.json")"
    echo "n_pool=$n: $status min_margin=$margin" | tee -a "$WORKBASE/npool_selection.txt"
  done
  "$PY" - "$LUNG_WORK" "$WORKBASE/npool_selection.txt" \
      "$NPOOL_CANDIDATES" "$results" <<'EOF'
import json, sys
wd, txt, candidates, out = sys.argv[1:5]
# Selection criterion (two-stage, recorded for reproducibility):
#  1. hard: strictly monotone vendor < border10 < border30 < dilate2 chain;
#  2. hard: every fixed pair detected in every variant (comparability) and
#     mean relative spread of the fixed pairs' marker-pool `coverage` across
#     variants <= 0.10 — the documented failure mode of small n_pool (pool
#     cutoff reshuffles under perturbation, 20 -> 60 stabilises it);
#  3. tie-break: largest minimal consecutive margin, then smaller n_pool.
statuses = {}
for line in open(txt):
    parts = line.split()
    statuses[int(parts[0].split("=")[1].rstrip(":"))] = parts[1]
rows = []
for n in sorted(int(x) for x in candidates.split()):
    summary = json.load(open(f"{wd}/summary_npool{n}.json"))
    r = summary["results"]
    chain = ["vendor", "border10", "border30", "dilate2"]
    vals = [r[k]["total_admixture_rate"] for k in chain]
    margins = [b - a for a, b in zip(vals, vals[1:])] if all(v is not None for v in vals) else []
    monotone = bool(margins) and all(m > 0 for m in margins)
    fixed = [(p["source"], p["target"])
             for p in json.load(open(f"{wd}/fixed_pairs_npool{n}.json"))]
    cov, undetected = {}, set()
    for variant in ("vendor_fixed", "border10", "border30", "dilate2"):
        audit = json.load(open(f"{wd}/{variant}_npool{n}.json"))
        top = {(p["source"], p["target"]): p["coverage"]
               for p in audit.get("pairs_top", [])}
        for pair in fixed:
            if pair in top:
                cov.setdefault(pair, {})[variant] = top[pair]
            else:
                undetected.add(pair)
    spread = []
    for pair, d in cov.items():
        if len(d) == 4:
            vs = list(d.values())
            spread.append((max(vs) - min(vs)) / max(sum(vs) / len(vs), 1e-9))
    mean_spread = sum(spread) / len(spread) if spread else float("inf")
    rows.append({
        "n_pool": n,
        "chain": dict(zip(chain, vals)),
        "monotone": monotone,
        "min_margin": min(margins) if margins else None,
        "n_fixed_pairs": len(fixed),
        "fixed_pairs_undetected_somewhere": len(undetected),
        "mean_coverage_spread": round(mean_spread, 4),
        "passes": monotone and not undetected and mean_spread <= 0.10,
    })
passing = [row for row in rows if row["passes"]]
if not passing:
    raise SystemExit("no n_pool candidate passes the stability criterion")
best = sorted(passing, key=lambda row: (-row["min_margin"], row["n_pool"]))[0]
json.dump({
    "selection_crop": "real/xenium_lung_cancer_quick",
    "criterion": [
        "strictly monotone vendor < border10 < border30 < dilate2 chain",
        "every fixed pair detected in every variant and mean relative "
        "marker-pool coverage spread across variants <= 0.10 "
        "(small n_pool reshuffles the pool under perturbation)",
        "tie-break: largest minimal consecutive margin, then smaller n_pool",
    ],
    "candidates": rows,
    "chosen_n_pool": best["n_pool"],
    "chosen_min_margin": best["min_margin"],
    "held_out_crop": "real/xenium_pancreas_377_full",
}, open(out, "w"), indent=2)
print(f"chosen n_pool = {best['n_pool']} (min margin {best['min_margin']}, "
      f"coverage spread {best['mean_coverage_spread']})")
EOF
}

validate_pancreas() {
  local n=${NPOOL_PANCREAS:-}
  if [[ -z "$n" ]]; then
    n=$("$PY" -c "import json;print(json.load(open('$WORKBASE/npool_selection.json'))['chosen_n_pool'])")
  fi
  degrade "$PANC" "$PANC_WORK"
  run_chain "$PANC" "$PANC_WORK" "$n"
  chain_margins "$PANC_WORK/summary_npool$n.json"
}

baysor_sd() {
  local ds_id=${SD_DATASET:-xenium_lung_cancer_quick}
  local base=${SD_BASELINE:-a2-admx-$(echo "$ds_id" | tr 'A-Z' 'a-z' | cut -c1-24)}
  local common=(--baysor "$BAYSOR_BIN" --datasets "$ds_id" --threads "$THREADS")
  echo "=== baseline run (1 replicate; anchors the typing)"
  "$PY" "$HARNESS/run.py" "${common[@]}" --run-id "${base}-src" --replicates 1
  "$PY" "$HARNESS/baseline.py" create --run-id "${base}-src" --name "$base" \
    --allow-incomplete --force
  echo "=== 3 replicates with baseline-transferred typing + baseline pair set"
  "$PY" "$HARNESS/run.py" "${common[@]}" --run-id "${base}-sd3" --replicates 3 \
    --celltypes-from "$base"
  echo "=== 3 replicates with per-replicate quick clustering (old behaviour)"
  "$PY" "$HARNESS/run.py" "${common[@]}" --run-id "${base}-qc3" --replicates 3
  "$PY" - "$DATA" "$base" "$ds_id" "$THREADS" <<'EOF'
import json, os, statistics, sys
data, base, ds_id, threads = sys.argv[1:5]
def rates(run_id):
    m = json.load(open(f"{data}/runs/{run_id}/{ds_id}/metrics.json"))
    cm = m["real"]["celladmix"]
    return cm
trans = rates(f"{base}-sd3")
quick = rates(f"{base}-qc3")
def stat(cm):
    vals = [v for v in cm["per_rep_total"] if v is not None]
    return {"per_rep": cm["per_rep_total"],
            "mean": statistics.mean(vals) if vals else None,
            "sd": statistics.stdev(vals) if len(vals) > 1 else None,
            "typing": cm.get("typing"), "status": cm["status"]}
out = {
    "dataset": ds_id, "threads": int(threads), "replicates": 3,
    "transferred_typing": stat(trans),
    "quick_cluster_per_rep": stat(quick),
    "tolerance_input": {
        "recommended_admixture_tolerance": (
            round(3 * stat(trans)["sd"], 6)
            if stat(trans)["sd"] is not None else None),
        "rule": "3 x SD of the baseline-transferred audit over 3 replicates",
    },
}
dst = os.path.join(os.environ["WORKBASE"], "baysor_admixture_sd.json")
json.dump(out, open(dst, "w"), indent=2)
print(json.dumps(out, indent=2))
EOF
}

copy_results() {
  local res="$ROOT/results"
  mkdir -p "$res"
  local n=${NPOOL_PANCREAS:-}
  [[ -n "$n" ]] || n=$("$PY" -c "import json;print(json.load(open('$WORKBASE/npool_selection.json'))['chosen_n_pool'])")
  cp "$WORKBASE/npool_selection.json" "$res/harness_npool_selection.json"
  cp "$LUNG_WORK/summary_npool$n.md" "$res/harness_validation_lung.md"
  cp "$LUNG_WORK/summary_npool$n.json" "$res/harness_validation_lung.json"
  cp "$PANC_WORK/summary_npool$n.md" "$res/harness_validation_pancreas.md"
  cp "$PANC_WORK/summary_npool$n.json" "$res/harness_validation_pancreas.json"
  cp "$WORKBASE/baysor_admixture_sd.json" "$res/harness_baysor_sd.json" 2>/dev/null || true
  echo "results copied to $res"
}

phase=""
case "${1:-}" in
  --select)   phase=select ;;
  --pancreas) phase=pancreas ;;
  --baysor-sd) phase=sd ;;
  --copy)     phase=copy ;;
  --all)      phase=all ;;
  *) echo "usage: $0 --select|--pancreas|--baysor-sd|--copy|--all" >&2; exit 2 ;;
esac

export WORKBASE
case "$phase" in
  select)  select_npool ;;
  pancreas) validate_pancreas ;;
  sd)      baysor_sd ;;
  copy)    copy_results ;;
  all)     select_npool; validate_pancreas; baysor_sd; copy_results ;;
esac
