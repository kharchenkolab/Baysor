# Baysor benchmark harness

Runner, metrics, baselines and comparison for the Baysor benchmark suite
(see the dataset contract in [`../README.md`](../README.md)).

| file | purpose |
|---|---|
| `run.py` | run a Baysor binary over datasets, normalize output, compute metrics |
| `metrics.py` | pure metric functions (unit-tested, no I/O) |
| `baseline.py` | create/list committed baselines from a run |
| `compare.py` | compare a run with a baseline, Markdown + JSON report, exit code |
| `bench.sh` | one command: run → compare → print report |
| `celladmix.py` | optional adapter for the cellAdmix audit (`../celladmix/audit.py`) |
| `tests/` | pytest suite incl. contract-conformant fixture datasets |

Python environment: `.deps/bench/bin/python` (see `../environment.yml`,
which now also lists `pytest`). All commands below use it as `$PY`.

## Quick start

```bash
export BAYSOR_BENCH_DATA=/path/to/.bench-data    # shared data root

# 1. run the quick tier, 3 replicates (noise floor), 6 threads
$PY benchmarks/harness/run.py --baysor /path/to/baysor \
    --datasets quick --run-id myrun --replicates 3

# 2. freeze it as a baseline (metrics committed, assignments kept in data root)
$PY benchmarks/harness/baseline.py create --run-id myrun --name mybase

# 3. later: check an unchanged algorithm
$PY benchmarks/harness/compare.py --run-id otherrun --baseline mybase --expect same

# or everything in one command
benchmarks/harness/bench.sh --baysor /path/to/baysor --baseline mybase \
    --expect same --replicates 3
```

`bench.sh --create-baseline NAME` bootstraps a baseline after the run instead
of comparing. `BENCH_PY` overrides the interpreter, `BAYSOR_BIN` the binary.

## Dataset selection (`--datasets`)

* `quick`, `full`, `all` — select by the `tier` field of each dataset's
  `meta.json` (the manifests under `../datasets/` mirror it);
* otherwise a comma-separated list of dataset ids and/or shell globs
  (`sim_tiled_*`, `strec_dense_s2_disjoint`, ...);
* `--kind sim|real` restricts to `$BAYSOR_BENCH_DATA/{sim,real}`.

## `run.py`

```
run.py --baysor PATH --datasets SPEC --run-id ID
       [--kind sim|real] [--threads 6] [--replicates 1] [--timeout S]
       [--data-root PATH] [--label SHA] [--no-celladmix] [--skip-existing]
       [--scale-factor F]
```

For each dataset × replicate the command is built from `meta.json`:

| meta field | flag |
|---|---|
| `molecules.parquet` | positional `coordinates`; `-x x -y y [-z z] -g gene` (`-z` only when the input has a `z` column) |
| `baysor.prior: "column"` | second positional `:prior` |
| `baysor.prior: "image:<rel>"` | second positional `<dataset>/<rel>` |
| `baysor.prior_confidence` | `--prior-segmentation-confidence` (only when a prior is used) |
| `baysor.scale_um` | `-s` (times `--scale-factor`) |
| `baysor.scale_std` | `--scale-std` |
| `baysor.min_molecules_per_cell` | `-m` |
| `baysor.config` | `-c` resolved **relative to the repository root** |
| `baysor.extra_args` | appended **verbatim after `-c`** (e.g. the Xenium datasets' `-x x -y y -g gene --qv-column qv --unassigned-prior-label 0`, which override `configs/xenium.toml`'s vendor column mapping and `unassigned_label="UNASSIGNED"` — without `--unassigned-prior-label 0` prior label 0 would become a real segment). The builder **omits its own occurrence** of any flag `extra_args` already provides, so each option appears exactly once (the CLI rejects repeated scalar options). Baysor loads the config as *defaults* before CLI11 parsing, so explicit CLI flags win regardless of position (`src/cli/main.cpp`), and `extra_args` after `-c` makes the intent explicit. |

The exact flag names are probed once per invocation from
`baysor run --help` of the given binary (`--output-style parquet` is requested
when supported; otherwise the legacy `segmentation.csv` is parsed).
`--threads N` has no CLI equivalent in Baysor, so the runner sets
`OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`/`NUMEXPR_NUM_THREADS`
to `N`. `--scale-factor` exists solely to build deliberately degraded runs
for validation; keep it at 1.0.

Each replicate executes under `/usr/bin/time -v` in its own process group
(killed on `--timeout`) and stores:

```
$BAYSOR_BENCH_DATA/runs/<run_id>/_binary.json          # sha256, baysor --help, label
$BAYSOR_BENCH_DATA/runs/<run_id>/<dataset>/rep<k>/
    seg/                     # raw Baysor output (parquet style)
    assignment.parquet       # normalized per-molecule assignment
    run.json                 # command, exit code, wall time, peak RSS,
                             # binary sha256, version info, --label git SHA
    celladmix.json           # real datasets only, when the audit exists
$BAYSOR_BENCH_DATA/runs/<run_id>/<dataset>/metrics.json
```

**Assignment table.** `assignment.parquet` has exactly three columns:
`mol_index` (int64, row order of the input `molecules.parquet`), `cell`
(int64, **0 = unassigned/noise**) and `confidence` (float64, Baysor's
`assignment_confidence`, NaN when absent). Molecules Baysor's loader drops
(`min_molecules_per_gene`, `exclude_genes`, `min_qv`, coordinate bounds —
e.g. `configs/xenium.toml` drops genes with < 10 molecules) keep cell 0 and
NaN confidence; `run.json` records `n_loader_filtered`. Row mapping is
positional (verified against gene + coordinates) with a coordinate/gene join
as fallback; output rows that match no input molecule abort the run.

**`metrics.json`** aggregates per dataset: provenance (binary sha256, label,
threads, replicates), runtime (wall/RSS per rep + mean/SD), sim metrics
(`sim.per_rep/mean/sd`), real replicate-vs-replicate agreement over all
replicate pairs (`real.rep_agreement`), the same metrics against
`cell_vendor` (`real.vs_vendor`, information only), the cellAdmix audit
(`real.celladmix`) and any failures/timeouts.

## Metrics (`metrics.py`)

### Sim vs truth — computed over molecules with `interior == True`

* **`matched_accuracy`** — each predicted cell (including 0 = noise) is
  mapped to the true cell it shares most molecules with (ties → smallest
  true label); a molecule is correct when its predicted cell maps to its own
  true cell. Predicted noise is only ever correct against true noise
  (noise ↔ noise). Many-to-one: several predicted cells may map to the same
  true cell.
* **`ari` / `ami`** — adjusted Rand / adjusted mutual information with label 0
  as its own label.
* **`noise_precision` / `noise_recall`** — P(true noise | predicted noise)
  and P(predicted noise | true noise); vacuously 1.0 when the denominator is 0.
* **`cell_count_ratio`** — #predicted cells / #true cells (labels > 0).
* **`over_segmentation_rate`** — fraction of true cells whose largest
  predicted part holds < 80 % of the cell.
* **`under_segmentation_rate`** — fraction of predicted cells whose largest
  true source holds < 80 %.
* **`recovery_rate`** — fraction of true cells recovered with molecule-set
  Jaccard ≥ 0.5 (against predicted cells > 0; unmatched → 0).
* **`median_matched_jaccard`** — median, over true cells, of the best
  Jaccard against one predicted cell.
* **`oracle_gap`** — `meta.truth.oracle_accuracy − matched_accuracy`, NaN
  when the generator provides no oracle.

**Difference vs st-recoverability's `matched_accuracy`**
(`$BAYSOR_BENCH_DATA/cache/sim/st-recoverability/src/headroom_common.py`):
that implementation scores molecules under the *optimal one-to-one*
matching of method cells to true cells (Hungarian on the contingency table),
with unassigned method labels participating like any other label (their
docstring: background/unassigned counts as errors). Ours is *many-to-one
majority* with an explicit noise↔noise rule. Consequences:

* a pure split (two predicted parts of one true cell) scores 1.0 here but
  ≤ (largest part) under one-to-one matching — only one part can be matched
  (`tests/test_metrics.py::test_matched_accuracy_vs_st_recoverability_on_split`
  asserts 1.0 vs 34/38 on a hand-made case);
* merged cells are penalized by both definitions (majority mapping vs
  one-to-one) and they agree on clean one-to-one situations
  (`..._on_merge_agree`);
* predicted noise is correct against true noise here, never against a real
  cell; there, unassigned labels cannot earn credit except through the
  single slot the Hungarian assignment gives label 0.

Both operate on interior molecules only.

### Real vs baseline segmentation

Both assignments are aligned by `mol_index` over the same input molecules:

* **`molecule_ari`** — ARI of the two label vectors (0 = unassigned);
* **`assigned_agreement`** — fraction of molecules agreeing on
  assigned-vs-unassigned status;
* **`frac_cells_matched`** / **`median_jaccard`** — for each *baseline*
  cell the best molecule-set Jaccard in the run: fraction with J ≥ 0.5 and
  the median;
* **`cell_count_ratio`** — #run cells / #baseline cells;
* **`median_mpc_rel_change`** — relative change of the median
  molecules-per-cell, `(run − baseline) / baseline`.

The same five metrics against `cell_vendor` are stored in
`real.vs_vendor` for information only.

### Real admixture (cellAdmix)

Optional and pluggable: when `benchmarks/celladmix/audit.py` (BENCH-CELLADMIX)
exists, every real replicate is audited as
`python .../audit.py --molecules <parquet> --assignment <parquet> --out <json>
[--image ...]` and `total_admixture_rate` is stored per replicate and as a
mean. `--no-celladmix` disables it; a missing or crashing module degrades to
`{"status": "absent" | "failed"}` and comparisons skip the check with a
warning instead of failing.

## Baselines and the noise floor

```bash
$PY benchmarks/harness/baseline.py create --run-id R --name NAME [--force]
$PY benchmarks/harness/baseline.py list
```

* metrics are copied to `../baselines/NAME/<dataset>.json` (committed, ~12–17 KB);
* assignment tables are copied to
  `$BAYSOR_BENCH_DATA/baselines/NAME/<dataset>/rep<k>/assignment.parquet`
  (never committed; their sha256 is recorded in the JSON);
* a baseline normally requires ≥ 3 replicates (`--allow-incomplete`
  overrides, for deterministic 1-thread runs and fixtures), because the
  binary is stochastic above one thread (below).

### Determinism findings (code inspection + experiment, 2026-09-29)

* **There is no seed option.** `baysor run --help` exposes no seed/random-seed
  flag.
* The RNGs are hard-wired: one global `Xoshiro256pp` seeded with `1`
  (`src/utils/general.cpp`; `reset_global_xoshiro_rng` exists but is never
  called), per-thread RNGs seeded from the thread index
  (`1 ^ t·2654435761`, `src/processing/bmm_algorithm/bmm_algorithm.cpp`),
  `std::mt19937(42)` in molecule clustering and neighborhood composition.
  The E-step is *always* stochastic (`expect_dirichlet_spatial(data, /*stochastic=*/true)`).
* **1 thread → deterministic.** Verified empirically: three runs at
  `OMP_NUM_THREADS=1` on a 47k-molecule fixture and three runs at 1 thread on
  the 130k-molecule `xenium_pancreas_377_quick` (across two separate
  invocations) produced **bitwise-identical** assignment tables.
* **>1 thread → not deterministic.** The E-step is
  `#pragma omp parallel for schedule(dynamic, 1024)` over molecules, each
  thread drawing from its own fixed-seed stream; scheduling decides which
  molecule sees which stream, and scheduling is timing-dependent. Measured
  on identical repeated runs at 6 threads:
  * fixture (47k mol): replicate ARI 0.69–0.76, cell counts 12→14;
  * `xenium_pancreas_377_quick`: replicate ARI 0.79–0.81;
  * results also differ *between* thread counts (fixture ARI 0.62 between
    1 and 6 threads), so a baseline must be reproduced with the same
    `--threads` — `compare.py` warns on a thread mismatch.

Therefore baselines use ≥ 3 replicates at the comparison thread count and
store per-metric mean and SD; `compare.py` uses
`tolerance = max(3·SD_baseline, floor)`.

### Measured noise floor (6 threads, 3 replicates, this binary)

Sim (mean ± sample SD over replicates):

| dataset | matched_accuracy | ari | cell_count_ratio | recovery_rate |
|---|---|---|---|---|
| `sim_circles_gaps_g100` | 0.9872 ± 0.0002 | 0.9611 ± 0.0004 | 1.0353 ± 0.0020 | 1.000 ± 0.000 |
| `sim_sparse_noisy_g100` | 0.9292 ± 0.0003 | 0.6603 ± 0.0014 | 1.1372 ± 0.0095 | 1.000 ± 0.000 |
| `sim_tiled_distinct_g100` | 0.9898 ± 0.0002 | 0.9608 ± 0.0003 | 1.0217 ± 0.0023 | 1.000 ± 0.000 |
| `strec_dense_s2_disjoint` | 0.5606 ± 0.0039 | 0.3964 ± 0.0023 | 1.0874 ± 0.0168 | 0.224 ± 0.020 |

The largest SDs across all sim metrics are 0.02 (`recovery_rate` on the hard
`strec_dense` dataset); most are ≤ 0.01, i.e. comfortably inside the
per-metric absolute floors (0.01–0.05) used by `--expect same`.

Real, replicate-vs-replicate agreement (3 replicate pairs):

| dataset | molecule ARI | assigned agreement | frac cells matched | median Jaccard | cell-count ratio |
|---|---|---|---|---|---|
| `xenium_pancreas_377_quick` | 0.8009 ± 0.0058 | 0.9998 | 0.8719 ± 0.0047 | 0.8413 | 0.9981 ± 0.0017 |
| `xenium_breast_rep1_dense_quick` | 0.8458 ± 0.0105 | 0.9992 | 0.9155 ± 0.0102 | 0.8583 | 1.0011 ± 0.0057 |
| `xenium_lung_cancer_quick` | 0.8567 ± 0.0011 | 0.9996 | 0.8866 ± 0.0031 | 0.8295 | 0.9972 ± 0.0013 |

So at 6 threads the segmentations themselves are **not** run-to-run stable
(label-identity ARI ≈ 0.8) while assigned/unassigned status is essentially
stable (≥ 0.999). The `--expect same` margins on real data
(ARI −0.05, assigned −0.02, matched cells −0.10, Jaccard −0.10 below the
replicate agreement) are drawn around this floor.

## `compare.py`

```
compare.py --run-id R --baseline NAME --expect {same,improved}
           [--k 3] [--admixture-tolerance 0.01] [--report-md P] [--report-json P]
```

Exit code **0 = pass, 1 = fail, 2 = usage/setup error**. Reports are written
to `runs/<R>/compare_<NAME>_<expect>.{md,json}` and the Markdown is printed.

`--expect same` — unchanged algorithm must stay inside the noise floor:

* every sim metric: `|run_mean − baseline_mean| ≤ max(k·SD_baseline, floor)`
  with k = 3 and per-metric floors (matched accuracy/oracle gap 0.01, ARI/AMI
  0.02, noise P/R 0.02, cell-count ratio 0.05, over/under-segmentation 0.05,
  recovery 0.05, median Jaccard 0.05);
* real: run-vs-baseline agreement (all run-rep × baseline-rep pairs) must be
  ≥ baseline replicate-vs-replicate agreement − margin
  (ARI −0.05, assigned −0.02, matched cells −0.10, Jaccard −0.10);
  cell-count ratio within ±max(0.10, 3·SD) of 1 and median-molecules-per-cell
  change within ±max(0.10, 3·SD) of 0 (10 % either way). A single-replicate
  baseline without a replicate-agreement floor falls back to absolute minima
  (0.90/0.95/0.75/0.75) with a warning.

`--expect improved` — changed algorithm must get better:

* mean sim matched accuracy over the baseline's sim datasets must **increase**
  (strictly), and no individual sim dataset may drop by more than
  `max(3·SD, 0.01)`;
* every real dataset: `total_admixture_rate ≤ baseline + tolerance`
  (default +0.01); unavailable audits are skipped with a warning;
* a baseline without sim datasets fails this mode — improvement cannot be
  demonstrated without sim truth.

Both modes fail on: replicate failures/timeouts, datasets in the baseline
missing from the run, and missing baseline assignment tables. Runtime and
peak RSS are compared per dataset; a slowdown > 20 % or RSS growth > 20 %
adds a warning (never a failure). A binary sha256 change is noted (expected
for algorithm changes, suspicious for `same`).

## Validation performed

All on this machine with the Release binary at
`/home/vpetukhov/.bb/thread-storage/thr_cpwic2f6q3/baysor-bugfixes/build-rel/baysor`,
data root `/home/vpetukhov/Projects/Baysor/.bench-data`:

1. **Unit tests** (72): `cd benchmarks/harness/tests && $PY -m pytest -q` —
   hand-made metric cases (perfect, one split, one merge, all noise,
   permuted labels, st-recoverability divergence), command construction,
   output alignment (positional, shuffled, loader-dropped, legacy CSV),
   dataset selection, `/usr/bin/time` parsing, timeouts, baseline creation,
   both compare modes incl. admixture tolerances, and an end-to-end pipeline
   run against the real binary.
2. **Sim pipeline** — `harness-val1`: 4 quick sim datasets × 3 reps
   (`strec_dense_s2_disjoint`, `sim_sparse_noisy_g100`,
   `sim_tiled_distinct_g100`, `sim_circles_gaps_g100`, 48k–116k molecules,
   ~30 s per run at 6 threads) → baseline `harness-dev` →
   `--expect same` against itself: **PASS** (41 checks);
   `--expect same` against a scale-halved run (`harness-val2`,
   `--scale-factor 0.5`): **FAIL** with 21 flagged checks. Note the halved
   scale *raised* matched accuracy on the dense strec dataset (+0.086), so
   that run legitimately passes `--expect improved`; a 20 %-random-reassign
   degradation fails both modes (covered by `tests/test_pipeline.py`).
3. **Real pipeline** — `harness-val3`: 3 quick Xenium datasets × 3 reps
   (130k molecules, `prior: column`, ~35 s per run) → baseline
   `harness-real-dev` → `--expect same` against itself: **PASS**
   (18 checks); scale-halved run `harness-val4`: **FAIL** with 15 flagged
   checks (ARI 0.59–0.65 vs required ≥ 0.796, cell counts 2.2×).
   `--expect improved` on the real-only baseline fails with
   "baseline contains no sim datasets" (by design) and skips the absent
   cellAdmix audit gracefully.
4. **Determinism** — `harness-det1a`/`harness-det1b`: 3 runs at `--threads 1`
   on `xenium_pancreas_377_quick`, bitwise-identical assignments
   (see findings above).

Reproduce:

```bash
export BAYSOR_BENCH_DATA=/home/vpetukhov/Projects/Baysor/.bench-data
B=/home/vpetukhov/.bb/thread-storage/thr_cpwic2f6q3/baysor-bugfixes/build-rel/baysor
$PY benchmarks/harness/run.py --baysor $B \
    --datasets strec_dense_s2_disjoint,sim_sparse_noisy_g100,sim_tiled_distinct_g100,sim_circles_gaps_g100 \
    --run-id harness-val1 --replicates 3
$PY benchmarks/harness/baseline.py create --run-id harness-val1 --name harness-dev --force
$PY benchmarks/harness/compare.py --run-id harness-val1 --baseline harness-dev --expect same
$PY benchmarks/harness/run.py --baysor $B --kind real \
    --datasets xenium_pancreas_377_quick,xenium_breast_rep1_dense_quick,xenium_lung_cancer_quick \
    --run-id harness-val3 --replicates 3
$PY benchmarks/harness/baseline.py create --run-id harness-val3 --name harness-real-dev --force
$PY benchmarks/harness/compare.py --run-id harness-val3 --baseline harness-real-dev --expect same
```

`harness-dev` and `harness-real-dev` are committed under
`../baselines/` as working examples; replace them once the final
baseline run is agreed on.

## Tests

```bash
cd benchmarks/harness/tests
/home/vpetukhov/Projects/Baysor/.deps/bench/bin/python -m pytest -q
```

Tests generate their own contract-conformant fixtures in tmp dirs
(`fixtures.py`); the end-to-end test additionally needs the Baysor binary
(`BAYSOR_BIN` overrides the default path) and is skipped when it is absent.
