# Baysor benchmark harness

Runner, metrics, baselines and comparison for the Baysor benchmark suite
(see the dataset contract in [`../README.md`](../README.md)).

| file | purpose |
|---|---|
| `run.py` | run a Baysor binary over datasets, normalize output, compute metrics |
| `metrics.py` | pure metric functions (unit-tested, no I/O) |
| `baseline.py` | create/list committed baselines from a run |
| `compare.py` | compare a run with a baseline (`identical`/`same`/`improved`), Markdown + JSON report, exit code |
| `recompute_metrics.py` | recompute `metrics.json` from stored `assignment.parquet` files (no Baysor rerun) |
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

# refactor gate (default when --expect is omitted): bitwise-identical at 1 thread
$PY benchmarks/harness/run.py --baysor /path/to/baysor --datasets quick \
    --run-id refactor1 --replicates 1 --threads 1
$PY benchmarks/harness/baseline.py create --run-id refactor1 --name mybase-1t \
    --allow-incomplete
$PY benchmarks/harness/compare.py --run-id refactor2 --baseline mybase-1t

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

Primary metrics (these gate `--expect same` and `--expect improved`):

* **`accuracy_1to1`** — **PRIMARY** one-to-one (Hungarian) assignment
  accuracy: `scipy.optimize.linear_sum_assignment` on the predicted × true
  overlap matrix over all labels (0 included). A pair earns credit only when
  the two sides agree on noise status — predicted noise matches true noise
  (and nothing else), and a real predicted cell never earns credit from
  true-noise molecules; unmatched labels contribute 0. A pure split scores
  below 1.0 (only one part can be matched). This is the metric the
  st-recoverability oracle is defined against, so `oracle_gap` is computed
  from it.
* **`ari_assigned`** / **`ami_assigned`** — **PRIMARY** ARI/AMI over the
  molecules assigned in *both* sides. Unassigned molecules are excluded
  instead of collapsing into one giant cluster (on
  `sim_sparse_noisy_g100` this moves ARI from 0.661 to 0.993; on the vendor
  comparison for pancreas from 0.02 to 0.549).
* **`recovery_rate`** — fraction of true cells recovered with molecule-set
  Jaccard ≥ 0.5 (against predicted cells > 0; unmatched → 0).
* **`cell_count_ratio`** — #predicted cells / #true cells (labels > 0).

Reported separately (informational):

* **`matched_accuracy`** — the many-to-one majority metric: each predicted
  cell maps to the true cell it shares most molecules with (ties → smallest
  true label); a molecule is correct when its predicted cell maps to its own
  true cell. Predicted noise is only ever correct against true noise, and a
  real predicted cell whose **majority is true noise earns no credit** for
  its noise molecules. Several predicted cells may map to the same true cell
  (a pure split still scores 1.0 — that is why it is *secondary*).
* **`ari` / `ami`** — all-molecule values with label 0 as its own label.
* **`assigned_agreement`** — fraction of molecules agreeing on
  assigned-vs-unassigned status vs truth.
* **`assigned_fraction_pred` / `assigned_fraction_truth`** — the assigned
  fraction of each side.
* **`noise_precision` / `noise_recall`** — P(true noise | predicted noise)
  and P(predicted noise | true noise). **NaN (skipped) when the relevant
  noise count is below 50 molecules** (`NOISE_MIN_COUNT`), so they no longer
  fire on noise-free st-recoverability data; with `min_count=0` the vacuous
  1.0 behaviour is available for tiny hand-made cases.
* **`over_segmentation_rate`** — fraction of true cells whose largest
  predicted part holds < 80 % of the cell.
* **`under_segmentation_rate`** — fraction of predicted cells whose largest
  true source holds < 80 %.
* **`median_matched_jaccard`** — median, over true cells, of the best
  Jaccard against one predicted cell.
* **`oracle_gap`** — `meta.truth.oracle_accuracy − accuracy_1to1`, NaN when
  the generator provides no oracle.

**Difference vs st-recoverability's `matched_accuracy`**
(`$BAYSOR_BENCH_DATA/cache/sim/st-recoverability/src/headroom_common.py`):
that implementation scores molecules under the optimal *one-to-one*
matching of method cells to true cells (Hungarian on the contingency table),
with unassigned method labels participating like any other label. Our
`accuracy_1to1` follows the same one-to-one principle but zeroes the credit
of noise-status-mismatched pairs (a real cell never earns credit from noise
molecules and vice versa); the two agree exactly when there is no noise
(`tests/test_metrics.py::test_matched_accuracy_vs_st_recoverability_on_*`
asserts this on split/merge hand-made cases). The old many-to-one
`matched_accuracy` is kept only as a secondary metric.

### Real vs baseline segmentation

Both assignments are aligned by `mol_index` over the same input molecules.

Primary (gating):

* **`ari_assigned`** — ARI over molecules assigned in *both* segmentations
  (the primary agreement metric; the all-molecule `molecule_ari` is kept for
  reference only);
* **`frac_cells_matched`** — for each *baseline* cell the best molecule-set
  Jaccard in the run: fraction with J ≥ 0.5;
* **`cell_count_ratio`** — #run cells / #baseline cells.

Reported separately (informational):

* **`molecule_ari`** — ARI of the two full label vectors (0 = unassigned);
* **`assigned_agreement`** — fraction of molecules agreeing on
  assigned-vs-unassigned status;
* **`assigned_fraction_candidate` / `assigned_fraction_reference`** — the
  assigned fraction of each side;
* **`noise_precision` / `noise_recall`** — with the reference playing the
  role of truth; same 50-molecule rule (NaN below it);
* **`median_jaccard`** — the median baseline-cell Jaccard;
* **`median_mpc_rel_change`** — relative change of the median
  molecules-per-cell, `(run − baseline) / baseline`.

The same metrics against `cell_vendor` are stored in `real.vs_vendor` for
information only.

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
    `--threads` — `compare.py` **fails** on a thread mismatch in
    `same`/`improved` (and requires 1 thread in `identical`).

Therefore baselines use ≥ 3 replicates at the comparison thread count and
store per-metric mean and SD; `compare.py` uses
`tolerance = max(k·SD_pooled, floor)` with k = 3 and the SD **pooled per
metric across the baseline's datasets of the same kind** (a single
3-replicate SD has a 95% CI of [0.52σ, 6.3σ] and cannot be trusted
alone).

### Measured noise floor (6 threads, 3 replicates, this binary)

Numbers below are from `recompute_metrics.py` with the current metric
definitions (2026-09-29), i.e. what the committed baselines `harness-dev`
and `harness-real-dev` actually store.

Sim (mean ± sample SD over replicates):

| dataset | accuracy_1to1 | matched_accuracy | ari_assigned | recovery | cell_count_ratio |
|---|---|---|---|---|---|
| `sim_circles_gaps_g100` | 0.9865 ± 0.0002 | 0.9870 ± 0.0002 | 0.9928 ± 0.0005 | 1.0000 ± 0.0000 | 1.0353 ± 0.0020 |
| `sim_sparse_noisy_g100` | 0.9207 ± 0.0006 | 0.9286 ± 0.0003 | 0.9928 ± 0.0005 | 1.0000 ± 0.0000 | 1.1372 ± 0.0095 |
| `sim_tiled_distinct_g100` | 0.9896 ± 0.0002 | 0.9896 ± 0.0002 | 0.9994 ± 0.0003 | 1.0000 ± 0.0000 | 1.0217 ± 0.0023 |
| `strec_dense_s2_disjoint` | 0.5481 ± 0.0043 | 0.5606 ± 0.0039 | 0.3964 ± 0.0023 | 0.2240 ± 0.0200 | 1.0874 ± 0.0168 |

Note `ari_assigned` vs the all-molecule `ari` on `sim_sparse_noisy_g100`:
0.993 vs 0.66 — the gap is exactly the giant-unassigned-cluster artefact the
assigned-only primary metric removes.

Real, replicate-vs-replicate agreement (3 replicate pairs):

| dataset | ari_assigned | molecule_ari | assigned agreement | frac cells matched | median Jaccard | cell-count ratio |
|---|---|---|---|---|---|---|
| `xenium_pancreas_377_quick` | 0.7912 ± 0.0059 | 0.8009 ± 0.0058 | 0.9998 ± 0.0000 | 0.8719 ± 0.0047 | 0.8413 ± 0.0014 | 0.9981 ± 0.0017 |
| `xenium_breast_rep1_dense_quick` | 0.8254 ± 0.0123 | 0.8458 ± 0.0105 | 0.9992 ± 0.0000 | 0.9155 ± 0.0102 | 0.8583 ± 0.0088 | 1.0011 ± 0.0057 |
| `xenium_lung_cancer_quick` | 0.8281 ± 0.0017 | 0.8567 ± 0.0011 | 0.9996 ± 0.0001 | 0.8866 ± 0.0031 | 0.8295 ± 0.0042 | 0.9972 ± 0.0013 |

So at 6 threads the segmentations themselves are **not** run-to-run stable
(label-identity ARI ≈ 0.8) while assigned/unassigned status is essentially
stable (≥ 0.999).

## `compare.py`

```
compare.py --run-id R --baseline NAME [--expect {identical,same,improved}]
           [--k 3] [--admixture-tolerance FLOOR] [--report-md P] [--report-json P]
```

Exit code **0 = pass, 1 = fail, 2 = usage/setup error**. Reports are written
to `runs/<R>/compare_<NAME>_<expect>.{md,json}` and the Markdown is printed.
`--expect` defaults to **`identical`**, the default refactor gate.

### `--expect identical` (default; the refactor gate)

Baysor at 1 thread is bitwise-deterministic (findings above), so a change
that must not alter behaviour can be verified exactly:

* the run **and** the baseline must both have `threads == 1` and exactly
  1 replicate — each violation fails;
* for every dataset the `assignment_sha256` recorded per replicate in
  `metrics.json` is compared against the baseline's; any mismatch fails;
* on a mismatch the report shows the metric deltas (sim: every metric of
  `sim.mean`; real: `real.rep_agreement` deltas plus pair metrics between
  run rep0 and baseline rep0);
* dataset content hashes are checked as in `same`.

### `--expect same`

Provenance gates (fail, not warn):

* thread count must match between run and baseline, per dataset;
* the binary sha256 must match per dataset (`same` means *unchanged
  binary*; use `identical` for rebuilt-but-equivalent binaries);
* `inputs.molecules_sha256` / `inputs.meta_sha256` must match when
  recorded on both sides (recorded by the runner; `recompute_metrics.py`
  backfills them). Missing hashes → `skip` + warning, never a pass silently
  waved through.

Metric gates — only the **primary** metrics can fail the run, everything
else is informational:

| kind | primary metrics | floors |
|---|---|---|
| sim | `accuracy_1to1`, `ari_assigned`, `recovery_rate`, `cell_count_ratio` | 0.01 / 0.02 / 0.05 / 0.03 |
| real | `ari_assigned`, `frac_cells_matched`, `cell_count_ratio` | 0.02 / 0.05 / 0.03 |

* tolerance = `max(k·SD_pooled, floor)`, k = 3, SD pooled per metric across
  the baseline's datasets of the same kind (sim: replicate SDs of the
  means; real: replicate-pair agreement SDs);
* real checks measure the run-vs-baseline agreement (all run-rep ×
  baseline-rep pairs) against the baseline's replicate agreement (and the
  cell-count ratio against 1.0); a single-replicate baseline without
  replicate agreement falls back to absolute minima (ARI ≥ 0.90, matched ≥
  0.75) with a warning;
* **false-alarm budget**: every gated check carries its normal-approximation
  tail probability at the used tolerance; the report sums them
  (`false-alarm budget: ~0.009 expected false failures across 16 gated
  checks`) and warns when the budget exceeds 0.5.

### `--expect improved`

Improvement must exceed the noise:

* the mean gain in `accuracy_1to1` over the baseline's sim datasets must
  exceed **both** 2 standard errors (SE of the mean gain computed from the
  per-dataset replicate SDs, `SE = sqrt(Σ(sd_r²/n_r + sd_b²/n_b))/D`) **and**
  a minimum effect of 0.005;
* no individual sim dataset may regress beyond `max(k·SD_pooled, floor)` in
  `accuracy_1to1`, `ari_assigned`, `recovery_rate` or
  `over_segmentation_rate`;
* on real data: `total_admixture_rate ≤ baseline + max(k·SD, floor)`, where
  SD is taken over the *baseline's* audit replicates (`--admixture-tolerance`
  sets the floor, default 0). The check is gated **only** when the audit
  status is `ok` in both run and baseline *and* the baseline crop has
  ≥ 2000 cells; otherwise it is reported as `unavailable` (`skip`), never as
  0;
* an unchanged binary therefore cannot pass: +0.0001 mean gain ≪ 0.005
  (verified on `rev-same-sim`, see calibration);
* a baseline without sim datasets fails this mode — improvement cannot be
  demonstrated without sim truth.

All modes fail on: replicate failures/timeouts, datasets in the baseline
missing from the run, and (in `same`) missing baseline assignment tables.
Runtime and peak RSS are compared per dataset; a slowdown > 20 % or RSS
growth > 20 % adds a warning (never a failure). A binary sha256 change is
expected for `improved`, informational otherwise.

## Calibration (2026-09-29)

Runs used (all exist under `$BAYSOR_BENCH_DATA/runs/`, all recorded with the
same binary sha `68b1b505…`, 6 threads, 3 replicates; metrics recomputed
with the current definitions, baselines regenerated from `harness-val1` /
`harness-val3`):

| run | role |
|---|---|
| `harness-val1`, `harness-val3` | baseline runs → `harness-dev`, `harness-real-dev` |
| `rev-same-sim`, `rev-same-real` | unchanged binary, independent rerun |
| `rev-scale09-sim`, `rev-scale09-real` | subtle change (`--scale-factor 0.9`) |
| `harness-val2`, `harness-val4` | strong change (`--scale-factor 0.5`) |

Reproduce:

```bash
$PY benchmarks/harness/recompute_metrics.py --run harness-val1 --run harness-val2 \
    --run harness-val3 --run harness-val4 --run rev-same-sim --run rev-same-real \
    --run rev-scale09-sim --run rev-scale09-real
$PY benchmarks/harness/baseline.py create --run-id harness-val1 --name harness-dev --force
$PY benchmarks/harness/baseline.py create --run-id harness-val3 --name harness-real-dev --force
$PY benchmarks/harness/compare.py --run-id rev-same-sim --baseline harness-dev --expect same
```

Calibrated tolerances (`k = 3`; pooled SD from `harness-dev` /
`harness-real-dev`):

| kind | metric | pooled SD | floor | tolerance |
|---|---|---|---|---|
| sim | accuracy_1to1 | 0.00216 | 0.01 | 0.0100 |
| sim | ari_assigned | 0.00123 | 0.02 | 0.0200 |
| sim | recovery_rate | 0.01002 | 0.05 | 0.0500 |
| sim | cell_count_ratio | 0.00978 | 0.03 | 0.0300 |
| real | ari_assigned | 0.00795 | 0.02 | 0.0239 |
| real | frac_cells_matched | 0.00674 | 0.05 | 0.0500 |
| real | cell_count_ratio | 0.00350 | 0.03 | 0.0300 |

Measured outcomes with these numbers:

| comparison | expect | result |
|---|---|---|
| `rev-same-sim` vs `harness-dev` | `same` | **PASS** (32 gates, 0 fail) |
| `rev-same-real` vs `harness-real-dev` | `same` | **PASS** (21 gates, 0 fail) |
| `rev-same-sim` vs `harness-dev` | `improved` | **FAIL** — mean gain +0.00009 ≪ 0.005 (unchanged binary must not pass) |
| `rev-same-real` vs `harness-real-dev` | `improved` | **FAIL** — no sim datasets in the baseline |
| `rev-scale09-sim` vs `harness-dev` | `same` | **FAIL** — 2/4 datasets: `strec` (accuracy +0.023, cells +8.3 %), `sim_sparse_noisy` (cells +3.4 %) |
| `rev-scale09-real` vs `harness-real-dev` | `same` | **FAIL** — 3/3 datasets (assigned-ARI −0.028…−0.046, cells +7…+11 %) |
| `harness-val2` (scale 0.5) vs `harness-dev` | `same` | **FAIL** — 4/4 datasets, 14 gates |
| `harness-val4` (scale 0.5) vs `harness-real-dev` | `same` | **FAIL** — 3/3 datasets, 9 gates |

Scale × 0.9 fails `same` on **5 of 7 datasets** (all three real ones plus
two sim). The three trivial sim datasets show no reaction to scale × 0.9 at
all (deltas ≤ 0.0017, i.e. below the 6-thread replicate noise of 0.0006–0.003)
— catching those would require tolerances smaller than the noise floor
itself; use `--expect identical` at 1 thread to detect changes that small.
False-alarm budgets: ~0.009 over 16 sim gates, ~0.008 over 9 real gates.

## Validation performed

All on this machine with the Release binary at
`/home/vpetukhov/.bb/thread-storage/thr_cpwic2f6q3/baysor-bugfixes/build-rel/baysor`,
data root `/home/vpetukhov/Projects/Baysor/.bench-data`:

1. **Unit/fixture tests** (94): `cd benchmarks/harness/tests && $PY -m pytest -q` —
   hand-made metric cases (perfect, one split, one merge, all noise,
   permuted labels, majority-noise rule, 50-molecule noise rule,
   assigned-only ARI, one-to-one vs many-to-one, st-recoverability
   divergence), command construction, output alignment (positional,
   shuffled, loader-dropped, legacy CSV), dataset selection,
   `/usr/bin/time` parsing, timeouts, baseline creation, `recompute_metrics`,
   and end-to-end `compare.py` tests for **all three modes** on synthetic
   metric JSONs (identical sha pass/mismatch/thread/replicate gates;
   same-mode primary gates, provenance/content-hash gates, pooled tolerances,
   false-alarm budget; improved-mode 2-SE + min-effect aggregate, regression
   gates, admixture gating conditions) plus an end-to-end pipeline run
   against the real binary.
2. **Calibration** — the table in the section above: unchanged runs pass
   `same` and fail `improved`, scale × 0.9 fails `same` on most datasets,
   scale × 0.5 fails `same` everywhere.
3. **Determinism** — `harness-det1a`/`harness-det1b`: runs at `--threads 1`
   on `xenium_pancreas_377_quick`, bitwise-identical assignments (see
   findings above); this is what the `identical` mode gates on.

## Tests

```bash
cd benchmarks/harness/tests
/home/vpetukhov/Projects/Baysor/.deps/bench/bin/python -m pytest -q
```

Tests generate their own contract-conformant fixtures in tmp dirs
(`fixtures.py`); the end-to-end test additionally needs the Baysor binary
(`BAYSOR_BIN` overrides the default path) and is skipped when it is absent.
