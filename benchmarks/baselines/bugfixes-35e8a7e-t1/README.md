# Baseline: `bugfixes-35e8a7e-t1` (exact, `identical` flavour)

Official 1-thread reference baseline of the Baysor benchmark suite for the
**Release build of branch `bugfixes` at commit `35e8a7e`**
(`--label bugfixes-35e8a7e`, binary
`/home/vpetukhov/.bb/thread-storage/thr_cpwic2f6q3/baysor-bugfixes/build-rel/baysor`),
created by task BENCH-BASELINE.

Use it with `compare.py --expect identical` (the refactor gate): at
`OMP_NUM_THREADS=1` Baysor is bitwise-deterministic, so any changed
`assignment_sha256` is a real behaviour change.

## Contents

* **65 quick-tier datasets** (49 simulated + 16 real), 1 thread,
  1 replicate, **no cellAdmix audit**;
* one `<dataset>.json` per dataset — the full `metrics.json` of run
  `benchbase-t1` plus a `baseline` block (flavour `identical`, the sha256 of
  the stored assignment table, binary sha256/label, input content hashes);
* assignment tables (not committed) under
  `$BAYSOR_BENCH_DATA/baselines/bugfixes-35e8a7e-t1/<dataset>/rep0/`.

## Reproduction

```bash
export BAYSOR_BENCH_DATA=/home/vpetukhov/Projects/Baysor/.bench-data
PY=.deps/bench/bin/python
B=/home/vpetukhov/.bb/thread-storage/thr_cpwic2f6q3/baysor-bugfixes/build-rel/baysor

# run (phase 1 sampled 4 datasets, phase 2 completed the tier with
# --skip-existing; a single invocation with --datasets quick works too)
$PY benchmarks/harness/run.py --baysor $B --datasets quick \
    --run-id benchbase-t1 --replicates 1 --threads 1 --timeout 1800 \
    --no-celladmix --label bugfixes-35e8a7e
# freeze
$PY benchmarks/harness/baseline.py create --run-id benchbase-t1 \
    --name bugfixes-35e8a7e-t1 --identical
```

Per-dataset timeouts: 1800 s (30 min); no dataset timed out or failed
(`benchbase-t1` finished `OK`, 65/65 replicates ok).

## Self-check

A fresh 1-thread, 1-replicate run of a 10-dataset subset
(`selfcheck-t1-fresh`) must be **bitwise identical** to this baseline:

```bash
$PY benchmarks/harness/compare.py --run-id selfcheck-t1-fresh \
    --baseline bugfixes-35e8a7e-t1 --expect identical     # expect PASS
```

A `--scale-factor 0.9` run of the same subset (`selfcheck-t1-scale09`) must
**fail** the same comparison. Both outcomes are recorded in the task report.

See [`bugfixes-35e8a7e/`](bugfixes-35e8a7e/) for the 6-thread noise-floor
baseline and its summary.
