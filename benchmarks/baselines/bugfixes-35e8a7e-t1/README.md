# Baseline: `bugfixes-35e8a7e-t1` (exact, `identical` flavour)

> **Note.** Branch `cpp-dev-llm` includes the canonical graph edge order fix
> (`60902fe`), which changes 1-thread results (e.g. ISS 10024 -> 10038 cells).
> Builds of `cpp-dev-llm` are therefore not expected to pass `--expect identical`
> against this baseline; recreate the baselines with the `release` suite first.


Official 1-thread reference baseline of the Baysor benchmark suite for the
**Release build of branch `bugfixes` at commit `35e8a7e`**
(`--label bugfixes-35e8a7e`, binary
`/home/vpetukhov/Projects/Baysor/.bench-data/binaries/baysor-bugfixes-35e8a7e`),
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
B=/home/vpetukhov/Projects/Baysor/.bench-data/binaries/baysor-bugfixes-35e8a7e

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

## Self-check (outcomes)

| check | run | expect | outcome |
|---|---|---|---|
| fresh 1-thread, 1-replicate run of 10 quick datasets (5 sim + 5 real) vs this baseline | `chk-t1a` | `identical` **PASS** | ✅ exit 0, 70 passed / 0 failed, all 10 `assignment_sha256` bitwise equal |
| `--scale-factor 0.9` run of the same subset | `selfcheck-t1-scale09` | `identical` **FAIL** | ✅ exit 1 (assignments differ on the degraded run) |

**Known binary issue discovered while validating this baseline** (details in
[`../../harness/README.md`](../../harness/README.md) → "Determinism
findings"): the first two self-check attempts used run-ids of 18–19
characters and failed reproducibly on `iss_mouse_hippocampus_quick` and
`osmfish_somatosensory_quick` — at 1 thread the binary's `-o` output-path
length participates in its (layout-sensitive) behaviour, and both datasets
flip deterministically once the run-id reaches 18 characters. Evidence: 7
runs with run-ids of 11–17 characters reproduce this baseline bitwise, 4
runs with 18–19-character run-ids produce a second stable outcome. Use a
run-id of ≤ 17 characters for `--expect identical` (the failing attempts
`selfcheck-t1-fresh{,2,3}` are kept under `$BAYSOR_BENCH_DATA/runs/` as
evidence; `compare.py` now warns about long run-ids).

See [`bugfixes-35e8a7e/`](bugfixes-35e8a7e/) for the 6-thread noise-floor
baseline and its summary.
