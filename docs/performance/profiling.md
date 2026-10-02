# Profiling

Where Baysor cpp-0.9.0 spends its time and memory, how that grows with the
size of the data, and what is still slow. The [benchmarks](benchmarks.md)
compare versions; this page looks inside one run.

!!! success "At a glance"

    - **Run time is linear in the number of molecules**: CPU time grows with
      exponent 1.08 from 0.1 to 10.6 million molecules. A whole Xenium lung
      slide (10.6M molecules, 377 genes) needs 168 CPU-minutes and 5.3 GiB.
    - **Memory: about 0.5 KB per molecule** at slide scale (539 bytes on the
      whole lung slide); peak RSS is 31–48 % lower than before the
      optimisation work on slides of 1M molecules and more.
    - **95 % of the work runs in parallel** at real sizes, so 8 threads can
      give up to 5.9× (Amdahl bound). On small crops a single-threaded step,
      the colour embedding, caps the gain at about 1.4×.
    - **The bottlenecks now are** the BMM E-step (the core of the algorithm),
      molecule clustering, which grows faster than linearly, and gene-rich
      panels.

!!! info "Setup"

    Same host as the [benchmarks](benchmarks.md): Intel Xeon E5-2670
    (Sandy Bridge, 8 cores / 16 threads, no AVX2), 60 GB RAM, Ubuntu 26.04.
    Code: `perf-optimization` @ `20bc45c` (2026-10-01), which has all the
    performance work of cpp-0.9.0 (the later commits fix bugs and simplify
    code), compared with `e45fddc` (2026-09-30, the start of the optimisation
    work). Build: `cmake --preset profiling` = Release
    (`-O3 -DNDEBUG`) plus `-g`, conda GCC 15.3; its machine code is identical
    to the Release binary. Unlike the benchmarks, the profiled runs compute the
    neighbourhood colour embedding (the default; `--skip-ncv-color` is off).
    The host was shared: wall times of the real-size runs were taken at a
    1-minute load of 5–53, so they are pessimistic; CPU time, instruction
    counts and memory do not depend on the load.

## What was measured

| tool | what it gives | data |
|---|---|---|
| Valgrind callgrind | instructions per pipeline phase and function, 1 thread; exact and load-independent | 11 crops of 10k–40k molecules and 35–8,407 genes, cut from 9 real and simulated datasets |
| Valgrind DHAT, heaptrack | heap peak, allocation counts and sites | the crops; 2M-molecule and whole-slide runs |
| gperftools + `/proc` sampling | CPU time per phase and function, peak RSS, parallel share | real-size ladders: a Xenium lung slide (0.1M → 10.6M molecules) and a Xenium Prime 5K slide (0.1M → 8M), 1 and 8 threads; a 3M-molecule CosMx whole-transcriptome slide |
| native timing | wall and CPU time, 1–16 threads | the crops, 3 runs each |

The suite and its HTML report generator live in
[`profiling/`](https://github.com/VPetukhov/baysor-benchmarks/tree/main/profiling)
of baysor-benchmarks (`profile.py`, `scaling.py`, `summarize.py`,
`report_html.py`).

## Where the time goes

![Stacked bars: instructions per phase on 20k-molecule crops and CPU share per phase at real sizes](img/phases-light.svg#only-light)
![Stacked bars: instructions per phase on 20k-molecule crops and CPU share per phase at real sizes](img/phases-dark.svg#only-dark)

**Figure 1.** Left: instructions per phase on 20k-molecule crops (callgrind,
1 thread). Right: share of the CPU time per phase at real sizes (gperftools,
8 threads). "BMM iterations" is the segmentation itself; "molecule
clustering" is the initial assignment of molecules to cell types (the MRF
clustering or, for whole-transcriptome panels, the neighbourhood-graph
clustering); "NCV colours" is the colour embedding used by the plots and the
`ncv_color` output column.

- On **small crops** the colour embedding dominates: 66–84 % of the
  instructions on every panel below 1,000 genes, mostly umappp's
  single-threaded layout optimisation. It costs a fixed amount (at most
  20,000 anchors), so it fades at real sizes: 8 % of the CPU time on the whole
  lung slide. Skip it with `--skip-ncv-color` if you do not need the colours.
- On **real slides** (1M molecules and more) the BMM iterations take 37–74 %
  of the CPU time and molecule clustering 13–43 %; on the 18,935-gene CosMx slide molecule
  clustering takes 79 %.
- The **dense ICA** behind the 2,793-gene crop's huge clustering bar is the
  1,000–3,000-gene case of [Benchmarks › Gene panel size](benchmarks.md#gene-panel-size).

## Scaling with the size of the data

![CPU time, wall time and peak memory against the number of molecules for the real-size ladders](img/scaling-light.svg#only-light)
![CPU time, wall time and peak memory against the number of molecules for the real-size ladders](img/scaling-dark.svg#only-dark)

**Figure 2.** Real-size ladders: random subsets of a Xenium lung slide
(blue) and a Xenium Prime 5K slide (orange), up to the whole slide, and a
CosMx whole-transcriptome slide (green). Left: total CPU time at 8 threads
(dotted: linear). Middle: wall time at 8 threads (solid) and 1 thread (dotted;
1-thread runs stop at 2M molecules). Right: peak memory at 8 threads. Dashed
open marks are the code before the optimisation (`e45fddc`). Wall times are
load-sensitive: the 8M lung rung ran at a load of 20, which explains its jump.

| | lung 1M | lung 2M | lung 10.6M (slide) | Prime 5K 1M | Prime 5K 8M | CosMx WTx 3M |
|---|---:|---:|---:|---:|---:|---:|
| genes | 377 | 377 | 377 | 4,358 | 5,078 | 18,935 |
| CPU time, 8 threads | 9.7 min | 24.3 min | 168 min | 13.7 min | 161 min | 98 min |
| CPU time before | 14.0 min | 30.1 min | 195 min | 23.1 min | 221 min | 257 min |
| wall time, 8 threads | 122 s | 297 s | 63 min | 152 s | 31 min | 19.5 min |
| peak RSS | 676 MiB | 1.15 GiB | 5.32 GiB | 843 MiB | 4.24 GiB | 2.16 GiB |
| peak RSS before | 1.13 GiB | 1.95 GiB | 8.17 GiB | 1.47 GiB | 6.87 GiB | 3.11 GiB |
| RSS per molecule | 709 B | 620 B | 539 B | 888 B | 569 B | 779 B |

??? note "Scaling exponents per phase"

    Fitted exponent *b* in value ∝ molecules<sup>*b*</sup> for the CPU time
    of the whole run and of each phase, and for the peak RSS (gperftools,
    8 threads; lung 0.1M–10.6M, Prime 5K 0.1M–8M). 1 = linear.

    | | lung | Prime 5K |
    |---|---:|---:|
    | whole run, CPU time | 1.08 | 1.07 |
    | whole run, peak RSS | 0.73 | 0.59 |
    | BMM iterations | 1.15 | 1.10 |
    | molecule clustering | **1.32** | **1.23** |
    | polygons | 1.24 | 1.38 |
    | molecule graph | 1.25 | 1.02 |
    | confidence | 1.04 | 1.10 |
    | NCV colours (capped) | 0.71 | 0.86 |

    Molecule clustering is the only large phase that is clearly super-linear:
    its iteration count grows with the slide (181 → 592 iterations on the lung
    ladder), so its share rises from 6 % (0.1M) to 23 % (whole slide) of the
    lung CPU time.

## Threads

At real sizes 95 % of the CPU time is in parallel code (lung 2M: serial
share 5.1 %, was 12.7 % before), and the thread pool hardly waits (≤ 0.4 % of
the CPU time at 8 threads on rungs of 1M molecules and more). That bounds the
speed-up at 5.9× on 8 threads and 9.1× on 16. The serial remainder is mostly
graph construction, confidence estimation and polygons.

Small crops behave differently: 62 % of the instructions of a 20k-molecule
crop run serially, most of them in the colour embedding, so it speeds up only
1.38× on 8 threads (18.9 s → 13.7 s). With `--skip-ncv-color`, as in the
[benchmark thread sweep](benchmarks.md#threads), the same kind of data scales
much further. Running more threads than physical cores gains little.

??? note "Crop thread series (with the colour embedding)"

    Minimum wall time and median CPU time of 3 runs.

    | threads | 1 | 2 | 4 | 8 | 16 |
    |---|---:|---:|---:|---:|---:|
    | Xenium pancreas 20k, wall | 18.9 s | 15.9 s | 14.4 s | 13.7 s | 13.8 s |
    | Xenium pancreas 20k, CPU | 21.3 s | 19.2 s | 20.0 s | 21.2 s | 24.3 s |
    | Xenium Prime 5K 20k, wall | 14.5 s | 12.5 s | 12.1 s | 10.2 s | 11.0 s |
    | Xenium Prime 5K 20k, CPU | 14.5 s | 16.5 s | 16.7 s | 16.6 s | 18.0 s |

    Before the optimisation the 20k pancreas crop used 97 s of CPU time on
    16 threads (OpenMP spin-waiting); now 24 s.

## Current bottlenecks

Ranked by their share of the CPU time at real sizes (whole lung slide and
Prime 5K 8M, 8 threads).

1. **BMM E-step** — about 35 % of the CPU on the whole lung slide (the E-step
   chunk, `CategoricalSmoothed::pdf` and `exp`). Linear and 99 % parallel:
   this is the algorithm's core work, not overhead.
2. **Molecule clustering** — 23 % (lung) and 42 % (Prime 5K 8M) of the CPU,
   super-linear (exponent 1.32), because the MRF clustering needs more
   iterations on larger slides.
3. **BMM bookkeeping** (splitting and grouping components, hash-map updates)
   — about 22 % of the CPU on the whole slide, mildly super-linear
   (exponents 1.2–1.4).
4. **NCV colour embedding** on small data — single-threaded, 61 % of the
   instructions on a 20k crop; negligible on slides. Use `--skip-ncv-color`
   when the colours are not needed.
5. **Gene-rich panels** — the neighbourhood k-NN with k = genes / 10 takes
   65 % of the CPU on the 18,935-gene CosMx slide, and panels of
   1,000–3,000 genes still run the dense O(genes³) ICA (87 % of the
   instructions on the 2,793-gene crop, single-threaded).

Memory at the whole-slide peak (5.06 GiB heap) is spread over the molecule
adjacency list (787 MiB), the assignment history (728 MiB), the MRF clustering
state (383 MiB) and three copies of the molecule positions (3 × 170 MiB).

## Reproduce

The full report — every table, the per-function and per-line hotspots, and
the method — is produced by the
[profiling suite](https://github.com/VPetukhov/baysor-benchmarks/tree/main/profiling)
(`profile.py` for the crops, `scaling.py` for the ladders, `report_html.py`
for the HTML report). The figures and tables on this page are regenerated
from its summaries by
[`docs_figures/make_figures.py`](https://github.com/VPetukhov/baysor-benchmarks/tree/main/docs_figures).
