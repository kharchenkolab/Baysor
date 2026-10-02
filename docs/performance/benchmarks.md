# Benchmarks

Run time, memory and segmentation accuracy of Baysor cpp-0.9.0 compared with
cpp-0.8.3. Both versions used the same commands on 77 datasets: 26 real-data
crops from six platforms and 51 simulations with ground truth.

!!! success "At a glance"

    - **About 2× faster** at 6 threads: median wall time 0.46× that of
      cpp-0.8.3 over all 77 datasets (CPU time 0.44×).
    - **Faster large-panel initialisation.** Above 3,000 genes the
      initialisation uses a truncated solver: a 4,963-gene simulation takes
      3.7 s instead of 237 s and needs 164 MiB instead of 901 MiB.
    - **Less memory**: median peak RSS 0.82× of cpp-0.8.3; the 2-million-molecule
      crops need about 1.1 GiB instead of 1.4–2.0 GiB.
    - **Same accuracy**: on the 51 simulations the one-to-one accuracy changes
      by a median of +0.001.
    - **Reproducible multi-threaded runs**: repeated runs give identical output
      on 70 of 70 datasets (cpp-0.8.3: 0 of 70).
    - **One slower case**: the 130k-molecule Xenium Prime 5K crop, whose
      2,672 used genes still go through the dense initialisation, takes
      73 s instead of 66 s at 6 threads.

??? info "Measurement setup"

    | | |
    |---|---|
    | CPU | Intel Xeon E5-2670 @ 2.60 GHz (Sandy Bridge, 2012), 1 socket, 8 cores / 16 threads, AVX but no AVX2 |
    | memory, OS | 60 GB RAM, Ubuntu 26.04 |
    | Baysor | cpp-0.9.0 release candidate (`perf-optimization` @ `1636267`, output-identical to the release) and cpp-0.8.3 |
    | build | both from source with the conda GCC 15.3 toolchain: `./configure.sh --deps=conda`, CMake `Release` (`-O3 -DNDEBUG`, generic x86-64, no `-march=native`) |
    | runs | 2026-10-01; [baysor-benchmarks](https://github.com/VPetukhov/baysor-benchmarks) harness, the dataset's command line plus `--skip-ncv-color --output-style parquet`; threads set with `OMP_NUM_THREADS`; wall time and peak RSS from `/usr/bin/time -v` |
    | replicates | 3 per version at 6 threads (1 for the 2M-molecule crops), 1–2 at 1 thread; the two versions alternate run by run |

    The host was shared with other jobs. During the comparison runs the
    1-minute load average was 5–21 on 16 logical CPUs, so the wall times
    below include contention for both versions. CPU times and alternating
    version comparisons help distinguish computational cost from contention.
    The [thread sweep](#threads) ran one process at a time at a load of 2–13.

## Run time and memory on real data

![Wall time and peak memory against the number of molecules for all real datasets, cpp-0.9.0 filled, cpp-0.8.3 open](img/runtime_vs_molecules-light.svg#only-light)
![Wall time and peak memory against the number of molecules for all real datasets, cpp-0.9.0 filled, cpp-0.8.3 open](img/runtime_vs_molecules-dark.svg#only-dark)

**Figure 1.** Wall time (left) and peak resident memory (right) at 6 threads
against the number of molecules, for the 26 real-data crops (log–log). Filled
marks are cpp-0.9.0, open marks cpp-0.8.3, joined per dataset. Run time grows
roughly linearly with the number of molecules. The labelled crops with
thousands of genes sit above the trend.

| dataset | molecules | genes | wall, 6 threads | peak RSS, 6 threads | wall, 1 thread |
|---|---:|---:|---:|---:|---:|
| Xenium pancreas | 130,001 | 376 | 8.3 s (17.9 s) | 187 MiB (251 MiB) | 36 s (47 s) |
| Xenium lung cancer | 550,000 | 377 | 50 s (90 s) | 393 MiB (533 MiB) | |
| Xenium breast cancer | 1,991,578 | 321 | 138 s (273 s) | 1.06 GiB (1.86 GiB) | |
| Xenium pancreas | 1,999,991 | 377 | 194 s (342 s) | 1.06 GiB (1.81 GiB) | |
| Xenium Prime 5K ovarian | 130,000 | 4,603 | 73 s (66 s) | 353 MiB (418 MiB) | 83 s (94 s) |
| Xenium Prime 5K ovarian | 1,999,997 | 5,088 | 242 s (10.2 min) | 1.06 GiB (1.95 GiB) | |
| CosMx lung cancer | 1,967,437 | 960 | 243 s (388 s) | 1.07 GiB (1.42 GiB) | |
| CosMx WTx colon | 149,771 | 17,533 | 35 s (91 s) | 1.01 GiB (1.71 GiB) | |
| MERFISH ileum (3D) | 819,665 | 241 | 67 s (127 s) | 534 MiB (876 MiB) | |
| ISS hippocampus | 82,077 | 84 | 8.9 s (14.3 s) | 158 MiB (187 MiB) | 29 s (37 s) |
| osmFISH cortex | 149,997 | 35 | 8.1 s (18.1 s) | 169 MiB (233 MiB) | |
| STARmap cortex (3D) | 149,994 | 1,020 | 10.8 s (22 s) | 218 MiB (268 MiB) | |

cpp-0.9.0 first, cpp-0.8.3 in brackets; medians over the replicates. "Genes"
is the panel size of the crop. All datasets are described in the benchmark
repository's [DATASETS.md](https://github.com/VPetukhov/baysor-benchmarks/blob/main/DATASETS.md).

??? note "Second sample: the five largest real crops, one process at a time"

    The same command lines, re-run on 2026-10-02 one process at a time, the
    two versions alternating. Other jobs kept the host busy (1-minute load
    6–22), so the wall times are not lower than above. The wall-time ratios
    cpp-0.9.0 / cpp-0.8.3 are similar: 0.32–0.56 here, 0.39–0.63 above.

    | dataset | wall, 6 threads | CPU time | peak RSS | load at start |
    |---|---:|---:|---:|---:|
    | MERFISH ileum (820k) | 70 s (137 s) | 352 s (10.3 min) | 544 MiB (885 MiB) | 11, 8 |
    | Xenium breast cancer (2M) | 150 s (314 s) | 12.1 min (22.2 min) | 1.06 GiB (1.83 GiB) | 6, 8 |
    | Xenium pancreas (2M) | 191 s (449 s) | 15.1 min (25.6 min) | 1.07 GiB (1.74 GiB) | 9, 11 |
    | CosMx lung cancer (2M) | 374 s (11.0 min) | 25.7 min (33.1 min) | 1.08 GiB (1.46 GiB) | 15, 22 |
    | Xenium Prime 5K ovarian (2M) | 377 s (19.8 min) | 27.0 min (40.9 min) | 1.06 GiB (1.83 GiB) | 15, 17 |

    cpp-0.9.0 first, cpp-0.8.3 in brackets.

## Threads

![Wall time and CPU time against the number of threads for two datasets](img/threads-light.svg#only-light)
![Wall time and CPU time against the number of threads for two datasets](img/threads-dark.svg#only-dark)

**Figure 2.** Wall time (left) and total CPU time (right) for 1–16 threads on
the 130k-molecule Xenium pancreas crop and the 820k-molecule 3D MERFISH crop
(log–log; one run per point, one process at a time). The dotted line is ideal
scaling from the 1-thread time of cpp-0.9.0. The grey vertical line marks the
8 physical cores; 16 threads use hyper-threading.

cpp-0.9.0 keeps scaling up to the 8 physical cores: 8 threads are 6.0×
faster than 1 thread on the Xenium crop (42.8 s → 7.2 s) and 5.8× on the
MERFISH crop (346 s → 60 s), while its total CPU time stays nearly flat
(+7 % and +11 %). cpp-0.8.3 reaches 3.5× and 3.7× and spends far more CPU
time doing it (54 s → 101 s on the Xenium crop), largely OpenMP threads
spin-waiting (see [Profiling](profiling.md#threads)). At 8 threads cpp-0.9.0 is 2.1× (Xenium) and 2.0× (MERFISH)
faster than cpp-0.8.3; at 1 thread it is 1.27× faster on both. Hyper-threads
(16 threads) add about 20 %.

By default cpp-0.9.0 uses one thread per physical core (8 on this host;
cpp-0.8.3's OpenMP default was all 16 logical CPUs). Set the count with
`-t/--threads`, the `threads` config key or `OMP_NUM_THREADS`; see
[Threading](../run.md#threading).

## Gene panel size

![Wall time and peak memory against the number of genes in the panel for three simulated families](img/time_vs_genes-light.svg#only-light)
![Wall time and peak memory against the number of genes in the panel for three simulated families](img/time_vs_genes-dark.svg#only-dark)

**Figure 3.** Wall time (left) and peak memory (right) at 6 threads against
the gene-panel size, for three families of simulations that differ only in the
number of genes (64k–116k molecules each). Up to about 1,000 genes the panel
size hardly matters. cpp-0.8.3 then grows steeply (dense ICA initialisation,
cubic in the number of genes); cpp-0.9.0 switches to a truncated solver above
3,000 genes. Real panels between 1,000 and 3,000 genes still take the dense
path, which is the one case where cpp-0.9.0 is not faster (see
[Profiling](profiling.md#current-bottlenecks)).

## Accuracy and reproducibility

![Scatter of one-to-one accuracy, cpp-0.9.0 against cpp-0.8.3, on 51 simulated datasets](img/accuracy_sim-light.svg#only-light){ width="420" }
![Scatter of one-to-one accuracy, cpp-0.9.0 against cpp-0.8.3, on 51 simulated datasets](img/accuracy_sim-dark.svg#only-dark){ width="420" }

**Figure 4.** One-to-one accuracy against the simulated ground truth (fraction
of true cells matched to exactly one segmented cell), cpp-0.9.0 against
cpp-0.8.3, mean of 3 runs at 6 threads, for 51 simulations. The dotted line is
equal accuracy. The median change is +0.001 (range −0.008 to +0.019); 11
datasets gain more than 0.005 and 4 lose more than 0.005. The tissue-like
simulations (squares) are hard by design and score 0.40–0.73 with both
versions.

![Cell counts of repeated runs on the ISS dataset for both versions at 6 threads and 1 thread](img/determinism_iss-light.svg#only-light)
![Cell counts of repeated runs on the ISS dataset for both versions at 6 threads and 1 thread](img/determinism_iss-dark.svg#only-dark)

**Figure 5.** Cells found on the ISS mouse hippocampus crop by repeated runs.
cpp-0.8.3 gives a different segmentation on every multi-threaded run, and its
6-thread results drift away from the 1-thread result (10,634–10,736 vs 10,051
cells): all threads re-seeded their random streams identically. cpp-0.9.0 keys
the random streams by iteration and work chunk, so multi-threaded outputs are
reproducible across multi-threaded counts (three identical runs, 9,996 cells)
and stay close to the 1-thread result (10,037 cells). The single-threaded
result can differ; reproducibility does not mean identical output between
versions.

On real data, where there is no ground truth, the two versions agree with
each other as well as two runs of cpp-0.8.3 agree among themselves (e.g. the
550k-molecule Xenium lung crop: ARI 0.860 between versions, 0.859 between
cpp-0.8.3 runs), and cell counts differ by at most 2 % on every dataset
except ISS.

## What the segmentations look like

![Four 100 µm crops of cpp-0.9.0 segmentations: Xenium, MERFISH, CosMx and ISS](img/segmentation_examples-light.png#only-light)
![Four 100 µm crops of cpp-0.9.0 segmentations: Xenium, MERFISH, CosMx and ISS](img/segmentation_examples-dark.png#only-dark)

**Figure 6.** cpp-0.9.0 segmentations of four platforms. Each panel is the
densest 100 × 100 µm window of a benchmark crop (150 × 150 µm for the sparse
ISS data). Molecules are coloured by cell (neighbouring cells get different
colours), noise molecules are grey, and the lines are Baysor's cell polygons.
MERFISH is 3-D: the molecules of one z-plane are shown with the outlines of
the cells that have at least 8 molecules in it; neighbouring 3-D cells can
overlap in this projection.

![UMAP of cells from two Xenium datasets, coloured by cell type](img/umap-light.png#only-light)
![UMAP of cells from two Xenium datasets, coloured by cell type](img/umap-dark.png#only-dark)

**Figure 7.** UMAP of the cells segmented by cpp-0.9.0 on two Xenium crops,
computed from the cell × gene count matrix (normalised, log-transformed,
30 principal components). Colours and numbers are the cell types used by the
cellAdmix audit below (transferred from a reference segmentation of the same
crop); grey cells could not be typed.

![cellAdmix total admixture rate per dataset and the top cell-type pairs on the Xenium lung crop](img/celladmix-light.svg#only-light)
![cellAdmix total admixture rate per dataset and the top cell-type pairs on the Xenium lung crop](img/celladmix-dark.svg#only-dark)

**Figure 8.** Admixture measured with
[cellAdmix](https://github.com/kharchenkolab/cellAdmix-core): the share of
molecules in a cell that leaked in from a neighbouring cell of another type
(lower is cleaner). Top: total rate per real dataset with at least 2,000
cells, cpp-0.9.0 (filled) and cpp-0.8.3 (open), mean of the replicates.
Bottom: the eight largest source → target cell-type pairs on the 550k Xenium
lung crop; the type numbers are those of Figure 7 (left). The two versions are
close on every dataset. On all five 550k-molecule crops built for this audit
cpp-0.9.0 is slightly cleaner, e.g. 5.3 % vs 5.4 % on the Xenium lung crop;
the 130k Xenium Prime 5K crop goes the other way (3.7 % vs 3.4 %).

## Reproduce

Figures and tables are regenerated from the benchmark results by
[`docs_figures/make_figures.py`](https://github.com/VPetukhov/baysor-benchmarks/tree/main/docs_figures)
in baysor-benchmarks; `docs_figures/generated/tables.md` there lists the
source file of each number. To benchmark your own build, see
[Development › Benchmarks](../development.md#benchmarks).
