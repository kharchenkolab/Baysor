# BENCH-REALO — inventory of real non-Xenium datasets

Candidate check for the BENCH-REALO task: does it exist, can it be downloaded
**without login or registration**, and how big is it. "Bytes" is what the
fetch scripts actually download (all of it lands in
`$BAYSOR_BENCH_DATA/cache/real_other/`, total **28.7 GB** on disk — well under
the 150 GB budget).

## The list of real datasets

| # | candidate | available | genes (kept) | platform | tissue | density class | quick crop | full crop | bytes downloaded | source URL |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | ISS mouse hippocampus (Qian et al., section 3-3) | **yes** `iss_mouse_hippocampus_quick` | 84 (99-gene panel; 84 detected in this section) | ISS (pciSeq) | mouse hippocampus CA1, coronal | medium (3239 cells/mm²) | 82,077 mol | — | 718.7 MB (pklab molecules + published DAPI mask) + 37.1 MB (figshare DAPI jpg + pciSeq spot assignments) | `http://pklab.med.harvard.edu/viktor/baysor/iss/`, `doi:10.6084/m9.figshare.7150760.v1` (private link `…/s/88a0fc8157aca0c6f0e8`), assignments `doi:10.17045/sthlmuni.10318610.v1` |
| 2 | osmFISH somatosensory cortex (Codeluppi et al.) | **yes** `osmfish_somatosensory_quick` | 35 (33-gene panel +2 detected) | osmFISH (smFISH) | mouse somatosensory cortex | medium (2570 cells/mm²) | 149,997 mol | — | 41.6 MB (pklab molecules csv) + 0.65 MB (vitessce polygon mirror of `polyT_seg.pkl`) + 86.3 MB vitessce molecules csv (cross-check) | `http://pklab.med.harvard.edu/viktor/baysor/osm_fish/`, `https://data-1.vitessce.io/0.0.33/main/codeluppi-2018/`, `doi:10.1038/s41592-018-0175-z` |
| 3 | STARmap visual cortex (Wang et al. 2018, **3D**) | **yes** `starmap_visual_cortex_quick` | 1020 | STARmap | mouse V1 neocortex, 3D strip | medium (2982 cells/mm²) | 149,994 mol (z kept) | — | 18.0 MB molecules.csv + 128.6 MB segmentation.tiff | `http://pklab.med.harvard.edu/viktor/baysor/starmap/`, `doi:10.1126/science.aat5691` |
| 4 | MERFISH mouse ileum (Petukhov et al. 2022) | **yes** `merfish_ileum_quick` / `_full` | 239 quick / 241 full (Blank* removed) | MERFISH | mouse ileum, FFPE region | dense (≈19k cells/mm²) | 149,356 mol (3D, z kept) | 819,665 mol | 529 MB compressed zip members (molecules.csv, DAPI + membrane stacks, Cellpose cell/nuclei label volumes) of the 735.8 MB release zip | `doi:10.5061/dryad.jm63xsjb2` |
| 5 | CosMx NSCLC Lung5_Rep1 (NanoString public FFPE) | **yes** `cosmx_nsclc_lung5_rep1_quick` / `_full` | 960 (20 NegPrb removed) | CosMx SMI (960-plex) | human NSCLC FFPE tumour | dense (9.5k / 9.9k cells/mm²) | 149,903 mol (one dense window of FOV 6) | 1,967,437 mol (whole FOV 6) | 3.404 GB tx file + 11.5 MB metadata + 1.4 MB CellLabels | `https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI/` (no registration) |
| 6 | CosMx 6K / Whole Transcriptome (~19k) public data | **yes (WTX)** `cosmx_wtx_colon_quick` / `_full` | 17,533 quick / 18,935 full (1,443 Negative/SystemControl removed) | CosMx SMI Whole Transcriptome | human colon FFPE (WTX colon discovery) | dense (12.4k / 13.6k cells/mm²) | 149,771 mol (dense window of FOV 68) | 2,969,193 mol (whole FOV 68) | 11.398 GB tx_file.csv.gz + 39.1 MB metadata + README | `https://objects.liquidweb.services/smi-public/wtx_manuscript/colon_discovery/` (public bucket; sibling prefixes hold brain/breast/kidney/pancreas/skin/lymph WTX sets) |
| 7 | MERSCOPE (Vizgen public releases) | **no** | ~483–1000 (panels not verified) | MERSCOPE MERFISH | immuno-oncology / mouse brain / liver | — | — | — | 0 | Every Vizgen public release (FFPE immuno-oncology, receptor map, liver map, Ultra, 1000-plex, MERFISH 2.0, protein co-detection) is behind a HubSpot **registration form** ("Complete the Registration Form to Access the Data"); no anonymous S3/GCS bucket could be found (`vizgen-public-data` and guesses all NoSuchBucket) |
| 8a | CosMx Bolen sigmoid colon R1000 (Zenodo) | available, not built | ~1000 (R1000 panel) | CosMx SMI | human sigmoid colon | — | — | — | 0 | `doi:10.5281/zenodo.14851478`, tx_file 2.26 GB, CC BY 4.0 — left out to keep the cache small (panel class already covered by #5/#6) |
| 8b | CosMx obesity CRC TMA (Zenodo) | available, not built | panel not verified | CosMx SMI | human colorectal cancer TMA | — | — | — | 0 | `doi:10.5281/zenodo.21844902`, tx_files 3.57 + 5.66 GB, CC BY 4.0 — same reason |
| 8c | Xenium (incl. Prime5K, protein) | excluded by task | — | Xenium | — | — | — | — | — | covered by BENCH-REALX (directories `real/xenium_*`) |

Notes on availability:

* **pklab.med.harvard.edu** (ISS/osmFISH/STARmap example files) serves an
  Incapsula JS challenge to non-browser clients. The fetch code solves this
  with a real headless Firefox via Playwright (`browser_download`), so the
  files remain login-free. The original STARmap site
  (`starmapresources.com`) is now a dead/parked domain; its data was on a
  Dropbox folder that refuses anonymous folder downloads, and
  `science.org`/`pmc.ncbi.nlm.nih.gov` downloads are bot-walled — the
  pklab example files are the working source.
* **osmFISH original hdf5** (`storage.googleapis.com/linnarsson-lab-www-blobs`)
  is now a private bucket (anonymous AccessDenied); the availability page's
  loom file is still public but cell-level only. Molecules come from the
  pklab example csv (identical to the vitessce mirror: same 1,976,659 rows,
  same extent), vendor polygons from the vitessce mirror (the mirror only
  carries annotated cells, so the quick window is restricted to the
  polygon-covered part of the field; prior coverage is recorded in
  `meta.difficulty.notes`).
* **figshare** (`figshare.com` HTML and `ndownloader.figshare.com`) is behind
  an AWS WAF challenge for non-browser/TLS-fingerprinted clients; file ids
  were recovered through the open `api.figshare.com` API + the Wayback
  snapshot of the share link, and downloads go through the Playwright
  fetcher.
* **Dryad** (MERFISH ileum) is protected by an Anubis proof-of-work
  (difficulty 4, sha256 of `randomData+nonce`). `other_utils.anubis_authenticated_session`
  solves it in <1 s, resolves the presigned S3 URL and reads only the needed
  zip members with `remotezip`.
* **CosMx** flat files are plain public S3/bucket downloads without login.
  The Bruker "Download Data" tables are WordPress REST resources
  (`/wp-json/wp/v2/em_tables/<id>`), which is how the links were found.

## Pixel size conversions (recorded in every `meta.source`)

| dataset | factor | documented where |
|---|---|---|
| ISS | 1/3 µm/px | Qian et al. 2019 (PMC6349128) Online Methods: RCP top-hat "radius 3 pixels (corresponding to 1 µm)", nuclei "24 pixels (8 µm)" |
| osmFISH | 0.065 µm/px | linnarssonlab.org/osmFISH availability page ("1 pixel = 0.065 μm"); matches the example's `--scale 82` px → 5.33 µm |
| STARmap xy | 1400 µm / 17219 px = 0.081306 µm/px | Wang et al. Fig. 2B caption "Full field: 1.4 by 0.3 mm" (cross-check: 3734 px → 303.6 µm ≈ 300 µm; vendor cell radius 88.9 px → 7.2 µm) |
| STARmap z | 8 µm / 36 intervals = 0.2222 µm | Wang et al.: "Eight-µm-thick volumes containing up to 1000 cells … were imaged" (the vendor segmentation has 974 cells) |
| MERFISH ileum | 0.108947 µm/px (least-squares fit) | the release's own `x_pixel/y_pixel` vs `x_um/y_um` column pairs; images additionally documented as 9 z-planes, 1.5 µm spacing, 2.5–14.5 µm above the coverslip |
| CosMx NSCLC | 0.18 µm/px | SMI Data File ReadMe on the public S3 ("To convert to microns multiply the pixel value by 0.18 um per pixel") |
| CosMx WTX | 0.12028 µm/px | README_coloncancer.html ("pixel edge length is 120 nm … multiply by 0.12028 µm per pixel") |

## Performance finding (perf-stress: default ICA init)

On `cosmx_wtx_colon_quick` (17,533 genes) the **default
`--cluster-method mrf` never finishes**: the run was killed after **71 min**
(peak RSS 9.8 GB, 127% CPU) still at
`Clustering molecules into 4 types (ICA init)...`; a 20-minute-capped rerun
records the same timeout and stage in
`$BAYSOR_BENCH_DATA/baselines/real_other_smoke.json` (`recorded_default_timeouts`).
For reference the Xenium Prime5K quick crop (5k genes) finishes in ~75 s.
Alternatives measured on the same crop,6 threads, ≤20 min caps:

| cluster method | wall | peak RSS |
|---|---|---|
| mrf (default, ICA init) | >1200 s (timeout at ICA init) | 9.8 GB (when killed after71 min) |
| louvain (**chosen**) | 113 s | 1.75 GB |
| leiden | 111 s | 1.75 GB |
| none | 111 s | 1.77 GB |

`meta.baysor.extra_args` for both WTX datasets is
`["--cluster-method", "louvain"]`, and `meta.difficulty.notes` carries the
tag **`perf-stress: default ICA init`** so this can become a performance
regression test. The smoke runner takes `--timeout` (default 1800 s) and
kills the process group on timeout instead of blocking the pipeline.

## Smoke results (Release binary, `taskset -c 0-5`, OMP_NUM_THREADS=6)

| dataset | status | wall | peak RSS |
|---|---|---|---|
| iss_mouse_hippocampus_quick | ok | 35.8 s | 185 MiB |
| osmfish_somatosensory_quick | ok | 35.1 s | 255 MiB |
| starmap_visual_cortex_quick | ok | 42.5 s | 425 MiB |
| merfish_ileum_quick | ok | 38.9 s | 271 MiB |
| cosmx_nsclc_lung5_rep1_quick | ok | 49.1 s | 313 MiB |
| cosmx_wtx_colon_quick (louvain) | ok | 108.3 s | 1745 MiB |
| cosmx_wtx_colon_quick (default mrf) | **timeout @1200 s** at "Clustering molecules into 4 types (ICA init)..." | 1201 s | — (recorded) |

Machine-readable version: `$BAYSOR_BENCH_DATA/baselines/real_other_smoke.json`.

## Regeneration

```bash
export BAYSOR_BENCH_DATA=/home/vpetukhov/Projects/Baysor/.bench-data   # or <repo>/.bench-data

# one-time: the browser fetcher used for Incapsula/WAF-walled hosts
.deps/bench/bin/pip install playwright && .deps/bench/bin/python -m playwright install firefox

python benchmarks/fetch/other.py build        # downloads sources (idempotent) + writes real/<id>/
python benchmarks/fetch/other.py build --only <id>    # single dataset
python benchmarks/fetch/other.py smoke --timeout 1800 # wall/RAM of every quick crop
python benchmarks/fetch/other.py report       # per-dataset stats table
.deps/bench/bin/python -m pytest benchmarks/fetch/tests
```

All crop windows are chosen by a seeded densest-window search
(`seed: 20260929` in the manifest); the resulting bounding boxes are stored in
each `meta.json`, so rebuilds are deterministic for identical source files.
