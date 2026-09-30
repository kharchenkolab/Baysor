# Xenium workflow

This page describes the recommended workflow for 10x Xenium data.

## Recommended input

Use the Xenium manifest as the main input:

```bash
baysor run -c configs/xenium.toml -o out data/experiment.xenium :cell_id
```

Why:

- Baysor resolves the transcript table automatically
- the Xenium column mapping and filtering are captured in
  [configs/xenium.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/xenium.toml)
- Xenium `transcript_id` is preserved in `legacy` output

Passing `transcripts.parquet` directly still works, but `experiment.xenium` is
the preferred Xenium-aware entrypoint.

## Common prior modes

### Transcript-native prior (recommended)

```bash
baysor run -c configs/xenium.toml -o out data/experiment.xenium :cell_id
```

### Cell boundary prior

```bash
baysor run -c configs/xenium.toml -o out \
  data/experiment.xenium data/cell_boundaries.parquet
```

### Nucleus boundary prior

```bash
baysor run -c configs/xenium.toml -o out \
  data/experiment.xenium data/nucleus_boundaries.parquet
```

Boundary tables are long-format vertex tables
([Prior segmentation](priors.md#boundary-tables-csv-parquet)).

## Recommended output style

For Xenium runs that may be handed off to Xenium Explorer, use `legacy`
output (the default):

```bash
baysor run -c configs/xenium.toml \
  --output-style legacy \
  -o out \
  data/experiment.xenium :cell_id
```

For Xenium-origin inputs this automatically produces Ranger-friendly
`segmentation.csv` and `segmentation_polygons_2d.json` (see
[Output files](output_files.md#legacy-bundle)).

## Very large / 5K panel runs

For very large Xenium runs, particularly high-gene-panel datasets such as 5K
panels, prefer Louvain clustering with about 10 final coarse clusters:

```bash
baysor run \
  -c configs/xenium.toml \
  --cluster-method louvain \
  --n-clusters 10 \
  -o out \
  data/experiment.xenium :cell_id
```

The Louvain path clusters NCV basis anchors and transfers labels back to all
molecules, which is usually a better large-run starting point than the legacy
MRF clustering prior. Keep `--n-clusters` near `10` unless the diagnostic
report shows clear under- or over-clustering.

## Xenium Explorer handoff

The recommended Explorer path is:

1. run Baysor on the original Xenium bundle in `legacy` mode
2. run `xeniumranger import-segmentation`

The two Baysor files used for the handoff are `segmentation.csv` and
`segmentation_polygons_2d.json`. Run the conversion from the directory that
contains the original Xenium bundle, or provide an absolute bundle path:

```bash
xeniumranger import-segmentation \
  --id baysor_xenium \
  --xenium-bundle data \
  --transcript-assignment out/segmentation.csv \
  --viz-polygons out/segmentation_polygons_2d.json \
  --units microns
```

You can add normal Xenium Ranger execution flags such as `--localcores` and
`--localmem`. This is the preferred route instead of direct Baysor-side Xenium
bundle generation.

## Full runs vs crops

Use full runs for Xenium Ranger handoff. Cropped runs (`--x-min` and friends)
are useful for development, debugging, visualization, and performance
profiling, but they are not the right input to `xeniumranger
import-segmentation`.

## Example

For a full runnable Xenium example, see
[examples/Xenium_pancreas_membrane_377](https://github.com/kharchenkolab/Baysor/tree/HEAD/examples/Xenium_pancreas_membrane_377).
