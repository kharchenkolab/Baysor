# Xenium workflow

Download the preset, then run on the original Xenium output bundle:

```bash
curl -fL -o xenium.toml \
  https://raw.githubusercontent.com/kharchenkolab/Baysor/cpp-0.9.0/configs/xenium.toml
baysor run -c xenium.toml -o out --threads 8 data/experiment.xenium :cell_id
```

`data/` is the directory containing `experiment.xenium` and its transcript
table. Baysor locates the table automatically; `:cell_id` uses its existing
cell assignments as the prior. The preset maps the Xenium columns, filters
control genes, sets `-m` to `50`, uses prior confidence `0.5` and treats
`UNASSIGNED` as unassigned. It does not enable quality-value filtering; add
`--min-qv` if needed. See [Configuration](configuration.md#protocol-presets).

Passing `transcripts.parquet` directly also works, but the manifest is the
recommended entrypoint. Add `-p` for [HTML diagnostics](outputs.md).

## Alternative priors

To use cell boundaries instead of molecule labels:

```bash
baysor run -c xenium.toml -o out data/experiment.xenium data/cell_boundaries.parquet
```

Use `data/nucleus_boundaries.parquet` instead for a nucleus prior. These are
[vertex tables](priors.md#boundary-tables-csv-parquet), assigned by x/y. For
no prior, omit the second input and provide `-s` / `--scale` in microns. See
[Prior segmentation](priors.md) for confidence and scale estimation.

## Very large / 5K panel runs

For large, high-gene-panel datasets, start with Louvain and about `10` coarse
molecule clusters:

```bash
baysor run -c xenium.toml --cluster-method louvain --n-clusters 10 \
  -o out --threads 8 data/experiment.xenium :cell_id
```

This uses a neighborhood-composition graph instead of the default MRF
clustering prior. Inspect the diagnostic report before adjusting the cluster
count; it describes coarse cell types, not the number of segmented cells.

## Xenium Explorer handoff

Use the full dataset and the default `legacy` output style. It writes
`segmentation.csv` with transcript IDs and
`segmentation_polygons_2d.json` with matching cell IDs. Cells whose boundary
estimation fails get fallback rectangles instead of being omitted.

Choose polygon format for the importing Xenium Ranger version:

| Xenium Ranger | Baysor setting |
| --- | --- |
| 4.0 and later | Default `FeatureCollection`; no extra flag. |
| 3.1 and earlier | Add `--polygon-format GeometryCollectionLegacy` to the Baysor command for integer polygon cell IDs. |

Then run:

```bash
xeniumranger import-segmentation \
  --id baysor_xenium \
  --xenium-bundle data \
  --transcript-assignment out/segmentation.csv \
  --viz-polygons out/segmentation_polygons_2d.json \
  --units microns
```

`data` must point to the original bundle. Do not use cropped Baysor runs for
this import. For Python / R analysis without Ranger, consider
[`--output-style parquet`](outputs.md#parquet-output) instead.

For data downloads and a full example, see
[Xenium pancreas](https://github.com/kharchenkolab/Baysor/tree/cpp-0.9.0/examples/Xenium_pancreas_membrane_377).
