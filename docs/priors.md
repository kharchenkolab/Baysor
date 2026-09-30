# Prior segmentation

Baysor can take an optional prior segmentation into account — for example a
DAPI-based nuclear segmentation or published cell boundaries. The prior is
passed as the second positional argument:

```bash
baysor run [OPTIONS] coordinates prior_segmentation
```

## Prior types

### Transcript-native labels: `:column_name`

If prior assignments already exist per molecule row, pass the column name
with a `:` prefix:

```bash
baysor run -m 50 -c configs/xenium.toml -o out data/transcripts.parquet :cell_id
```

This is the fastest prior mode and the recommended one for Xenium. Values
equal to `--unassigned-prior-label` (default `0`) count as unassigned.

### Image masks: TIFF

A labeled (integer) or binary TIFF mask where pixel values identify segments:

```bash
baysor run -m 30 -c configs/iss.toml -o out molecules.csv dapi_mask.tif
```

Molecules are assigned to the segment covering their x/y position. Both
binary masks (single foreground component per cell) and integer-labeled masks
are supported.

### Boundary tables: CSV / Parquet

A long-format vertex table with one row per polygon vertex and columns:

| Column | Type | Meaning |
| --- | --- | --- |
| `vertex_x` | float | vertex x coordinate |
| `vertex_y` | float | vertex y coordinate |
| `label_id` | int | segment id (alternative: `cell_id`) |
| `cell_id` | string | segment label; used when `label_id` is absent |

```bash
baysor run -m 50 -c configs/xenium.toml -o out data/experiment.xenium data/cell_boundaries.parquet
```

Polygons outside the molecule bounds are ignored. Molecules not covered by
any polygon are unassigned.

!!! note
    Image and boundary priors assign molecules by their x/y position. For 3D
    segmentations, prefer transcript-native `:column_name` priors.

## Prior behavior

- **Segment filtering.** Prior segments with fewer than
  `min_molecules_per_segment` assigned molecules are treated as unassigned.
  The default is `max(min_molecules_per_cell / 4, 2)`.
- **Confidence.** `--prior-segmentation-confidence` (default `0.2`, range
  `[0, 1]`) controls how strongly the segmentation must adhere to the prior:
  `0` ignores the prior, `1` forbids contradicting it. For sparse protocols
  (ISS, DARTFISH) or when the prior is high quality, values above `0.7` are
  recommended; otherwise the default works well.
- **Scale estimation.** Unless `--scale` is given, Baysor estimates the cell
  scale and `--scale-std` from the prior segments
  (`estimate_scale_from_prior`, on by default). Passing an explicit `--scale`
  disables this.
- **Unassigned molecules.** Molecules without a prior segment are segmented
  freely; their cell ids can be influenced with `--unassigned-prior-label`
  for transcript-native priors.

## When to use which prior

Use `:column_name` when prior assignments exist per transcript already — it is
the fastest mode and gives a clean `xeniumranger import-segmentation` handoff.

Use boundary tables when you have published cell or nucleus outlines.

Use image masks when you have labeled pixels, e.g. DAPI/watershed masks.

Use no prior for fully de novo segmentation — then `--scale` must be set
explicitly.
