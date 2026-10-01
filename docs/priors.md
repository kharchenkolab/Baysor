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
  `0` ignores the prior, `1` treats it as a hard constraint (see
  [What confidence 1 guarantees](#what-confidence-1-guarantees) below). For
  sparse protocols (ISS, DARTFISH) or when the prior is high quality, values
  above `0.7` are recommended; otherwise the default works well.
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

## What confidence 1 guarantees

At confidence 1 the prior is a hard constraint on grouping: every prior
segment is assigned as a whole. All molecules that share a prior label end
up in one final cell. In particular:

- A prior cell is never split across several final cells and is never
  partially dropped to noise.
- Exact duplicate molecules — rows with identical coordinates, as produced
  by the CosMx duplicate-transcript bug behind
  [#117](https://github.com/kharchenkolab/Baysor/issues/117) — always stay
  together with the rest of their prior cell.
- A prior cell may be renamed (final cell ids are run-local, they are not
  the prior labels), may grow, and may end up merged with another prior
  cell or with unassigned molecules. Confidence 1 forbids *splitting* a
  prior cell, not merging or expanding it.
- Molecules without a prior label (the `--unassigned-prior-label` value,
  `0` by default, or molecules outside the prior) are segmented by the
  model alone: they may join any cell or go to noise and are not part of
  the constraint.
- A whole prior segment is reported as unassigned when it has no molecules
  left in the algorithm: its label was filtered at input because the
  segment has fewer than `min_molecules_per_segment` molecules (default
  `max(min_molecules_per_cell / 4, 2)`) or equals
  `--unassigned-prior-label`. A segment whose molecules are all isolated in
  the molecule graph (no neighbors at all) also stays entirely unassigned,
  again together rather than split.

The constraint is implemented as a projection: when the model would scatter
a prior segment, the component that already holds the largest share of the
segment's molecules (ties broken by component id) receives all of them. This
can move individual molecules further than the model alone would, which is
the point of confidence 1 — use it only when the prior is trusted at the
single-molecule level. Values below 1 keep the original soft penalty, in
which a prior cell may be split.
