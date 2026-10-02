# Prior segmentation

A prior, such as a DAPI nuclear mask or existing cell labels, helps when gene
expression alone is not enough. Pass it after the molecule table:

```bash
baysor run -m 30 --prior-segmentation-confidence 0.5 \
  -o out molecules.csv :cell_id
```

Here `:cell_id` selects a column of prior labels; it is not a filename.
Baysor estimates scale from the prior unless you supply `-s` / `--scale`.

## Prior behavior

`--prior-segmentation-confidence` controls trust in the prior (default `0.2`):

- `0` removes the prior's assignment penalty. The prior can still affect scale
  estimation and initialization; omit it and set `-s` for a fully de novo run.
- Values between `0` and `1` allow Baysor to revise the prior. For sparse
  protocols such as ISS / DARTFISH, or a high-quality prior, values above
  `0.7` are a useful starting point.
- `1` keeps each retained prior segment together, but does not preserve its
  label or forbid merging. See the [guarantee below](#what-confidence-1-guarantees).

Before fitting, prior segments with fewer than `min_molecules_per_segment`
molecules are removed; their molecules become unassigned. The default is
`max(min_molecules_per_cell / 4, 2)` (integer division). Set it in `[prior]`
in the [config](configuration.md#prior) if needed.

Without an explicit positive scale, Baysor estimates scale and scale
variation from the retained prior. If estimation fails, provide `-s`.
Unassigned molecules are segmented freely and may join cells or go to noise.

## Prior types

### Molecule labels: `:column_name`

Use a column when assignments already exist per molecule, as in the example
above. This is the recommended mode for Xenium; see its
[preset and workflow](xenium.md).

`--unassigned-prior-label` identifies the value meaning no assignment
(default `0`; the Xenium preset uses `UNASSIGNED`). It does not control final
cell names.

### Image masks: TIFF

```bash
baysor run -m 30 --prior-segmentation-confidence 0.5 \
  -o out molecules.csv dapi_mask.tif
```

Use a single-channel integer-labeled mask, or a binary mask with one connected
foreground component per segment. Zero is background. Molecule x/y positions
must use the mask's pixel coordinate system. Non-segmented staining images
must be segmented first, for example with watershed in ImageJ or Cellpose.

### Boundary tables: CSV / Parquet

```bash
baysor run -m 30 --prior-segmentation-confidence 0.5 \
  -o out molecules.csv cell_boundaries.parquet
```

Use a long-format table with one row per polygon vertex:

| Column | Type | Meaning |
| --- | --- | --- |
| `vertex_x`, `vertex_y` | float | Vertex coordinates in the molecule coordinate system. |
| `label_id` | int | Segment ID. |
| `cell_id` | string | Alternative segment label, used when `label_id` is absent. |

Molecules outside the polygons are unassigned. Image and boundary priors use
x/y only; prefer molecule-label priors for 3D data.

## What confidence 1 guarantees

Each prior segment kept by `min_molecules_per_segment` is assigned as a
whole: its molecules are never split between final cells or partially
assigned to noise. This includes exact duplicate coordinates within the
same prior segment.

- Final cell IDs are run-local, not the prior labels. A prior cell may grow
  or merge with another prior cell or unassigned molecules.
- Molecules with no prior assignment are unconstrained.
- Segments removed by the size threshold are unconstrained; lower
  `min_molecules_per_segment` to keep smaller segments.
- A segment whose molecules are all isolated in the molecule graph stays
  wholly unassigned, rather than being split.

The grouping constraint can move molecules farther than the model alone
would. Use confidence `1` only when you trust the prior at the molecule level;
values below `1` use a soft penalty and may split prior cells.

??? note "How the hard constraint is applied"

    Baysor projects each prior segment onto the component holding the largest
    share of its molecules (ties use component ID). A segment with no assigned
    molecule remains wholly unassigned. This prevents the partial splitting
    reported for duplicate molecules in
    [issue #117](https://github.com/kharchenkolab/Baysor/issues/117).
