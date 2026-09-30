# Input data

Baysor consumes a molecule table plus an optional
[prior segmentation](priors.md) input.

## Supported molecule inputs

The `coordinates` positional argument of `run`, `preview`, and `segfree` may
be:

- a CSV molecule table
- a Parquet molecule table
- a Xenium `experiment.xenium` manifest

When `experiment.xenium` is passed, Baysor resolves the adjacent Xenium
transcript table automatically and keeps enough source context to make the
`legacy` outputs compatible with `xeniumranger import-segmentation` (see the
[Xenium workflow](xenium.md)).

## Molecule table columns

Required columns:

- `x`, `y` — molecule coordinates
- `gene` — gene name

Optional columns, used when present:

- `z` — third coordinate; enables 3D segmentation (ignored with `--force-2d`)
- `qv` — per-molecule quality value, filtered by `--min-qv` (column name set
  by `--qv-column`)
- `transcript_id` — preserved and written back to `legacy` output, and makes
  `segmentation.csv` compatible with `xeniumranger import-segmentation`
- `confidence` — precomputed molecule confidence, reused instead of estimating
  the noise model
- `cluster` — precomputed molecule-cluster labels
- the prior column referenced by `:column_name` (see
  [Prior segmentation](priors.md))

Column names are configurable through CLI flags (`-x`, `-y`, `-z`, `-g`,
`--qv-column`) or the config file ([Configuration](configuration.md)). For
Xenium, [configs/xenium.toml](https://github.com/kharchenkolab/Baysor/blob/HEAD/configs/xenium.toml)
already maps `x_location`, `y_location`, `z_location`, `feature_name`, and
`qv`.

If a numeric z column has only one unique value, it is dropped and the run
becomes 2D.

## Filtering during input load

Molecules are filtered while loading:

- `--min-qv` drops molecules with quality value below the threshold (no
  filtering by default);
- `--x-min` / `--x-max` / `--y-min` / `--y-max` / `--z-min` / `--z-max` crop
  the dataset spatially;
- `--min-molecules-per-gene` drops genes with too few molecules;
- `--exclude-genes` drops genes by name or glob pattern (`*`, `?`), e.g.
  `--exclude-genes 'Blank*,MALAT1'`.

Cropping and filtering are useful for quick development runs, protocol
debugging, and testing large datasets without loading the full field of view.
Note that a cropped run is not appropriate input for
`xeniumranger import-segmentation` (see [Xenium workflow](xenium.md)).

## 2D vs 3D data

Baysor segments in 3D when a z column with varying values is present (and
`--force-2d` is not set), in 2D otherwise. 3D runs write per-layer polygons
(`segmentation_polygons_3d.json` in the `legacy` style), 2D runs write joined
polygons (`segmentation_polygons_2d.json`).

## Protocol notes

### Xenium

Preferred input is the manifest:

```bash
baysor run -c configs/xenium.toml -o out data/experiment.xenium :cell_id
```

### ISS / osm-FISH / STARmap

These usually use CSV molecule tables plus either no prior with an explicit
`--scale`, an image mask prior, or a boundary table prior. See
[Examples](examples.md) for runnable commands.
