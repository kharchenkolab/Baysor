# Input data

Use a molecule table with one row per detected molecule. `run`, `preview` and
`segfree` accept CSV, Parquet or a Xenium `experiment.xenium` manifest:

```bash
baysor run -m 30 -s 8 -o out molecules.csv
```

## Table columns

| Column | Required? | Meaning |
| --- | --- | --- |
| `x`, `y` | yes | Numeric molecule coordinates. |
| `gene` | yes | Gene name. |
| `z` | no | Numeric third coordinate; varying values enable 3D segmentation. |
| `qv` | no | Quality value; filtered only when `--min-qv` is nonnegative. |
| `transcript_id` | no | Numeric transcript ID, preserved in legacy output for Xenium Ranger. |
| A prior-label column | no | Existing cell assignments, selected with `:column_name`. |

Coordinates and scale must use the same units. Column names can be changed
with `-x`, `-y`, `-z`, `-g` and `--qv-column`, or in a
[config file](configuration.md). For example:

```bash
baysor run -m 30 -s 8 -x x_location -y y_location -g feature_name \
  -o out transcripts.parquet
```

Input `confidence` and `cluster` columns do not bypass fitting: the CLI
recalculates molecule confidence and computes its own clustering prior.

## Xenium

Prefer `experiment.xenium` to a transcript table alone. Baysor locates the
adjacent transcript table automatically. The
[Xenium workflow](xenium.md) includes the preset for `x_location`,
`y_location`, `z_location`, `feature_name` and `qv`, prior labels and the
Xenium Ranger handoff.

## Filtering and crops

Input filtering applies before segmentation:

- `--min-qv` drops low-quality molecules when the quality column is present
  (disabled by default).
- `--x-min` / `--x-max`, `--y-min` / `--y-max` and `--z-min` / `--z-max` keep
  a coordinate range.
- `--min-molecules-per-gene` drops genes with too few molecules.
- `--exclude-genes 'Blank*,MALAT1'` drops names or glob patterns (`*`, `?`).

The gene-filter flags are `run`-only; use the corresponding config keys for
`preview` or `segfree`. Crops are useful for choosing parameters, but use the
full dataset for [Xenium Ranger import](xenium.md#xenium-explorer-handoff).

## 2D and 3D

A varying z column enables 3D segmentation; missing or constant z gives a 2D
run. `--force-2d` ignores z. 3D runs write both per-layer polygons and joined
2D polygons pooled across the z-stack; see [Outputs](outputs.md).

For optional TIFF masks, boundary tables and molecule-label priors, see
[Prior segmentation](priors.md).
