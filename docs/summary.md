# Summary

Baysor can be used in several ways.

## Cell segmentation

A minimal command for cell segmentation:

```bash
baysor run -m MIN_MOLECULES [-s SCALE -x X_COL -y Y_COL -z Z_COL -g GENE_COL -c config.toml -o OUTPUT_DIR] MOLECULES_FILE [PRIOR_SEGMENTATION]
```

See [Cell segmentation](run.md).

## Dataset preview

A full run takes some time, so a quick preview helps to understand the data
and to choose the parameters of the full run:

```bash
baysor preview -m MIN_MOLECULES [-x X_COL -y Y_COL -g GENE_COL -c config.toml -o OUTPUT_FILE] MOLECULES_FILE
```

See [Dataset preview](preview.md).

## Segmentation-free analysis

Many analyses don't need a segmentation and can use local neighborhoods
instead, the Neighborhood Composition Vectors (NCVs):

```bash
baysor segfree -m MIN_MOLECULES [-k K_NEIGHBORS -x X_COL -y Y_COL -g GENE_COL -c config.toml -o OUTPUT_FILE] MOLECULES_FILE
```

See [Segmentation-free analysis](segfree.md) and `baysor segfree --help`.
