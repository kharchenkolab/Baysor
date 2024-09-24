Baysor can be used in several ways:

## Cell segmentation

A minimal command for cell segmentation:

```bash
baysor run [-s SCALE -x X_COL -y Y_COL -z Z_COL --gene GENE_COL -c config.toml -o OUTPUT_PATH] MOLECULES_FILE [PRIOR_SEGMENTATION]
```

## Dataset preview

As a full run takes some time, it can be useful to run a quick preview to get meaning from the data and to get some guesses about the parameters of the full run.

```bash
baysor preview [-x X_COL -y Y_COL --gene GENE_COL -c config.toml -o OUTPUT_PATH] MOLECULES_FILE
```


## Segmentation-free analysis

Many analyses don't require segmentation, and can be run on local neighborhoods instead. In the paper, we call them Neighborhood Composition Vectors (NCVs). To obtain them from Baysor, you may run `baysor segfree`. For more information, see `baysor segfree --help`. Minimal command:

```bash
baysor segfree [-k K_NEIGHBORS -x X_COL -y Y_COL --gene GENE_COL -c config.toml -o OUTPUT_PATH] MOLECULES_FILE
```

