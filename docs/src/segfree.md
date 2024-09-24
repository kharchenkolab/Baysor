## Segmentation-free analysis

Many analyses don't require segmentation, and can be run on local neighborhoods instead. In the paper, we call them Neighborhood Composition Vectors (NCVs). To obtain them from Baysor, you may run `baysor segfree`. This would output a [loom file](https://linnarssonlab.org/loompy/format/index.html) with NCVs of size `k` in (`/matrix`), storing `ncv_color` and `confidence` as column attributes (`/col_attrs/`).

```bash
baysor segfree <args> [options]
```

CLI parameters:

```@docs
Baysor.CLI.segfree
```
