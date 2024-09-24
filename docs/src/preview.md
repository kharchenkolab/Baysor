## Preview

As a full run takes some time, it can be useful to run a quick preview to get meaning from the data and to get some guesses about the parameters of the full run. The only output of this step is `preview.html`, which visualizes the dataset and provides some diagnostics. Be careful, as this can get too large to be rendered for large datasets, so it's better to run on a subset of the data.

```bash
baysor preview <args> [options]
```

CLI parameters:

```@docs
Baysor.CLI.preview
```