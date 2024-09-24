# Advanced configuration

The pipeline options are described in the CLI help. Run `baysor --help` or the corresponding command like `baysor run --help` for the list of main options.

However, there are additional parameters that can be specified through the TOML config. See [example_config.toml](https://github.com/kharchenkolab/Baysor/blob/master/configs/example_config.toml) for their description. All parameters from the config can also be passed through the command line. For example, to set `exclude_genes` from the `data` section you need to pass `--config.data.exclude-genes='Blank*,MALAT1'` parameter. Please, keep in mind that the CLI parameters require replacing all underscores (`_`) with `-`.

For more details on the syntax for CLI arguments, see the [Comonicon documentation](https://comonicon.org/stable/conventions/#Syntax-and-Conventions). TL;DR, possible spelling options are: `--x-column X` or `-x X`, also `--x-column=X` or `-xX`. And if you're using strings with unusual symbols like `*` or `?`, it's better to have them in quotes: `--config.data.exclude-genes='Blank*'`.

## Multi-threading

All running options support some basic multi-threading. To enable it, set `JULIA_NUM_THREADS` environment variable before running the script You can either do it globally by running `export JULIA_NUM_THREADS=13` or for an individual command:

```bash
JULIA_NUM_THREADS=13 baysor run -m 30 -s 30 ./molecules.csv
```

The latest julia version is recommended, as multi-threading is being actively developed in julia.