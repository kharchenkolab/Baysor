# Documentation maintenance

The [user guide](index.md) is built with MkDocs + Material and versioned with
mike. Navigation is in [`mkdocs.yml`](../mkdocs.yml); tooling versions are
pinned in [`requirements.txt`](requirements.txt).

## Local checks

Run from the repository root:

```bash
python3 -m venv /tmp/mkdocs-venv
/tmp/mkdocs-venv/bin/python -m pip install -r docs/requirements.txt
/tmp/mkdocs-venv/bin/mkdocs build --strict
python3 docs/tools/check_cli_docs.py
```

The build writes `site/`. Use `/tmp/mkdocs-venv/bin/mkdocs serve` for a local
preview. Neither check requires a Baysor build.

## Editing

- Keep install-and-run commands at the top of `README.md` and `index.md`,
  and the key parameters at the top of `run.md`.
- Keep usage pages short; put full option tables in `cli.md`, config keys in
  `configuration.md`, and file schemas in `output_files.md`.
- Check CLI options against [`src/cli/main.cpp`](../src/cli/main.cpp), config
  behavior against [`src/utils/options.cpp`](../src/utils/options.cpp), and
  defaults against [`options.h`](../include/baysor/utils/options.h). Check
  output contracts in [`src/reporting/output.cpp`](../src/reporting/output.cpp).
  The CLI checker validates names, not defaults, semantics or output schemas.
- Keep performance measurements, tables and figures unchanged during prose
  edits. Verify release URLs when the release assets are available.

The [docs workflow](../.github/workflows/docs.yml) handles release publishing.
This maintainer README is excluded from the published user guide.
