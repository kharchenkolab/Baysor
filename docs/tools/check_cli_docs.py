#!/usr/bin/env python3
"""Check that the documentation matches the Baysor CLI and config parser.

Run from the repository root (no Baysor build needed):

    python docs/tools/check_cli_docs.py

The checker fails (exit code 1) when:

1. the docs mention a ``--long-option`` that none of the documented command
   line tools define (Baysor's CLI in ``src/cli/main.cpp`` and the build
   front-end ``configure.sh``; options of external tools referenced in the
   docs are listed in ``EXTERNAL_OPTIONS`` below);
2. the docs mention a config key that ``src/utils/options.cpp`` does not
   parse. A "config key mention" is a snake_case identifier inside inline
   code spans or a top-level ``key = value`` line in a fenced ``toml`` block;
3. a fenced ``toml`` block uses an unknown key.

It warns (exit code 0) about Baysor CLI options and config keys that exist in
the sources but are never mentioned in the docs.

Limitations: single-letter options are not checked, and snake_case identifiers
that are not config keys (column names, cell-id patterns, ...) must be listed
in ``KNOWN_IDENTIFIERS``. Regions can be excluded from all checks with

    <!-- check-cli-docs: ignore-start -->
    ...
    <!-- check-cli-docs: ignore-end -->
"""

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
DOCS = REPO / "docs"

# --- options of external tools that the docs legitimately reference -------
EXTERNAL_OPTIONS = {
    # xeniumranger import-segmentation (10x Xenium Ranger)
    "--id", "--xenium-bundle", "--transcript-assignment", "--viz-polygons",
    "--units", "--localcores", "--localmem",
    # apt-get / lsb_release snippets in the installation page
    "--no-install-recommends", "--short", "--codename",
    # mike (docs versioning, see docs/tools/migrate_gh_pages.py)
    "--push", "--update-aliases",
    # git / docker / mkdocs / sha256sum snippets
    "--no-local", "--rm", "--strict", "--check", "--ignore-missing",
}

# --- top-level flags shared by every CLI11 command -------------------------
# `--version` ships with the release-binary packaging; `--help` is generated
# by CLI11.
UNIVERSAL_OPTIONS = {"--help", "--version"}

# --- snake_case identifiers that are not config keys ----------------------
KNOWN_IDENTIFIERS = {
    # input / output column names
    "transcript_id", "cell_id", "label_id", "vertex_x", "vertex_y",
    "feature_name", "x_location", "y_location", "z_location",
    "is_noise", "ncv_color", "assignment_confidence", "max_cluster_frac",
    "avg_confidence", "avg_assignment_confidence", "n_transcripts",
    "n_vertices",
    # CLI positionals and accepted enum values, not config keys
    "prior_segmentation", "ica_mrf",
    # GitHub Actions event names referenced from the docs pages
    "workflow_dispatch",
    # gene-name pattern examples
    "antisense_",
}
KNOWN_IDENTIFIER_PATTERNS = [
    re.compile(p) for p in (r"cell_\d+", r"V\d+", r"z_\d+")
]

IGNORE_START = "<!-- check-cli-docs: ignore-start -->"
IGNORE_END = "<!-- check-cli-docs: ignore-end -->"

OPTION_TOKEN = re.compile(r"(?<![\w-])--[a-z0-9][a-z0-9-]*")
SNAKE_TOKEN = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)+$")
TOML_KEY_LINE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=")


def extract_baysor_options():
    """Map each CLI11 option of src/cli/main.cpp to its long name."""
    src = (REPO / "src" / "cli" / "main.cpp").read_text()
    longs = {}
    for kind in ("add_option", "add_flag"):
        for m in re.finditer(r"->%s\(\s*\"([^\"]+)\"" % kind, src):
            for name in m.group(1).split(","):
                name = name.strip()
                if name.startswith("--"):
                    longs[name] = kind
    return longs


def extract_configure_sh_options():
    """Long options understood by configure.sh (from its usage/case text)."""
    src = (REPO / "configure.sh").read_text()
    return set(OPTION_TOKEN.findall(src))


def extract_config_keys():
    """Config sections and keys parsed by src/utils/options.cpp."""
    src = (REPO / "src" / "utils" / "options.cpp").read_text()
    keys = set(re.findall(r"toml_get(?:_\w+)?\(\s*\w+\s*,\s*\"([^\"]+)\"", src))
    sections = set(re.findall(r"doc\.count\(\"([^\"]+)\"\)", src))
    sections |= set(re.findall(r"doc\[\"([^\"]+)\"\]", src))
    return sections, keys


def strip_ignored(text):
    out = []
    ignoring = False
    for line in text.splitlines():
        if line.strip() == IGNORE_START:
            ignoring = True
            continue
        if line.strip() == IGNORE_END:
            ignoring = False
            continue
        if not ignoring:
            out.append(line)
    return "\n".join(out)


def iter_md_blocks(text):
    """Yield (is_code_block, language, content) for prose and fenced blocks."""
    lines = text.splitlines()
    buf, lang, in_fence = [], "", False
    for line in lines:
        fence = line.strip().startswith("```")
        if fence and not in_fence:
            yield (False, "", "\n".join(buf))
            buf, lang, in_fence = [], line.strip()[3:].strip(), True
        elif fence and in_fence:
            yield (True, lang, "\n".join(buf))
            buf, lang, in_fence = [], "", False
        else:
            buf.append(line)
    yield (in_fence, lang, "\n".join(buf))


def is_known_identifier(token):
    if token in KNOWN_IDENTIFIERS:
        return True
    return any(p.fullmatch(token) for p in KNOWN_IDENTIFIER_PATTERNS)


def main():
    errors = []
    warnings = []

    baysor_options = extract_baysor_options()
    configure_options = extract_configure_sh_options()
    sections, config_keys = extract_config_keys()
    known_options = (set(baysor_options) | configure_options
                     | EXTERNAL_OPTIONS | UNIVERSAL_OPTIONS)

    mentioned_options = set()
    mentioned_keys = set()

    md_files = sorted(DOCS.rglob("*.md"))
    if not md_files:
        print("no markdown files found under docs/", file=sys.stderr)
        return 1

    for md in md_files:
        rel = md.relative_to(REPO)
        text = strip_ignored(md.read_text())
        for is_code, lang, content in iter_md_blocks(text):
            # 1. long options everywhere
            for tok in OPTION_TOKEN.findall(content):
                mentioned_options.add(tok)
                if tok not in known_options:
                    errors.append(f"{rel}: unknown option '{tok}'")

            if is_code:
                # 3. top-level keys of toml blocks must be config keys
                if lang == "toml":
                    for line in content.splitlines():
                        m = TOML_KEY_LINE.match(line)
                        if m:
                            key = m.group(1)
                            mentioned_keys.add(key)
                            if key not in config_keys:
                                errors.append(
                                    f"{rel}: unknown config key '{key}' "
                                    f"in toml block")
                continue

            # 2. snake_case identifiers inside inline code spans
            for span in re.findall(r"`([^`\n]+)`", content):
                span = span.strip()
                if span in config_keys:
                    mentioned_keys.add(span)
                if SNAKE_TOKEN.fullmatch(span):
                    mentioned_keys.add(span)
                    if not (span in config_keys or is_known_identifier(span)):
                        errors.append(f"{rel}: unknown config key '{span}'")

    # warnings for undocumented options / config keys
    for opt in sorted(set(baysor_options) - mentioned_options):
        warnings.append(f"undocumented CLI option: {opt}")
    for key in sorted(config_keys - mentioned_keys):
        warnings.append(f"undocumented config key: {key}")

    for w in warnings:
        print(f"WARNING: {w}")
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)

    print(f"\nchecked {len(md_files)} pages: "
          f"{len(baysor_options)} CLI options, {len(config_keys)} config keys "
          f"({len(errors)} errors, {len(warnings)} warnings)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
