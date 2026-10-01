#!/usr/bin/env python3
"""One-time migration of the docs gh-pages branch to a versioned layout.

The migration moves the Julia (Documenter.jl) site from ``dev/`` to ``0.7.1/``,
registers it in mike's ``versions.json`` as ``0.7.1 (Julia)``, and replaces
every ``dev/**.html`` page with a redirect stub to the same page under
``0.7.1/`` so old deep links keep working. It creates a commit on the local
gh-pages checkout and NEVER pushes.

Usage, from a scratch clone of the repository (the script itself is run from
any checkout that has the docs sources, e.g. the main one):

    # scratch clone whose working tree shows the gh-pages site
    git clone --no-local /path/to/Baysor baysor-ghpages
    cd baysor-ghpages
    git branch gh-pages origin/gh-pages
    git checkout gh-pages

    # 1. migrate the Julia site (creates a commit, does NOT push)
    python /path/to/Baysor/docs/tools/migrate_gh_pages.py --site-dir .

    # 2. optional: add the current docs as version 0.9.0 the way
    #    .github/workflows/docs.yml does on a release, minus --push
    git checkout <docs-source-branch>
    mike deploy --update-aliases 0.9.0 latest
    mike set-default latest

    # 3. inspect, then publish yourself (or let the docs workflow do it on
    #    the next release)
    git checkout gh-pages
    python -m http.server 8000
"""

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

VERSION = "0.7.1"
TITLE = "0.7.1 (Julia)"

REDIRECT_STUB = """\
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Redirecting</title>
  <link rel="canonical" href="{href}">
  <noscript>
    <meta http-equiv="refresh" content="0; url={href}" />
  </noscript>
  <script>
    window.location.replace(
      "{href}" + window.location.search + window.location.hash
    );
  </script>
</head>
<body>
  Redirecting to <a href="{href}">{href}</a>&hellip;
</body>
</html>
"""

DOC_VERSIONS_JS = """\
// Documentation version selector (legacy Documenter.jl site).
// The Julia documentation lives in the {version}/ directory since the
// migration to the versioned docs site; see versions.json for the mike
// version list used by the current site.
var DOC_VERSIONS = [
  "{version}",
];
var DOCUMENTER_NEWEST = "{version}";
var DOCUMENTER_STABLE = "{version}";
"""


def git(site_dir, *args):
    return subprocess.run(["git", "-C", str(site_dir), *args],
                          check=True, text=True, capture_output=True).stdout.strip()


def target_href(rel_path):
    """Relative URL from the redirect stub at dev/<rel_path> to its target.

    Pages (index.html) redirect to their directory, so that both /dev/run/ and
    /dev/run/index.html keep working.
    """
    if rel_path.name == "index.html":
        target, trailing = Path(VERSION) / rel_path.parent, "/"
    else:
        target, trailing = Path(VERSION) / rel_path, ""
    return os.path.relpath(target, Path("dev") / rel_path.parent).replace(os.sep, "/") + trailing


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site-dir", type=Path, default=Path("."),
                    help="working tree of the gh-pages checkout (default: .)")
    site_dir = ap.parse_args().site_dir.resolve()

    for required in ("dev", "index.html", "versions.js"):
        if not (site_dir / required).exists():
            raise SystemExit(f"error: {site_dir} does not look like the gh-pages site "
                             f"({required} is missing)")
    if (site_dir / VERSION).exists():
        raise SystemExit(f"error: {site_dir / VERSION} already exists")
    versions_json = site_dir / "versions.json"
    entries = json.loads(versions_json.read_text()) if versions_json.exists() else []
    if any(e.get("version") == VERSION for e in entries):
        raise SystemExit(f"error: version {VERSION!r} already in versions.json")
    try:
        status = git(site_dir, "status", "--porcelain")
    except subprocess.CalledProcessError:
        raise SystemExit(f"error: {site_dir} is not a git working tree")
    if status:
        raise SystemExit("error: working tree is not clean:\n" + status)

    shutil.move(site_dir / "dev", site_dir / VERSION)
    pages = sorted((site_dir / VERSION).rglob("*.html"))
    for html in pages:
        rel_path = html.relative_to(site_dir / VERSION)
        stub = site_dir / "dev" / rel_path
        stub.parent.mkdir(parents=True, exist_ok=True)
        stub.write_text(REDIRECT_STUB.format(href=target_href(rel_path)))
    entries.append({"version": VERSION, "title": TITLE, "aliases": []})
    versions_json.write_text(json.dumps(entries, indent=2) + "\n")
    (site_dir / "versions.js").write_text(DOC_VERSIONS_JS.format(version=VERSION))

    git(site_dir, "add", "-A")
    git(site_dir, "commit", "-m", f"docs: archive Julia site as {VERSION} and redirect dev/ to it")
    print(f"migrated dev/ -> {VERSION}/ and wrote {len(pages)} redirect stubs; commit created "
          f"on branch {git(site_dir, 'rev-parse', '--abbrev-ref', 'HEAD')}; NOT pushed")


if __name__ == "__main__":
    main()
