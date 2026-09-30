#!/usr/bin/env python3
"""One-time migration of the docs gh-pages branch to a versioned layout.

The migration moves the Julia (Documenter.jl) site from ``dev/`` to ``0.7.1/``,
registers it in mike's ``versions.json`` as ``0.7.1 (Julia)``, and replaces
every ``dev/**.html`` page with a redirect stub to the same page under
``0.7.1/`` so old deep links keep working. It creates a commit on the local
gh-pages checkout and NEVER pushes.

Optionally it can also deploy the current source checkout's docs as a given
version through mike and make ``latest`` the default (again: local commit
only).

Usage, from a scratch clone of the repository (the script itself is run from
any checkout that has the docs sources, e.g. the main one):

    # scratch clone whose working tree shows the gh-pages site
    git clone --no-local /path/to/Baysor baysor-ghpages
    cd baysor-ghpages
    git branch gh-pages origin/gh-pages
    git checkout gh-pages

    # 1. migrate the Julia site (creates a commit, does NOT push)
    python /path/to/Baysor/docs/tools/migrate_gh_pages.py --site-dir .

    # 2. optional: deploy this checkout's docs as version 0.8.4 with the
    #    `latest` alias and set `latest` as the default (local commit, no push)
    git checkout <docs-source-branch>
    python docs/tools/migrate_gh_pages.py --skip-migrate \
        --site-dir . --deploy-version 0.8.4

    # 3. inspect, then publish yourself (or let .github/workflows/docs.yml
    #    do it on the next release)
    git checkout gh-pages
    python -m http.server 8000

Steps 2-3 are exactly what the docs workflow automates on a release
(`mike deploy --push --update-aliases <version> latest` and
`mike set-default --push latest`); this script deliberately omits --push.
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_VERSION = "0.7.1"
DEFAULT_TITLE = "0.7.1 (Julia)"

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


def run_git(site_dir, *args, check=True):
    return subprocess.run(["git", "-C", str(site_dir), *args],
                          check=check, text=True, capture_output=True)


def loose_version_key(version):
    """Numeric-aware sort key, so 0.10 sorts after 0.9."""
    return tuple(int(p) if p.isdigit() else p
                 for p in re.split(r"(\d+)", version))


def target_href(version, rel_path):
    """Relative URL from the redirect stub at dev/<rel_path> to its target.

    Page paths end with index.html; those redirect to the containing
    directory so that both /dev/run/ and /dev/run/index.html keep working.
    """
    rel = Path(rel_path)
    if rel.name == "index.html":
        target = Path(version) / rel.parent
        trailing = "/"
    else:
        target = Path(version) / rel
        trailing = ""
    stub_dir = Path("dev") / rel.parent
    import os
    href = os.path.relpath(target, stub_dir).replace(os.sep, "/")
    return href + trailing


def write_redirect_stubs(site_dir, version):
    moved = site_dir / version
    stubs = 0
    for html in sorted(moved.rglob("*.html")):
        rel_path = html.relative_to(moved)
        stub = site_dir / "dev" / rel_path
        stub.parent.mkdir(parents=True, exist_ok=True)
        stub.write_text(REDIRECT_STUB.format(href=target_href(version, rel_path)))
        stubs += 1
    return stubs


def update_versions_json(site_dir, version, title):
    path = site_dir / "versions.json"
    entries = []
    if path.exists():
        entries = json.loads(path.read_text())
    if any(e.get("version") == version for e in entries):
        raise SystemExit(f"error: version {version!r} already in versions.json")
    entries.append({"version": version, "title": title, "aliases": []})
    entries.sort(key=lambda e: loose_version_key(e["version"]), reverse=True)
    path.write_text(json.dumps(entries, indent=2) + "\n")
    return entries


def migrate(site_dir, version, title, message):
    for required in ("dev", "index.html", "versions.js"):
        if not (site_dir / required).exists():
            raise SystemExit(
                f"error: {site_dir} does not look like the gh-pages site "
                f"({required} is missing)")
    if (site_dir / version).exists():
        raise SystemExit(f"error: {site_dir / version} already exists")
    if version in ("dev", "."):
        raise SystemExit(f"error: bad version name {version!r}")

    try:
        run_git(site_dir, "rev-parse", "--is-inside-work-tree")
    except subprocess.CalledProcessError:
        raise SystemExit(f"error: {site_dir} is not a git working tree")
    status = run_git(site_dir, "status", "--porcelain").stdout.strip()
    if status:
        raise SystemExit("error: working tree is not clean:\n" + status)

    shutil.move(str(site_dir / "dev"), str(site_dir / version))
    stubs = write_redirect_stubs(site_dir, version)
    entries = update_versions_json(site_dir, version, title)
    (site_dir / "versions.js").write_text(DOC_VERSIONS_JS.format(version=version))

    run_git(site_dir, "add", "-A")
    run_git(site_dir, "commit", "-m", message)
    print(f"migrated dev/ -> {version}/ and wrote {stubs} redirect stubs")
    print("versions.json now lists:")
    for e in entries:
        print(f"  - {e['version']} (title: {e['title']}, aliases: {e['aliases']})")
    print(f"commit created on branch "
          f"{run_git(site_dir, 'rev-parse', '--abbrev-ref', 'HEAD').stdout.strip()}; "
          "NOT pushed")


def deploy(source_dir, version):
    if not (source_dir / "mkdocs.yml").exists():
        raise SystemExit(f"error: {source_dir} has no mkdocs.yml")
    for cmd in (["mike", "deploy", "--update-aliases", version, "latest"],
                ["mike", "set-default", "latest"]):
        print("$ " + " ".join(cmd))
        subprocess.run(cmd, cwd=source_dir, check=True)
    print(f"deployed docs of {source_dir} as version {version} with the "
          f"'latest' alias and set 'latest' as the default; "
          f"commit created on the local gh-pages branch; NOT pushed")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site-dir", type=Path, default=Path("."),
                    help="working tree of the gh-pages checkout (default: .)")
    ap.add_argument("--version", default=DEFAULT_VERSION,
                    help=f"version directory for the Julia site (default: {DEFAULT_VERSION})")
    ap.add_argument("--title", default=DEFAULT_TITLE,
                    help=f"versions.json title (default: {DEFAULT_TITLE!r})")
    ap.add_argument("--message", default=None,
                    help="commit message for the migration commit")
    ap.add_argument("--skip-migrate", action="store_true",
                    help="skip the gh-pages file migration")
    ap.add_argument("--deploy-version", default=None,
                    help="also deploy the source checkout's docs as this "
                         "version with the 'latest' alias (mike, no push)")
    ap.add_argument("--source-dir", type=Path,
                    default=Path(__file__).resolve().parent.parent.parent,
                    help="docs source checkout used for --deploy-version "
                         "(default: the repository containing this script)")
    args = ap.parse_args()

    site_dir = args.site_dir.resolve()
    if not args.skip_migrate:
        migrate(site_dir, args.version, args.title,
                args.message or
                f"docs: archive Julia site as {args.version} and redirect dev/ to it")
    if args.deploy_version:
        deploy(args.source_dir.resolve(), args.deploy_version)
    if args.skip_migrate and not args.deploy_version:
        raise SystemExit("error: nothing to do (--skip-migrate without --deploy-version)")


if __name__ == "__main__":
    main()
