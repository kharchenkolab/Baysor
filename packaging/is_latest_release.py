#!/usr/bin/env python3
"""Decide whether a release is the newest stable one and gets the `latest` tag.

    packaging/is_latest_release.py VERSION [--prerelease]
                                    [--repo OWNER/REPO] [--releases-file FILE]

Prints "true" when the release VERSION should also carry the `latest` tag on
its Docker images, "false" otherwise; exits non-zero only on errors (unreadable
release list, unparsable VERSION).

The rule mirrors the docs site's `latest` alias (the "Decide whether to move
the `latest` alias" step in .github/workflows/docs.yml and should_move_latest
in docs/tools/migrate_gh_pages.py): the release must not be a pre-release and
its version must be >= the version of every published, non-draft,
non-prerelease release of the repository. Backport releases therefore never
move `latest`.

Releases are read with `gh release list --json tagName,isDraft,isPrerelease`
(--repo defaults to $GITHUB_REPOSITORY); GH_TOKEN/GH_HOST as set by gh are
respected. `--releases-file` reads the same JSON from a file instead — used
by the local tests in packaging/is_latest_release_test.py.

Tags become versions the way packaging/release_version.sh does: strip a
leading "cpp-", then a leading "v" (cpp-0.8.3, cpp-v0.8.0, v0.9.0 and the
Julia-era v0.7.1 give 0.8.3, 0.8.0, 0.9.0 and 0.7.1). Julia-era `v0.*` tags
need no special casing: they are simply older versions. Tags that do not
parse to a version are ignored.
"""

import argparse
import json
import os
import re
import subprocess
import sys

VERSION_RE = re.compile(r"^(\d+(?:\.\d+)*)(?:[-.]?(.+))?$")


def version_key(version):
    """Comparable key of a bare version ('0.8.3', 'v0.7.1', '0.7').

    Returns None when the version does not parse. Suffixes sort before the
    same version without a suffix ('0.9.0-rc1' < '0.9.0'), matching the docs
    workflow. Tuples of different lengths compare correctly ('0.4' < '0.8.3').
    """
    m = VERSION_RE.match(version.strip())
    if not m:
        return None
    return (tuple(int(p) for p in m.group(1).split(".")),
            not m.group(2), m.group(2) or "")


def version_of_tag(tag):
    """version_key() of a release tag (cpp-X.Y.Z / cpp-vX.Y.Z / vX.Y.Z)."""
    t = tag.strip()
    if t.startswith("cpp-"):
        t = t[4:]
    if t.startswith("v"):
        t = t[1:]
    return version_key(t)


def is_latest(version, is_prerelease, releases):
    """True if VERSION (a bare version like '0.8.3') is the newest stable.

    releases: iterable of dicts with tagName/isDraft/isPrerelease (the JSON
    of `gh release list --json tagName,isDraft,isPrerelease`). Drafts and
    pre-releases do not count; the candidate itself may be in the list (it is
    when the release already exists).
    """
    new = version_key(version)
    if new is None:
        raise ValueError(f"version does not parse: {version!r}")
    if is_prerelease:
        return False
    for release in releases:
        if release.get("isDraft") or release.get("isPrerelease"):
            continue
        key = version_of_tag(release["tagName"])
        if key is None:
            continue
        if new < key:
            return False
    return True


def fetch_releases(repo):
    """The repository's releases as `gh release list --json ...` emits them."""
    cmd = ["gh", "release", "list", "--limit", "1000",
           "--json", "tagName,isDraft,isPrerelease"]
    if repo:
        cmd += ["--repo", repo]
    try:
        out = subprocess.run(cmd, check=True, capture_output=True,
                             text=True).stdout
    except FileNotFoundError:
        raise SystemExit("is_latest_release: gh CLI not found on PATH")
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        raise SystemExit(f"is_latest_release: `gh release list` failed "
                         f"(exit {e.returncode})")
    return json.loads(out)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Print "true" if VERSION is the newest stable release '
                    'and should get the Docker `latest` tag.')
    parser.add_argument("version", help="bare release version, e.g. 0.8.3")
    parser.add_argument("--prerelease", action="store_true",
                        help="the release is a GitHub pre-release")
    parser.add_argument("--repo", default=None,
                        help="OWNER/REPO to list releases of "
                             "(default: $GITHUB_REPOSITORY)")
    parser.add_argument("--releases-file", metavar="FILE",
                        help="read the `gh release list --json` output from "
                             "FILE instead of running gh (for tests)")
    args = parser.parse_args(argv)

    repo = args.repo
    if repo is None:
        repo = os.environ.get("GITHUB_REPOSITORY")

    if args.releases_file:
        with open(args.releases_file) as f:
            releases = json.load(f)
    else:
        releases = fetch_releases(repo)

    print("true" if is_latest(args.version, args.prerelease, releases)
          else "false")
    return 0


if __name__ == "__main__":
    sys.exit(main())
