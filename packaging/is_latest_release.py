#!/usr/bin/env python3
"""Print "true" if a release is the newest stable one and gets the `latest` tag.

    packaging/is_latest_release.py VERSION [--prerelease] [--repo OWNER/REPO]

Same rule as the docs site's `latest` alias (.github/workflows/docs.yml): the
release is not a pre-release and VERSION is >= the version of every
published, non-draft, non-prerelease release, as listed by `gh release list`.
Tags become versions as in packaging/release_version.sh (cpp-0.8.3, cpp-v0.8.0
and v0.9.0 give 0.8.3, 0.8.0 and 0.9.0), so the Julia-era v0.* releases are
simply older versions; tags that do not parse are ignored.
"""

import argparse
import json
import re
import subprocess

VERSION_RE = re.compile(r"^(\d+(?:\.\d+)*)(?:[-.]?(.+))?$")


def version_key(version):
    """Comparable key of a version ('0.8.3', '0.7', '0.9.0-rc1'), or None.

    A suffix sorts before the same version without one ('0.9.0-rc1' < '0.9.0').
    """
    m = VERSION_RE.match(version.strip())
    if not m:
        return None
    return tuple(int(p) for p in m.group(1).split(".")), not m.group(2), m.group(2) or ""


def version_of_tag(tag):
    return version_key(re.sub(r"^(cpp-)?v?", "", tag.strip()))


def is_latest(version, is_prerelease, releases):
    """releases: the JSON of `gh release list --json tagName,isDraft,isPrerelease`."""
    new = version_key(version)
    if new is None:
        raise ValueError(f"version does not parse: {version!r}")
    stable = (version_of_tag(r["tagName"]) for r in releases
              if not r.get("isDraft") and not r.get("isPrerelease"))
    return not is_prerelease and all(new >= key for key in stable if key is not None)


def fetch_releases(repo):
    cmd = ["gh", "release", "list", "--limit", "1000", "--json", "tagName,isDraft,isPrerelease"]
    if repo:
        cmd += ["--repo", repo]
    return json.loads(subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True).stdout)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("version", help="bare release version, e.g. 0.8.3")
    parser.add_argument("--prerelease", action="store_true",
                        help="the release is a GitHub pre-release")
    parser.add_argument("--repo", help="OWNER/REPO (default: the repository of the current directory)")
    args = parser.parse_args(argv)
    print("true" if is_latest(args.version, args.prerelease, fetch_releases(args.repo)) else "false")


if __name__ == "__main__":
    main()
