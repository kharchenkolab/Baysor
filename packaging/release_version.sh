#!/usr/bin/env bash
# Derive the release version from a tag and check it against the sources.
#
#   packaging/release_version.sh [<tag>]    # prints the version, e.g. 0.8.3
#
# The version is the tag with a leading "cpp-" and then "v" removed
# (cpp-0.8.3, cpp-v0.8.0 and v0.9.0 give 0.8.3, 0.8.0 and 0.9.0); without a
# tag (or with an empty one) it is the version of the sources. It must equal
# project(baysor VERSION ...) in CMakeLists.txt and "version-string" in
# vcpkg.json; otherwise the script fails.
set -euo pipefail

SRC_DIR="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
die() { printf 'release_version: error: %s\n' "$*" >&2; exit 1; }

[[ $# -le 1 ]] || die "usage: $0 [<tag>]"
cmake_version="$(sed -n 's/^project(baysor VERSION \([^ )]*\).*/\1/p' "$SRC_DIR/CMakeLists.txt")"
vcpkg_version="$(sed -n 's/^ *"version-string": "\([^"]*\)".*/\1/p' "$SRC_DIR/vcpkg.json")"
[[ -n "$cmake_version" ]] || die "cannot read project(baysor VERSION ...) from CMakeLists.txt"
tag="${1:-cpp-$cmake_version}"
version="${tag#cpp-}"
version="${version#v}"
[[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] \
    || die "tag '$tag' does not look like [cpp-][v]MAJOR.MINOR.PATCH"
[[ "$version" == "$cmake_version" ]] \
    || die "tag '$tag' gives version $version, but CMakeLists.txt has project(baysor VERSION $cmake_version)"
[[ "$version" == "$vcpkg_version" ]] \
    || die "tag '$tag' gives version $version, but vcpkg.json has version-string $vcpkg_version"
printf '%s\n' "$version"
