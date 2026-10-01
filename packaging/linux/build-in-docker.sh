#!/usr/bin/env bash
# Build the portable Linux x86_64 release archive in Docker.
#
#   packaging/linux/build-in-docker.sh [build_release.py options]
#
# Builds the image from packaging/linux/Dockerfile and runs
# packaging/build_release.py --platform linux-x86_64 inside it as the calling
# user. The archive is written to dist/. vcpkg, its downloads and its binary
# cache live in .release-cache/ (override with BAYSOR_RELEASE_CACHE), so
# rebuilds only recompile Baysor. Parallelism: BAYSOR_JOBS (default 8).
set -euo pipefail

SRC_DIR="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
CACHE_DIR="${BAYSOR_RELEASE_CACHE:-$SRC_DIR/.release-cache}"
IMAGE=baysor-release-builder:manylinux_2_28

mkdir -p "$CACHE_DIR"
CACHE_DIR="$(CDPATH='' cd -- "$CACHE_DIR" && pwd)"
docker build -t "$IMAGE" "$SRC_DIR/packaging/linux"

# The cache may live outside the source tree, so it is mounted separately.
exec docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache/home \
    -v "$SRC_DIR:/src" -v "$CACHE_DIR:/cache" -w /src "$IMAGE" \
    python3 packaging/build_release.py --platform linux-x86_64 --cache-dir /cache \
    --jobs "${BAYSOR_JOBS:-8}" "$@"
