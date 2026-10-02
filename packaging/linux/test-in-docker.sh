#!/usr/bin/env bash
# Test a Linux release archive where users will run it.
#
#   packaging/linux/test-in-docker.sh dist/baysor-<version>-linux-x86_64.tar.gz
#
# Extracts the archive and runs packaging/smoke_test.sh
#   1. natively on this machine,
#   2. in bare almalinux:8 and debian:10 containers (glibc 2.28, the floor;
#      nothing installed besides the base image),
#   3. under qemu-user with `-cpu qemu64` (x86-64 baseline: no SSSE3, SSE4,
#      AVX or FMA), in a Debian container that only adds qemu-user.
set -euo pipefail

SRC_DIR="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
[[ $# -eq 1 ]] || { echo "usage: $0 <baysor-*-linux-x86_64.tar.gz>" >&2; exit 2; }
name="$(basename "$1" .tar.gz)"               # baysor-<version>-linux-x86_64
version="${name#baysor-}"; version="${version%-linux-x86_64}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/baysor-release-test.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
tar -xzf "$1" -C "$WORK"
[[ -x "$WORK/$name/bin/baysor" ]] || { echo "archive has no $name/bin/baysor" >&2; exit 1; }
echo "==> archive contents"
(cd "$WORK" && find "$name" -maxdepth 2 | sort)

echo "==> native ($(uname -m), $(ldd --version 2>/dev/null | head -n1))"
"$SRC_DIR/packaging/smoke_test.sh" "$WORK/$name/bin/baysor" "$version"

mounts=(-v "$WORK/$name:/opt/baysor:ro" -v "$SRC_DIR/packaging:/packaging:ro")
for image in almalinux:8 debian:10; do
    echo "==> $image"
    docker run --rm "${mounts[@]}" "$image" bash /packaging/smoke_test.sh /opt/baysor/bin/baysor "$version"
done

echo "==> qemu-x86_64 -cpu qemu64"
docker run --rm "${mounts[@]}" -e BAYSOR_SMOKE_WRAPPER="qemu-x86_64 -cpu qemu64" debian:trixie bash -c '
    set -e
    apt-get update -qq && apt-get install -y -qq qemu-user > /dev/null
    qemu-x86_64 --version | head -n1
    echo "x86-64 ISA levels of the emulated CPU (glibc):"
    qemu-x86_64 -cpu qemu64 /lib64/ld-linux-x86-64.so.2 --help | grep -E "^ +x86-64-v[234]" || true
    bash /packaging/smoke_test.sh /opt/baysor/bin/baysor '"$version"
echo "==> all stages passed"
