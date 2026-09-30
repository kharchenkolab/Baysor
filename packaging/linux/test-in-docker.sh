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
#
# Set BAYSOR_SMOKE_DATA=/path/to/molecules.parquet (sim_circles_gaps_g100 from
# github.com/VPetukhov/baysor-benchmarks) to run it instead of the synthetic grid.
# BAYSOR_TEST_STAGES selects stages (default: "native old-distros qemu").
set -euo pipefail

SRC_DIR="$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
[[ $# -eq 1 ]] || { echo "usage: $0 <baysor-*-linux-x86_64.tar.gz>" >&2; exit 2; }
ARCHIVE="$(CDPATH='' cd -- "$(dirname -- "$1")" && pwd)/$(basename -- "$1")"
STAGES="${BAYSOR_TEST_STAGES:-native old-distros qemu}"

name="$(basename "$ARCHIVE" .tar.gz)"               # baysor-<version>-linux-x86_64
version="${name#baysor-}"; version="${version%-linux-x86_64}"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/baysor-release-test.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
tar -xzf "$ARCHIVE" -C "$WORK"
[[ -x "$WORK/$name/bin/baysor" ]] || { echo "archive has no $name/bin/baysor" >&2; exit 1; }
echo "==> archive contents"
(cd "$WORK" && find "$name" -maxdepth 2 | sort)

mounts=(-v "$WORK/$name:/opt/baysor:ro" -v "$SRC_DIR/packaging:/packaging:ro")
data_env=()
if [[ -n "${BAYSOR_SMOKE_DATA:-}" ]]; then
    data="$(CDPATH='' cd -- "$(dirname -- "$BAYSOR_SMOKE_DATA")" && pwd)/$(basename -- "$BAYSOR_SMOKE_DATA")"
    mounts+=(-v "$data:/data/molecules.parquet:ro")
    data_env=(-e BAYSOR_SMOKE_DATA=/data/molecules.parquet)
fi

for stage in $STAGES; do
    case "$stage" in
        native)
            echo "==> native ($(uname -m), $(ldd --version 2>/dev/null | head -n1))"
            "$SRC_DIR/packaging/smoke_test.sh" "$WORK/$name/bin/baysor" "$version"
            ;;
        old-distros)
            for image in almalinux:8 debian:10; do
                echo "==> $image"
                docker run --rm "${mounts[@]}" "${data_env[@]+"${data_env[@]}"}" "$image" \
                    bash /packaging/smoke_test.sh /opt/baysor/bin/baysor "$version"
            done
            ;;
        qemu)
            echo "==> qemu-x86_64 -cpu qemu64"
            docker run --rm "${mounts[@]}" "${data_env[@]+"${data_env[@]}"}" \
                -e BAYSOR_SMOKE_WRAPPER="qemu-x86_64 -cpu qemu64" \
                debian:trixie bash -c '
                    set -e
                    apt-get update -qq && apt-get install -y -qq qemu-user > /dev/null
                    qemu-x86_64 --version | head -n1
                    echo "x86-64 ISA levels of the emulated CPU (glibc):"
                    qemu-x86_64 -cpu qemu64 /lib64/ld-linux-x86-64.so.2 --help \
                        | grep -E "^ +x86-64-v[234]" || true
                    bash /packaging/smoke_test.sh /opt/baysor/bin/baysor '"$version"
            ;;
        *)
            echo "unknown stage '$stage'" >&2; exit 2 ;;
    esac
done
echo "==> all stages passed"
