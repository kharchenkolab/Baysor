#!/usr/bin/env bash
# Check that a macOS release binary is portable.
#
#   packaging/macos/check_binary.sh <bin/baysor> <deployment-target>
#
# Fails unless the binary is arm64, its minimum macOS version (LC_BUILD_VERSION
# minos) is at most <deployment-target>, and every linked dylib is a system
# library (/usr/lib, /System) or ships in ../lib (via @rpath/@loader_path).
set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: $0 <binary> <deployment-target>" >&2; exit 2; }
BIN="$1"; TARGET="$2"
LIB_DIR="$(cd "$(dirname "$BIN")/.." && pwd)/lib"
errors=0
fail() { echo "ERROR: $*" >&2; errors=$((errors + 1)); }
version_le() { [[ "$(printf '%s\n%s\n' "$1" "$2" | sort -t. -k1,1n -k2,2n -k3,3n | head -n1)" == "$1" ]]; }

archs="$(lipo -archs "$BIN")"
minos="$(otool -l "$BIN" | awk '/LC_BUILD_VERSION/ {f = 1} f && $1 == "minos" {print $2; exit}')"
echo "binary:  $BIN"
echo "archs:   $archs"
echo "minos:   ${minos:-unknown} (deployment target $TARGET)"
echo "dylibs:"
[[ " $archs " == *" arm64 "* ]] || fail "binary is not arm64"
if [[ -z "$minos" ]] || ! version_le "$minos" "$TARGET"; then
    fail "minimum macOS version ${minos:-unknown} is newer than $TARGET"
fi

while read -r lib _; do
    echo "  $lib"
    case "$lib" in
        /usr/lib/*|/System/Library/*) ;;
        @rpath/*|@loader_path/*|@executable_path/*)
            [[ -e "$LIB_DIR/$(basename "$lib")" ]] || fail "$lib is not shipped in $LIB_DIR" ;;
        *) fail "$lib is not a system library" ;;
    esac
done < <(otool -L "$BIN" | tail -n +2)

exit $(( errors > 0 ))
