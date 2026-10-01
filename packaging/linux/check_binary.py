#!/usr/bin/env python3
"""Check that a Linux release binary is portable.

    packaging/linux/check_binary.py <bin/baysor>

Fails unless
  * every NEEDED entry is a glibc library,
  * no required GLIBC_ symbol version is newer than GLIBC_FLOOR (the glibc of
    the manylinux_2_28 build image),
  * no libstdc++/libgcc symbol versions are required (they are linked
    statically),
  * the binary is not marked as needing an x86-64 ISA level above the
    baseline (glibc >= 2.33 refuses to start such binaries on older CPUs).

Uses readelf from binutils. Prints a short report either way.
"""

import re
import subprocess
import sys

GLIBC_FLOOR = "2.28"
GLIBC_LIBS = {"libc.so.6", "libm.so.6", "libpthread.so.0", "libdl.so.2", "librt.so.1",
              "ld-linux-x86-64.so.2"}


def readelf(*args):
    return subprocess.run(["readelf", "-W", *args], check=True,
                          capture_output=True, text=True).stdout


def version_key(v):
    return tuple(int(x) for x in v.split("."))


def main():
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <bin/baysor>")
    binary = sys.argv[1]
    needed = re.findall(r"\(NEEDED\)\s+Shared library: \[([^\]]+)\]", readelf("-d", binary))
    versions = set(re.findall(r"Name: ([A-Z_]+[0-9.]*[0-9])\s", readelf("-V", binary)))
    glibc = sorted((v[len("GLIBC_"):] for v in versions if re.fullmatch(r"GLIBC_[0-9.]+", v)),
                   key=version_key)
    cxx_versions = sorted(v for v in versions if v.startswith(("GLIBCXX_", "CXXABI_", "GCC_")))
    isa = [s.strip() for s in re.findall(r"x86 ISA needed: (.*)", readelf("-n", binary))]

    errors = [f"NEEDED {lib} is not a glibc library" for lib in needed if lib not in GLIBC_LIBS]
    newer = [v for v in glibc if version_key(v) > version_key(GLIBC_FLOOR)]
    if newer:
        errors.append(f"requires GLIBC_ versions above the floor {GLIBC_FLOOR}: {', '.join(newer)}")
    if cxx_versions:
        errors.append("requires libstdc++/libgcc symbol versions: " + ", ".join(cxx_versions))
    errors += [f"binary is marked as needing x86 ISA level(s): {levels}" for levels in isa
               if any(l.strip() not in ("x86-64-baseline", "") for l in levels.split(","))]

    print(f"binary:         {binary}")
    print(f"NEEDED:         {', '.join(needed) or '(none)'}")
    print(f"max GLIBC:      {glibc[-1] if glibc else 'none'} (floor {GLIBC_FLOOR})")
    print(f"x86 ISA needed: {'; '.join(isa) or '(not marked)'}")
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
