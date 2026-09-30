#!/usr/bin/env python3
"""Check that a Linux release binary is portable.

    packaging/linux/check_binary.py <bin/baysor> [--glibc-floor 2.28]

Fails unless
  * every NEEDED entry is a glibc library or a library shipped in ../lib,
  * no required GLIBC_ symbol version is newer than the glibc floor,
  * no libstdc++/libgcc symbol versions are required (they are linked
    statically),
  * the binary is not marked as needing an x86-64 ISA level above the
    baseline (glibc >= 2.33 refuses to start such binaries on older CPUs).

Uses readelf from binutils. Prints a short report either way.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

GLIBC_LIBS = {
    "libc.so.6", "libm.so.6", "libpthread.so.0", "libdl.so.2", "librt.so.1",
    "ld-linux-x86-64.so.2", "ld-linux-aarch64.so.1",
}


def readelf(*args):
    return subprocess.run(["readelf", "-W", *args], check=True,
                          capture_output=True, text=True).stdout


def version_key(v):
    return tuple(int(x) for x in v.split("."))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("binary", type=Path)
    ap.add_argument("--glibc-floor", default="2.28",
                    help="newest GLIBC_ symbol version allowed (default: %(default)s)")
    args = ap.parse_args()

    errors = []
    dynamic = readelf("-d", str(args.binary))
    needed = re.findall(r"\(NEEDED\)\s+Shared library: \[([^\]]+)\]", dynamic)
    rpaths = re.findall(r"\((?:RPATH|RUNPATH)\)\s+Library (?:rpath|runpath): \[([^\]]+)\]", dynamic)
    lib_dir = args.binary.resolve().parent.parent / "lib"
    bundled = {p.name for p in lib_dir.iterdir()} if lib_dir.is_dir() else set()
    for lib in needed:
        if lib not in GLIBC_LIBS and lib not in bundled:
            errors.append(f"NEEDED {lib} is neither a glibc library nor shipped in {lib_dir}")

    versions = set(re.findall(r"Name: ([A-Z_]+[0-9.]*[0-9])\s", readelf("-V", str(args.binary))))
    glibc = sorted((v.split("_", 1)[1] for v in versions if v.startswith("GLIBC_")
                    and re.fullmatch(r"GLIBC_[0-9.]+", v)), key=version_key)
    max_glibc = glibc[-1] if glibc else "none"
    if glibc and version_key(max_glibc) > version_key(args.glibc_floor):
        newer = [v for v in glibc if version_key(v) > version_key(args.glibc_floor)]
        errors.append(f"requires GLIBC_{max_glibc} > floor {args.glibc_floor} "
                      f"(versions above the floor: {', '.join(newer)})")
    cxx_versions = sorted(v for v in versions if v.startswith(("GLIBCXX_", "CXXABI_", "GCC_")))
    if cxx_versions:
        errors.append("requires libstdc++/libgcc symbol versions: " + ", ".join(cxx_versions))

    isa = re.findall(r"x86 ISA needed: (.*)", readelf("-n", str(args.binary)))
    for levels in isa:
        if any(l.strip() not in ("x86-64-baseline", "") for l in levels.split(",")):
            errors.append(f"binary is marked as needing x86 ISA level(s): {levels.strip()}")

    print(f"binary:        {args.binary}")
    print(f"NEEDED:        {', '.join(needed) or '(none)'}")
    print(f"RPATH/RUNPATH: {', '.join(rpaths) or '(none)'}")
    print(f"bundled libs:  {', '.join(sorted(bundled)) or '(none)'}")
    print(f"max GLIBC:     {max_glibc} (floor {args.glibc_floor})")
    print(f"x86 ISA needed: {'; '.join(s.strip() for s in isa) or '(not marked)'}")
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
