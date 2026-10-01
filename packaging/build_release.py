#!/usr/bin/env python3
"""Build a portable Baysor release archive for one platform.

    python3 packaging/build_release.py --platform linux-x86_64|macos-arm64|windows-x86_64

Steps: check out and bootstrap vcpkg at the ref in packaging/vcpkg-ref,
configure with the matching `release-<platform>` CMake preset (dependencies
come from vcpkg with the triplets in packaging/vcpkg-triplets), build and
install baysor, and pack

    dist/baysor-<version>-<platform>.tar.gz   (.zip on Windows)
      baysor-<version>-<platform>/bin/baysor[.exe]  (+ DLLs on Windows)
      baysor-<version>-<platform>/LICENSE
      baysor-<version>-<platform>/README.md

<version> is `project(baysor VERSION ...)` from CMakeLists.txt.

Linux builds must run in the packaging/linux/Dockerfile image (use
packaging/linux/build-in-docker.sh); the binary is then checked with
packaging/linux/check_binary.py.

vcpkg, its downloads and its binary cache live in --cache-dir (default
.release-cache/); a warm cache turns a rebuild into a Baysor-only compile.
"""

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent.parent

PLATFORMS = ("linux-x86_64", "macos-arm64", "windows-x86_64")

# Must match CMAKE_OSX_DEPLOYMENT_TARGET in the release-macos-arm64 preset and
# VCPKG_OSX_DEPLOYMENT_TARGET in packaging/vcpkg-triplets/arm64-osx-release.cmake.
MACOS_DEPLOYMENT_TARGET = "12.0"
# Installed by InstallRequiredSystemLibraries (BAYSOR_INSTALL_RUNTIME=ON).
WINDOWS_RUNTIME_DLLS = ["msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll"]


def log(msg):
    print(f"==> {msg}", flush=True)


def run(cmd, **kwargs):
    print("    $ " + " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def project_version():
    text = (SRC_DIR / "CMakeLists.txt").read_text()
    m = re.search(r"project\(\s*baysor\s+VERSION\s+([0-9][0-9A-Za-z.\-]*)", text)
    if not m:
        sys.exit("error: cannot find project(baysor VERSION ...) in CMakeLists.txt")
    return m.group(1)


# ---------------------------------------------------------------------------
# vcpkg
# ---------------------------------------------------------------------------

def ensure_vcpkg(cache_dir):
    """Clone vcpkg at packaging/vcpkg-ref into the cache and bootstrap it."""
    ref = (SRC_DIR / "packaging" / "vcpkg-ref").read_text().strip()
    root = cache_dir / "vcpkg"
    if not (root / ".git").is_dir():
        log(f"Cloning vcpkg {ref}")
        run(["git", "clone", "--quiet", "--depth", "1", "--branch", ref,
             "https://github.com/microsoft/vcpkg.git", root])
    else:
        head = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        want = subprocess.run(["git", "-C", root, "rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}"],
                              capture_output=True, text=True).stdout.strip()
        if head != want:
            log(f"Updating vcpkg checkout to {ref}")
            run(["git", "-C", root, "fetch", "--quiet", "--depth", "1", "origin", "tag", ref])
            run(["git", "-C", root, "checkout", "--quiet", "--detach", ref])

    exe = root / ("vcpkg.exe" if os.name == "nt" else "vcpkg")
    stamp = root / ".baysor-bootstrapped"
    if not exe.exists() or not stamp.exists() or stamp.read_text().strip() != ref:
        log("Bootstrapping vcpkg")
        if os.name == "nt":
            run(["cmd", "/c", root / "bootstrap-vcpkg.bat", "-disableMetrics"])
        else:
            run(["sh", root / "bootstrap-vcpkg.sh", "-disableMetrics"])
        stamp.write_text(ref + "\n")
    return root


# ---------------------------------------------------------------------------
# Packaging
# ---------------------------------------------------------------------------

def make_archive(stage_root, name, kind, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    archive = out_dir / f"{name}.{kind}"
    archive.unlink(missing_ok=True)
    top = stage_root / name
    if kind == "tar.gz":
        def normalize(info):
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            return info
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(top, arcname=name, filter=normalize)
    else:
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
            for path in sorted(top.rglob("*")):
                zf.write(path, (Path(name) / path.relative_to(top)).as_posix())
    return archive


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="\n".join(__doc__.splitlines()[2:]))
    ap.add_argument("--platform", required=True, choices=PLATFORMS)
    ap.add_argument("--out", type=Path, default=SRC_DIR / "dist",
                    help="directory for the archive (default: dist/)")
    ap.add_argument("--cache-dir", type=Path, default=SRC_DIR / ".release-cache",
                    help="vcpkg checkout, downloads and binary cache (default: .release-cache/)")
    ap.add_argument("--jobs", type=int, default=8, help="build parallelism (default: 8)")
    args = ap.parse_args()

    windows = args.platform.startswith("windows")
    version = project_version()
    name = f"baysor-{version}-{args.platform}"
    cache_dir = args.cache_dir.resolve()
    out_dir = args.out.resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    log(f"Building {name}")

    vcpkg_root = ensure_vcpkg(cache_dir)
    for d in ("vcpkg-bincache", "downloads"):
        (cache_dir / d).mkdir(exist_ok=True)
    os.environ.update({
        "VCPKG_ROOT": str(vcpkg_root),
        "VCPKG_DEFAULT_BINARY_CACHE": str(cache_dir / "vcpkg-bincache"),
        "VCPKG_DOWNLOADS": str(cache_dir / "downloads"),
        "VCPKG_DISABLE_METRICS": "1",
        "VCPKG_MAX_CONCURRENCY": str(args.jobs),
        "CMAKE_BUILD_PARALLEL_LEVEL": str(args.jobs),
    })

    preset = f"release-{args.platform}"
    build_dir = SRC_DIR / "build" / preset
    log(f"Configuring ({preset})")
    run(["cmake", "--preset", preset], cwd=SRC_DIR)
    log("Building")
    run(["cmake", "--build", build_dir, "--config", "Release", "--target", "baysor",
         "--parallel", str(args.jobs)])

    stage_root = build_dir / "package"
    shutil.rmtree(stage_root, ignore_errors=True)
    stage = stage_root / name
    log(f"Staging {stage}")
    strip = [] if windows else ["--strip"]
    run(["cmake", "--install", build_dir, "--config", "Release", "--prefix", stage, *strip])
    for f in ("LICENSE", "README.md"):
        shutil.copy2(SRC_DIR / f, stage / f)
    exe = stage / "bin" / ("baysor.exe" if windows else "baysor")
    if not exe.is_file():
        sys.exit(f"error: {exe} was not installed")

    if args.platform == "linux-x86_64":
        log("Checking portability")
        run([sys.executable, SRC_DIR / "packaging" / "linux" / "check_binary.py", exe])
    elif args.platform == "macos-arm64":
        log("Checking portability")
        run(["bash", SRC_DIR / "packaging" / "macos" / "check_binary.sh", exe, MACOS_DEPLOYMENT_TARGET])
    else:
        # Runners have the VC++ runtime in System32, so a smoke test alone
        # would not notice if it were missing from the archive.
        dlls = sorted(p.name.lower() for p in exe.parent.glob("*.dll"))
        print("    bundled DLLs: " + ", ".join(dlls))
        missing = [d for d in WINDOWS_RUNTIME_DLLS if d not in dlls]
        if missing:
            sys.exit("error: MSVC runtime DLLs missing from bin/: " + ", ".join(missing))

    archive = make_archive(stage_root, name, "zip" if windows else "tar.gz", out_dir)
    log(f"Wrote {archive}")
    print(f"{hashlib.sha256(archive.read_bytes()).hexdigest()}  {archive.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
