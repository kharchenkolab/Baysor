# Releasing Baysor (C++ line)

Publishing a GitHub release is the only manual step that produces artifacts.
The `release` workflow (`.github/workflows/release.yml`) then builds, tests and
attaches the binaries. The documentation site is rebuilt by its own workflow
(`.github/workflows/docs.yml`) when the release is published, and again when
a pre-release is promoted to a full release (only full releases become the
default `latest` docs version).

## 1. Prepare the release commit

1. Set the new version `X.Y.Z` in both places (they must agree, and the
   release workflow fails if either differs from the tag):
   - `CMakeLists.txt`: `project(baysor VERSION X.Y.Z LANGUAGES C CXX)`
   - `vcpkg.json`: `"version-string": "X.Y.Z"`
2. In `CHANGELOG.md`, rename `## Unreleased` to
   `## [cpp-X.Y.Z] — YYYY-MM-DD`. Start a new `## Unreleased` section above it
   in the next commit that changes behaviour.
3. Check the tag you are going to use:

   ```bash
   packaging/release_version.sh cpp-X.Y.Z    # prints X.Y.Z or explains the mismatch
   ```

4. Optionally build and test the Linux binary locally (see
   [Reproducing the Linux build](#reproducing-the-linux-build)).
5. Commit, merge into the release branch (`cpp`) and push.

## 2. Tag and publish

Tags have the form `cpp-X.Y.Z` (the forms `cpp-vX.Y.Z` and `vX.Y.Z` are also
accepted: the version is the tag with a leading `cpp-` and then `v` removed).

```bash
git tag -a cpp-X.Y.Z -m "Baysor cpp-X.Y.Z"
git push origin cpp-X.Y.Z
gh release create cpp-X.Y.Z --verify-tag --title "Baysor X.Y.Z" --notes-file notes.md
```

`notes.md` is the release's section of `CHANGELOG.md`. You can also create the
release in the GitHub web UI. A draft release does not trigger the workflow;
publishing the draft does.

## 3. What the workflow does

On `release: published` the workflow

1. checks out the tag and fails if the tag version differs from
   `CMakeLists.txt` or `vcpkg.json`;
2. builds three archives in parallel with `packaging/build_release.py`:

   | Asset | Runner | Build |
   | --- | --- | --- |
   | `baysor-X.Y.Z-linux-x86_64.tar.gz` | `ubuntu-24.04` | `packaging/linux/build-in-docker.sh` (manylinux_2_28 container) |
   | `baysor-X.Y.Z-macos-arm64.tar.gz` | `macos-14` | `build_release.py --platform macos-arm64` |
   | `baysor-X.Y.Z-windows-x86_64.zip` | `windows-2022` | `build_release.py --platform windows-x86_64` |

3. smoke-tests every archive on its runner with `packaging/smoke_test.sh`
   (`--version`, `--help` and a full `baysor run` on a synthetic dataset). The
   Linux archive is additionally run in bare `almalinux:8` and `debian:10`
   containers and under `qemu-x86_64 -cpu qemu64`
   (`packaging/linux/test-in-docker.sh`);
4. writes `SHA256SUMS` and uploads the archives and `SHA256SUMS` to the
   release with `gh release upload --clobber`. If a platform fails, the
   archives that did build are still attached and the run fails, naming the
   missing platform; fix it and rerun the workflow for the tag
   (`workflow_dispatch`). vcpkg packages are cached even when a run fails, so
   the rerun does not rebuild all dependencies from scratch. When a build
   fails, the job also uploads vcpkg's per-port configure/build logs
   (`.release-cache/vcpkg/buildtrees/**/*.log`) as the `vcpkg-logs-<platform>`
   artifact, so the failing port's log (e.g. thrift's or gmp's) can be read
   without re-running;
5. in a parallel `docker` job, builds a small runtime image from the Linux
   archive, smoke-tests it and pushes it to GHCR and Docker Hub — see
   [Docker images](#docker-images) below. It never blocks step 4: each job
   fails independently, so a Docker problem does not keep the archives from
   being attached and vice versa.

Each archive contains one directory `baysor-X.Y.Z-<platform>/` with
`bin/baysor` (`bin/baysor.exe` plus the DLLs it needs on Windows), `LICENSE`
and `README.md`. Users can verify downloads with `sha256sum -c SHA256SUMS`
(`shasum -a 256 -c SHA256SUMS` on macOS).

vcpkg dependencies are built from source on the first run (roughly one to
three hours per platform) and cached with `actions/cache`; later releases only
compile Baysor unless `vcpkg.json`, `vcpkg-configuration.json`,
`packaging/vcpkg-ref`, `packaging/vcpkg-overlay-ports/`, the triplets, the
Linux Dockerfile or `build_release.py` change.

`vcpkg-configuration.json` registers `packaging/vcpkg-overlay-ports/`, which
currently overrides the `gmp` port: MSYS2 removes superseded package builds
from its mirrors, so the `autoconf2.71` package pinned inside vcpkg's gmp port
started to return 404 and broke every Windows build, and no vcpkg release
contains the upstream fix yet (microsoft/vcpkg#53437). The overlay is
gmp 6.3.0 from the `packaging/vcpkg-ref` baseline plus that one fix; remove
it (and `vcpkg-configuration.json`) once `packaging/vcpkg-ref` moves past
#53437.

## Rerunning the workflow

To rebuild the binaries of an existing release (for example after a failed
job), start the workflow by hand with the release tag:

```bash
gh workflow run release.yml -f tag=cpp-X.Y.Z
```

or use **Actions → release → Run workflow** in the web UI. The workflow
definition comes from the branch you run it on (the default branch unless you
pass `--ref`), but it always builds the sources of the tag, so the tag must
contain `packaging/` (every release after `cpp-0.8.3`). Existing assets are
replaced (`--clobber`). To rerun only failed jobs of a run, use **Re-run failed
jobs** on the run page.

## Dry-run builds

To exercise the full release build before tagging — for example to verify a
vcpkg bump or a dependency change — run the workflow in dry-run mode on any
ref; no release and no tag are needed:

```bash
gh workflow run release.yml --ref my-branch \
  -f dry_run=true [-f ref=my-branch] [-f platforms=macos-arm64]
```

- `ref` is a branch name or commit SHA to build; empty means the ref the
  workflow runs on. The version comes from `project(baysor VERSION ...)` in
  `CMakeLists.txt` (checked against `vcpkg.json`).
- `platforms` optionally limits the matrix to a comma-separated subset of
  `linux-x86_64`, `macos-arm64`, `windows-x86_64`, e.g. to iterate on the
  platform being fixed while its cache warms up. Empty builds all three.
- Nothing is uploaded to any release and no release is consulted. Each
  archive is uploaded as a workflow artifact (`baysor-X.Y.Z-<platform>`, kept
  for 7 days) and the run is green only when every selected platform builds
  and passes its smoke test. On failure the `vcpkg-logs-<platform>` artifact
  carries vcpkg's build logs.

Real releases (`release: published`, or `workflow_dispatch` with a `tag` and
without `dry_run`) are unaffected: the same version check, the same partial
upload that fails the run when a platform is missing, and the same
`SHA256SUMS`.

## Docker images

The `docker` job of the release workflow packages the Linux release archive
into a small runtime image and publishes it. It runs after all platform
builds and in parallel with `publish`, which it never waits on (see above).

**Image.** `packaging/docker/Dockerfile`: a pinned `debian:12-slim` base
(the release binary needs only glibc ≥ 2.28) with the extracted
`baysor-X.Y.Z-linux-x86_64/` tree installed — `bin/baysor` as
`ENTRYPOINT ["/usr/local/bin/baysor"]`, `LICENSE`/`README.md` under
`/usr/share/doc/baysor/`, OCI labels (`org.opencontainers.image.source`,
`.version`, `.licenses`, `.description`), `WORKDIR /data`, no compilers. It
runs as the non-root user `baysor` (UID/GID 1000); bind-mounted data
directories are then written as UID 1000 unless `docker run --user` is used
(documented in `docs/installation.md`). The build fails if
`baysor --version` inside the image does not print the release version.

**Registries and authentication.**

- GHCR: `ghcr.io/<repository owner in lower case>/baysor`
  (`ghcr.io/kharchenkolab/baysor` upstream), pushed with the workflow's
  `GITHUB_TOKEN`; `packages: write` is granted on the `docker` job only.
  A newly created GHCR package is **private until made public once** in the
  package settings (package page → *Package settings* → *Change visibility*)
  — do this after the first push if anonymous `docker pull` should work.
- Docker Hub: `vpetukhov/baysor`, or whatever the repository variable
  `DOCKERHUB_REPOSITORY` says (e.g. `yourname/baysor`), pushed with the
  secrets `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN` (a Docker Hub access
  token is fine; store both as repository secrets). If either secret is
  missing — forks, local runners — the job logs a warning and skips Docker
  Hub; GHCR still works, so forks are unaffected.

**Tags.** `X.Y.Z` always, plus `latest` when
`packaging/is_latest_release.py` decides the release is the newest stable
one: not a pre-release and its version is ≥ the version of every published
non-draft, non-prerelease release of the repository. This is the same rule
as the docs site's `latest` alias; Julia-era `v0.*` releases are simply
older versions and need no special handling. The script reads the release
list with `gh release list`; its local tests are
`python3 packaging/is_latest_release_test.py`.

**Smoke tests.** In the container, before anything is pushed: `--version`
must print the release version, `--help` must list its subcommands, and a
small synthetic dataset must segment via `packaging/smoke_test.sh` bind-
mounted into the container. Pushing happens only after all of that passes.

**Dry runs and subsets.** With `dry_run=true` the image is built and
smoke-tested but not pushed, and no release is consulted (no `latest`
decision). A `platforms` subset without `linux-x86_64` skips the job.

**Which ref provides the image files.** The job checks out the commit the
workflow runs on — for `release` events the tag (so tags must contain
`packaging/`, including `packaging/docker/` and
`packaging/is_latest_release.py`), for `workflow_dispatch` the branch the
workflow definition came from, which may be newer than the tag being rebuilt.

To build the same image locally:

```bash
gh release download vX.Y.Z -R <repo> -p 'baysor-X.Y.Z-linux-x86_64.tar.gz'
mkdir ctx && mv baysor-X.Y.Z-linux-x86_64.tar.gz ctx/
docker build --build-arg VERSION=X.Y.Z -f packaging/docker/Dockerfile -t baysor:X.Y.Z ctx
docker run --rm baysor:X.Y.Z --version
```

## Reproducing the Linux build

The Linux archive is built by exactly the command the workflow runs; it needs
only Docker:

```bash
packaging/linux/build-in-docker.sh
packaging/linux/test-in-docker.sh dist/baysor-X.Y.Z-linux-x86_64.tar.gz
```

- `build-in-docker.sh` builds the image from `packaging/linux/Dockerfile`
  (pinned `manylinux_2_28`, i.e. AlmaLinux 8, glibc 2.28, GCC 14) and runs
  `packaging/build_release.py --platform linux-x86_64` in it as the calling
  user. The result is `dist/baysor-X.Y.Z-linux-x86_64.tar.gz`. vcpkg, its
  downloads and binary cache are kept in `.release-cache/` (override with
  `BAYSOR_RELEASE_CACHE`), so only the first build compiles the dependencies.
  Parallelism is `BAYSOR_JOBS` (default 8). Extra arguments are passed to
  `build_release.py` (e.g. `--out DIR`).
- `test-in-docker.sh` runs the smoke test natively, in bare `almalinux:8` and
  `debian:10` containers, and under qemu with the `qemu64` CPU model.

`build_release.py` also runs `packaging/linux/check_binary.py`, which fails if
the binary needs a shared library other than glibc's, a `GLIBC_` symbol
version newer than 2.28, any libstdc++/libgcc symbol version, or an x86-64 ISA
level above the baseline.

macOS and Windows archives can be built on those systems with
`python3 packaging/build_release.py --platform macos-arm64` (needs Xcode
command line tools, CMake, Ninja, autoconf, autoconf-archive, automake,
libtool, pkg-config and bison — Homebrew's, since Apple's `/usr/bin/bison` 2.3
is too old for thrift, an Arrow/Parquet dependency) or
`python packaging/build_release.py --platform windows-x86_64` (needs Visual
Studio 2022 with the C++ workload and CMake).

## Portability of the binaries

| Platform | Requirement | How it is achieved |
| --- | --- | --- |
| Linux x86_64 | glibc ≥ 2.28 (RHEL/Alma/Rocky 8, Debian 10, Ubuntu 18.10 and newer), any x86-64 CPU | built in manylinux_2_28; libstdc++ and libgcc linked statically; all other libraries static from vcpkg |
| macOS arm64 | macOS ≥ 12 on Apple silicon | deployment target 12.0; static vcpkg libraries |
| Windows x64 | Windows 10 or newer, any x64 CPU | MSVC (SSE2 baseline); vcpkg and MSVC runtime DLLs shipped next to `baysor.exe` |

No part of the build uses `-march=native` or assumes SSE4/AVX. Everything is
compiled for the architecture baseline (the triplets in
`packaging/vcpkg-triplets/` set `-march=x86-64 -mtune=generic` for Linux and
build Arrow with `ARROW_SIMD_LEVEL=NONE`); faster code paths in Arrow, zstd,
OpenSSL, libjpeg-turbo and GMP (built with `--enable-fat`) are selected at run
time from the CPU's capabilities. CUDA is off in release builds.

To change a platform's floor, keep these in sync: the Dockerfile base image
and `GLIBC_FLOOR` in `check_binary.py` (Linux); `MACOS_DEPLOYMENT_TARGET` in
`build_release.py`, `VCPKG_OSX_DEPLOYMENT_TARGET` in
`packaging/vcpkg-triplets/arm64-osx-release.cmake` and
`CMAKE_OSX_DEPLOYMENT_TARGET` in the `release-macos-arm64` preset (macOS).
