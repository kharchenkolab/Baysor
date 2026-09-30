# syntax=docker/dockerfile:1

# Runtime image for a Baysor release: installs the portable Linux x86_64
# release binary (baysor-X.Y.Z-linux-x86_64.tar.gz, needs only glibc >= 2.28;
# Debian 12 ships glibc 2.36) — no compilers, no source build, nothing but
# the base image and the three files of the archive.
#
# Built by the `docker` job of .github/workflows/release.yml from the release
# archive artifact; see RELEASING.md, "Docker images". To build locally:
#
#   gh release download vX.Y.Z -R <repo> -p 'baysor-X.Y.Z-linux-x86_64.tar.gz'
#   mkdir ctx && mv baysor-X.Y.Z-linux-x86_64.tar.gz ctx/
#   docker build --build-arg VERSION=X.Y.Z \
#     -f packaging/docker/Dockerfile -t baysor:X.Y.Z ctx
#
# Build arguments:
#   VERSION     release version X.Y.Z; must match the archive's directory
#               name (the build verifies `baysor --version`)
#   SOURCE_URL  repository the image is built from (OCI .source label)
#   BASE_IMAGE  runtime base image, pinned by tag; override to rebuild
#               against another Debian/Ubuntu release

ARG BASE_IMAGE=debian:12-slim
FROM ${BASE_IMAGE}

ARG VERSION
ARG SOURCE_URL=https://github.com/kharchenkolab/Baysor

LABEL org.opencontainers.image.source="${SOURCE_URL}" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.description="Baysor: spatial single-cell segmentation from transcriptomics data"

# ADD unpacks the local tarball into baysor-X.Y.Z-linux-x86_64/; install the
# tree (bin/baysor -> /usr/local/bin, LICENSE/README.md -> /usr/share/doc)
# and drop the archive, all in one layer. The final `test` fails the build
# fast if the binary does not run in this image or prints another version.
ADD baysor-${VERSION}-linux-x86_64.tar.gz /tmp/baysor-archive/
RUN set -eu; \
    tree="/tmp/baysor-archive/baysor-${VERSION}-linux-x86_64"; \
    install -m 0755 "$tree/bin/baysor" /usr/local/bin/baysor; \
    mkdir -p /usr/share/doc/baysor; \
    install -m 0644 "$tree/LICENSE" "$tree/README.md" /usr/share/doc/baysor/; \
    rm -rf /tmp/baysor-archive; \
    test "$(/usr/local/bin/baysor --version)" = "baysor ${VERSION}"

# Non-root by default: baysor needs no privileges, and the user's home
# doubles as the working directory. The volume trade-off is documented in
# docs/installation.md: files written to a bind-mounted /data are owned by
# UID/GID 1000 on the host; override with --user "$(id -u):$(id -g)".
RUN groupadd --gid 1000 baysor \
    && useradd --uid 1000 --gid 1000 --home-dir /data --no-create-home \
        --shell /usr/sbin/nologin baysor \
    && mkdir -p /data \
    && chown baysor:baysor /data

USER baysor
WORKDIR /data

ENTRYPOINT ["/usr/local/bin/baysor"]
