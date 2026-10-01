# syntax=docker/dockerfile:1

# Runtime image for a Baysor release: the portable Linux x86_64 release binary
# (needs only glibc >= 2.28; Debian 12 ships 2.36) on a slim base, nothing else.
# The build context holds baysor-<VERSION>-linux-x86_64.tar.gz. Built by the
# `docker` job of .github/workflows/release.yml; for a local build see
# RELEASING.md, "Docker images". SOURCE_URL is the repository for the OCI
# source label.
FROM debian:12-slim

ARG VERSION
ARG SOURCE_URL=https://github.com/kharchenkolab/Baysor

LABEL org.opencontainers.image.source="${SOURCE_URL}" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.description="Baysor: spatial single-cell segmentation from transcriptomics data"

# ADD unpacks the archive; the final `test` fails the build if the binary does
# not run in this image or prints another version.
ADD baysor-${VERSION}-linux-x86_64.tar.gz /tmp/baysor-archive/
RUN set -eu; \
    tree="/tmp/baysor-archive/baysor-${VERSION}-linux-x86_64"; \
    install -m 0755 "$tree/bin/baysor" /usr/local/bin/baysor; \
    mkdir -p /usr/share/doc/baysor; \
    install -m 0644 "$tree/LICENSE" "$tree/README.md" /usr/share/doc/baysor/; \
    rm -rf /tmp/baysor-archive; \
    test "$(/usr/local/bin/baysor --version)" = "baysor ${VERSION}"

# Non-root by default (see docs/installation.md): files written to a
# bind-mounted /data belong to UID/GID 1000 unless `docker run --user` is used.
RUN groupadd --gid 1000 baysor \
    && useradd --uid 1000 --gid 1000 --home-dir /data --no-create-home \
        --shell /usr/sbin/nologin baysor \
    && mkdir -p /data \
    && chown baysor:baysor /data

USER baysor
WORKDIR /data

ENTRYPOINT ["/usr/local/bin/baysor"]
