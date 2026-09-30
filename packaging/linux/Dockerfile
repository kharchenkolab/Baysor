# Build environment for the portable Linux x86_64 release binary.
#
# manylinux_2_28 is AlmaLinux 8 (glibc 2.28) with GCC 14 from gcc-toolset.
# Linking against its glibc sets the binary's glibc floor to 2.28; the
# toolset ships static libstdc++, libgcc and libgomp, so nothing else from the
# build machine is needed at run time. Pinned by date tag for reproducibility.
FROM quay.io/pypa/manylinux_2_28_x86_64:2026.09.30-1

# Host tools the vcpkg ports need: nasm (libjpeg-turbo's runtime-dispatched
# SIMD), flex/bison (thrift), perl modules and kernel headers (OpenSSL),
# autoconf-archive (autotools ports: gmp, mpfr, libiconv), ninja and zip.
RUN dnf install -y --enablerepo=powertools \
        autoconf-archive bison flex kernel-headers nasm ninja-build \
        perl-IPC-Cmd perl-Time-Piece zip \
    && dnf clean all \
    && rm -rf /var/cache/dnf
