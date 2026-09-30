# Release triplet for the portable Linux x86_64 binary.
#
# Static libraries, release configuration only, compiled for the generic
# x86-64 baseline (SSE2) so the result runs on any x86_64 CPU. Code that needs
# newer instruction sets must select it at run time via cpuid (Arrow, zstd,
# OpenSSL and libjpeg-turbo do).
#
# vcpkg hashes this file into every package ABI, so keep the settings inline
# rather than include()-ing a shared file.
set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)
set(VCPKG_CMAKE_SYSTEM_NAME Linux)
set(VCPKG_BUILD_TYPE release)

set(VCPKG_C_FLAGS "-march=x86-64 -mtune=generic")
set(VCPKG_CXX_FLAGS "-march=x86-64 -mtune=generic")

if(PORT STREQUAL "arrow")
    # Arrow's default ARROW_SIMD_LEVEL on x86_64 compiles all of Arrow with
    # SSE4.2. Keep only the runtime-dispatched kernels (ARROW_RUNTIME_SIMD_LEVEL).
    list(APPEND VCPKG_CMAKE_CONFIGURE_OPTIONS -DARROW_SIMD_LEVEL=NONE -DARROW_RUNTIME_SIMD_LEVEL=MAX)
endif()
