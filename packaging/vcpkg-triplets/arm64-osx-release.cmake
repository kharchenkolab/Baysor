# Release triplet for the macOS arm64 binary.
#
# Static libraries, release configuration only. The deployment target must
# match CMAKE_OSX_DEPLOYMENT_TARGET in the release-macos-arm64 CMake preset.
# The default arm64 code generation (NEON, Apple M1 baseline) runs on every
# Apple silicon Mac.
#
# vcpkg hashes this file into every package ABI, so keep the settings inline.
set(VCPKG_TARGET_ARCHITECTURE arm64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)
set(VCPKG_CMAKE_SYSTEM_NAME Darwin)
set(VCPKG_OSX_ARCHITECTURES arm64)
set(VCPKG_OSX_DEPLOYMENT_TARGET 12.0)
set(VCPKG_BUILD_TYPE release)
