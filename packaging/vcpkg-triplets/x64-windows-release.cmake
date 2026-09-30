# Release triplet for the Windows x64 binary.
#
# DLLs with the dynamic MSVC runtime (the same linkage as the x64-windows
# triplet used by CI), release configuration only. The DLLs and the MSVC/OpenMP
# runtime are copied next to baysor.exe at install time. MSVC targets SSE2 by
# default (no /arch flag).
#
# vcpkg hashes this file into every package ABI, so keep the settings inline.
set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE dynamic)
set(VCPKG_BUILD_TYPE release)

if(PORT STREQUAL "arrow")
    # Arrow's default ARROW_SIMD_LEVEL on x64 assumes SSE4.2. Keep only the
    # runtime-dispatched kernels (ARROW_RUNTIME_SIMD_LEVEL).
    list(APPEND VCPKG_CMAKE_CONFIGURE_OPTIONS -DARROW_SIMD_LEVEL=NONE -DARROW_RUNTIME_SIMD_LEVEL=MAX)
endif()
