cmake_minimum_required(VERSION 3.20)
project(baysor VERSION 0.8.3 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Keep default installs user-owned. Presets and explicit command-line values can
# still override this with CMAKE_INSTALL_PREFIX or cmake --install --prefix.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install" CACHE PATH "Default Baysor install prefix." FORCE)
endif()

option(BAYSOR_EXPORT_COMPILE_COMMANDS "Export compile_commands.json for developer tooling" OFF)
if(BAYSOR_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()
include(GNUInstallDirs)

if(APPLE)
    foreach(_homebrew_prefix /opt/homebrew /usr/local)
        if(EXISTS "${_homebrew_prefix}")
            list(APPEND CMAKE_PREFIX_PATH "${_homebrew_prefix}")
        endif()
    endforeach()
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(BaysorDependencies)

# OpenMP is no longer used: all parallelism goes through Baysor's own
# std::thread pool (include/baysor/utils/thread_pool.h). The pool also backs
# the subpar-based FetchContent dependencies via their custom-parallel hooks.

# Prefer optimized single-config builds unless the user explicitly chose
# something else. This does not affect multi-config generators.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Options
option(BAYSOR_WITH_CUDA "Enable CUDA/GPU support" OFF)
option(BAYSOR_WITH_TESTS "Build tests" OFF)
option(BAYSOR_WITH_REPORTING "Build reporting/visualization support" ON)
option(BAYSOR_WITH_COVERAGE "Instrument Baysor targets for gcov coverage" OFF)

if(BAYSOR_WITH_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR
            "BAYSOR_WITH_COVERAGE requires GCC or Clang to produce gcov data, "
            "but the compiler is ${CMAKE_CXX_COMPILER_ID}. "
            "Reconfigure with a GCC/Clang toolchain or -DBAYSOR_WITH_COVERAGE=OFF.")
    endif()
endif()

# Apply gcov instrumentation to one of the project's own targets. Only called
# for baysor_lib, baysor_umap_wrappers, baysor and baysor_tests so that
# third-party and FetchContent code stays uninstrumented. No-op when coverage
# is off, so default builds are unaffected.
function(baysor_apply_coverage target)
    if(NOT BAYSOR_WITH_COVERAGE)
        return()
    endif()
    target_compile_options(${target} PRIVATE --coverage -O0 -g)
    target_link_options(${target} PRIVATE --coverage)
endfunction()

# Dependencies managed by vcpkg or find_package
# Core (always required)
baysor_find_package(Eigen3
    VERSION 3.3
    TARGET Eigen3::Eigen
    INSTALL_HINT "Ubuntu: sudo apt-get install libeigen3-dev\nmacOS: brew install eigen\nvcpkg: vcpkg install eigen3"
)
baysor_find_package(spdlog
    TARGET spdlog::spdlog
    INSTALL_HINT "Ubuntu: sudo apt-get install libspdlog-dev\nmacOS: brew install spdlog\nvcpkg: vcpkg install spdlog"
)

# UMAP and dependencies (all header-only, fetched via FetchContent).
# We only need the headers, so avoid adding upstream CMake projects here.
# Some upstream CMakeLists require newer CMake versions than Baysor itself.
include(FetchContent)
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()
FetchContent_Declare(aarand
    GIT_REPOSITORY https://github.com/LTLA/aarand
    GIT_TAG v1.0.2
)
# knncolle v2.3.0 transitively depends on kmeans v3.x and subpar v0.2.1
FetchContent_Declare(kmeans
    GIT_REPOSITORY https://github.com/LTLA/CppKmeans
    GIT_TAG v3.1.1
)
FetchContent_Declare(subpar
    GIT_REPOSITORY https://github.com/LTLA/subpar
    GIT_TAG v0.3.1
)
FetchContent_Declare(knncolle
    GIT_REPOSITORY https://github.com/knncolle/knncolle
    GIT_TAG v2.3.0
)
FetchContent_Declare(irlba
    GIT_REPOSITORY https://github.com/LTLA/CppIrlba
    GIT_TAG v2.0.2
)
FetchContent_Declare(umappp
    GIT_REPOSITORY https://github.com/LTLA/umappp
    GIT_TAG v2.0.1
)

function(baysor_fetch_header_only_dependency dependency target alias)
    FetchContent_GetProperties(${dependency})
    if(NOT ${dependency}_POPULATED)
        FetchContent_Populate(${dependency})
    endif()

    if(NOT TARGET ${target})
        add_library(${target} INTERFACE)
        target_include_directories(${target} INTERFACE "${${dependency}_SOURCE_DIR}/include")
        target_compile_features(${target} INTERFACE cxx_std_17)
    endif()
    if(alias AND NOT TARGET ${alias})
        add_library(${alias} ALIAS ${target})
    endif()
endfunction()

baysor_fetch_header_only_dependency(aarand aarand ltla::aarand)
baysor_fetch_header_only_dependency(subpar subpar ltla::subpar)
baysor_fetch_header_only_dependency(kmeans kmeans ltla::kmeans)
baysor_fetch_header_only_dependency(knncolle knncolle knncolle::knncolle)
baysor_fetch_header_only_dependency(irlba irlba ltla::irlba)
baysor_fetch_header_only_dependency(umappp umappp libscran::umappp)

target_link_libraries(kmeans INTERFACE ltla::aarand ltla::subpar)
target_link_libraries(knncolle INTERFACE ltla::kmeans ltla::subpar)
target_link_libraries(irlba INTERFACE Eigen3::Eigen ltla::aarand ltla::subpar)
target_link_libraries(umappp INTERFACE ltla::aarand ltla::irlba knncolle::knncolle)

# Patch irlba v2.0.2 bug: Eigen::VectorXd has no begin()/end().
# Replace std::max_element(BS.begin(), BS.end()) with BS.maxCoeff().
FetchContent_GetProperties(irlba)
set(_irlba_compute "${irlba_SOURCE_DIR}/include/irlba/compute.hpp")
if(EXISTS "${_irlba_compute}")
    file(READ "${_irlba_compute}" _irlba_content)
    string(REPLACE
        "*std::max_element(BS.begin(), BS.end())"
        "BS.maxCoeff()"
        _irlba_content "${_irlba_content}")
    file(WRITE "${_irlba_compute}" "${_irlba_content}")
endif()

# Spatial
baysor_find_package(CGAL
    TARGET CGAL::CGAL
    INSTALL_HINT "Ubuntu: sudo apt-get install libcgal-dev\nmacOS: brew install cgal\nvcpkg: vcpkg install cgal"
)
# nanoflann is header-only, vendored in third_party/ or found via find_package
# find_package(nanoflann REQUIRED)

# I/O
baysor_find_package(Parquet
    INSTALL_HINT "Ubuntu: install libparquet-dev from the Apache Arrow apt repository\nmacOS: brew install apache-arrow\nvcpkg: vcpkg install 'arrow[parquet]'"
)
if(NOT TARGET Arrow::arrow_shared AND NOT TARGET Arrow::arrow_static)
    baysor_find_package(Arrow
        INSTALL_HINT "Ubuntu: install libarrow-dev from the Apache Arrow apt repository\nmacOS: brew install apache-arrow\nvcpkg: vcpkg install 'arrow[compute,csv,parquet]'"
    )
endif()
baysor_find_package(HDF5
    COMPONENTS C
    INSTALL_HINT "Ubuntu: sudo apt-get install libhdf5-dev\nmacOS: brew install hdf5\nvcpkg: vcpkg install 'hdf5[zlib]'"
)
baysor_find_package(nlohmann_json
    TARGET nlohmann_json::nlohmann_json
    INSTALL_HINT "Ubuntu: sudo apt-get install nlohmann-json3-dev\nmacOS: brew install nlohmann-json\nvcpkg: vcpkg install nlohmann-json"
)
# zlib: PNG encoding of the HTML report images (already a dependency of
# HDF5, TIFF and Arrow/Parquet).
baysor_find_package(ZLIB
    TARGET ZLIB::ZLIB
    INSTALL_HINT "Ubuntu: sudo apt-get install zlib1g-dev\nmacOS: part of the SDK (or brew install zlib)\nvcpkg: vcpkg install zlib"
)
baysor_find_package(TIFF
    TARGET TIFF::TIFF
    INSTALL_HINT "Ubuntu: sudo apt-get install libtiff-dev\nmacOS: brew install libtiff\nvcpkg: vcpkg install tiff"
)
baysor_select_target(BAYSOR_ARROW_TARGET Arrow::arrow_shared Arrow::arrow_static)
baysor_select_target(BAYSOR_PARQUET_TARGET Parquet::parquet_shared Parquet::parquet_static)
message(STATUS "Using Arrow target: ${BAYSOR_ARROW_TARGET}")
message(STATUS "Using Parquet target: ${BAYSOR_PARQUET_TARGET}")
if(TARGET hdf5::hdf5)
    set(BAYSOR_HDF5_LIBRARIES hdf5::hdf5)
elseif(TARGET HDF5::HDF5)
    set(BAYSOR_HDF5_LIBRARIES HDF5::HDF5)
else()
    set(BAYSOR_HDF5_LIBRARIES ${HDF5_C_LIBRARIES})
    set(BAYSOR_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
endif()

# CLI
# CLI11 and toml++ are header-only
# find_package(CLI11 REQUIRED)
# find_package(tomlplusplus REQUIRED)

if(BAYSOR_WITH_CUDA)
    enable_language(CUDA)
    baysor_find_package(CUDAToolkit
        TARGET CUDA::cudart
        INSTALL_HINT "Install the NVIDIA CUDA Toolkit and reconfigure with -DBAYSOR_WITH_CUDA=ON only on CUDA-capable systems."
    )
    # find_package(faiss REQUIRED)
    add_compile_definitions(BAYSOR_CUDA_ENABLED)
endif()

# Threads: Baysor's own std::thread pool (and the thread pools of the
# FetchContent dependencies routed through it) need pthreads on POSIX.
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Eigen dense products used to be parallelized by OpenMP. Without OpenMP we
# enable Eigen's own GEMM thread pool when the Eigen in use supports it
# (EIGEN_GEMM_THREADPOOL, Eigen >= 3.4.90/5.x); otherwise dense products fall
# back to single-threaded GEMM. The macro must be defined consistently for
# every translation unit that includes Eigen, so it is added globally here.
set(_eigen_gemm_threadpool OFF)
get_target_property(_eigen_includes Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
foreach(_dir IN LISTS _eigen_includes)
    if(EXISTS "${_dir}/Eigen/src/Core/products/Parallelizer.h")
        file(STRINGS "${_dir}/Eigen/src/Core/products/Parallelizer.h"
             _eigen_tp REGEX "EIGEN_GEMM_THREADPOOL")
        if(_eigen_tp)
            set(_eigen_gemm_threadpool ON)
        endif()
    endif()
endforeach()
if(_eigen_gemm_threadpool)
    add_compile_definitions(EIGEN_GEMM_THREADPOOL)
    message(STATUS "Eigen GEMM: parallelized through Eigen's own thread pool (EIGEN_GEMM_THREADPOOL)")
else()
    message(STATUS "Eigen GEMM: this Eigen version has no thread-pool GEMM; dense products run single-threaded")
endif()

# --- Library target ---
add_library(baysor_lib
    # Utils
    src/utils/general.cpp
    src/utils/options.cpp
    src/utils/thread_pool.cpp
    src/utils/xenium.cpp

    # Data loading
    src/data_loading/data.cpp
    src/data_loading/prior_segmentation.cpp

    # Processing - distributions
    src/processing/distributions/mv_normal.cpp
    src/processing/distributions/categorical_smoothed.cpp

    # Processing - models
    src/processing/models/adj_list.cpp
    src/processing/models/component.cpp
    src/processing/models/bmm_data.cpp

    # Processing - data_processing
    src/processing/data_processing/triangulation.cpp
    src/processing/data_processing/initialization.cpp
    src/processing/data_processing/neighborhood_composition.cpp
    src/processing/data_processing/noise_estimation.cpp
    src/processing/data_processing/boundary_estimation.cpp
    # Processing - bmm_algorithm
    src/processing/bmm_algorithm/bmm_algorithm.cpp
    src/processing/bmm_algorithm/molecule_clustering.cpp
    src/processing/bmm_algorithm/molecule_clustering_ica.cpp
    src/processing/bmm_algorithm/molecule_clustering_louvain.cpp
    src/processing/bmm_algorithm/compartment_segmentation.cpp
    src/processing/bmm_algorithm/tracing.cpp
    src/processing/bmm_algorithm/history_analysis.cpp

    # Processing - utils
    src/processing/utils/utils.cpp
    src/processing/utils/convex_hull.cpp

    # Reporting
    src/reporting/color_utils.cpp
    src/reporting/gene_structure.cpp
    src/reporting/output.cpp
    src/reporting/preview_report.cpp
    src/reporting/run_report.cpp
)

add_library(baysor_umap_wrappers STATIC
    src/processing/data_processing/umap_wrappers.cpp
)

target_include_directories(baysor_umap_wrappers PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(baysor_umap_wrappers PUBLIC
    Eigen3::Eigen
    umappp
    knncolle::knncolle
)
baysor_apply_coverage(baysor_umap_wrappers)

target_include_directories(baysor_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${BAYSOR_HDF5_INCLUDE_DIRS}
)

target_link_libraries(baysor_lib PUBLIC
    Eigen3::Eigen
    Threads::Threads
    spdlog::spdlog
    CGAL::CGAL
    ${BAYSOR_ARROW_TARGET}
    ${BAYSOR_PARQUET_TARGET}
    ${BAYSOR_HDF5_LIBRARIES}
    nlohmann_json::nlohmann_json
    TIFF::TIFF
    ZLIB::ZLIB
    baysor_umap_wrappers
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(baysor_lib PRIVATE -Wno-stringop-overflow)
endif()
baysor_apply_coverage(baysor_lib)

if(BAYSOR_WITH_CUDA)
    target_link_libraries(baysor_lib PUBLIC
        CUDA::cudart
        CUDA::cusparse
    )
endif()

# --- Executable ---
add_executable(baysor
    src/cli/main.cpp
)

target_link_libraries(baysor PRIVATE baysor_lib)
target_compile_definitions(baysor PRIVATE BAYSOR_VERSION="${PROJECT_VERSION}")
install(TARGETS baysor RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Release binaries (packaging/) ship the compiler runtime with the executable.
# On MSVC this installs the VC++ runtime DLLs next to baysor.exe;
# vcpkg's X_VCPKG_APPLOCAL_DEPS_INSTALL copies the dependency DLLs. Linux and
# macOS release builds link their runtimes statically, so nothing to do there.
option(BAYSOR_INSTALL_RUNTIME "Install the compiler runtime libraries next to the executable" OFF)
if(BAYSOR_INSTALL_RUNTIME AND MSVC)
    set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ${CMAKE_INSTALL_BINDIR})
    include(InstallRequiredSystemLibraries)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(baysor PRIVATE -Wno-stringop-overflow)
endif()
baysor_apply_coverage(baysor)

# --- Tests ---
if(BAYSOR_WITH_TESTS)
    enable_testing()
    baysor_find_package(GTest
        TARGET GTest::gtest_main
        INSTALL_HINT "Ubuntu: sudo apt-get install libgtest-dev\nmacOS: brew install googletest\nvcpkg: enable the 'tests' manifest feature or install gtest"
    )

    # Every tests/test_*.cpp file is compiled into the single test binary, so
    # new test files are picked up without editing CMake (CONFIGURE_DEPENDS
    # re-globs on build).
    file(GLOB BAYSOR_TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.cpp"
    )
    if(NOT BAYSOR_TEST_SOURCES)
        message(FATAL_ERROR "BAYSOR_WITH_TESTS is ON but no tests/test_*.cpp files were found.")
    endif()

    add_executable(baysor_tests
        ${BAYSOR_TEST_SOURCES}
    )
    target_link_libraries(baysor_tests PRIVATE baysor_lib GTest::gtest_main)
    # End-to-end CLI tests spawn the instrumented `baysor` executable as a
    # subprocess (its own .gcda counters cover src/cli/main.cpp, which is not
    # linked into the test binary). Build it first and pass its path in.
    add_dependencies(baysor_tests baysor)
    target_compile_definitions(baysor_tests PRIVATE
        BAYSOR_CLI_PATH="$<TARGET_FILE:baysor>"
        BAYSOR_VERSION="${PROJECT_VERSION}"
    )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(baysor_tests PRIVATE -Wno-stringop-overflow)
    endif()
    baysor_apply_coverage(baysor_tests)
    add_test(NAME baysor_tests COMMAND baysor_tests)
    # The suite spawns the CLI as a subprocess; a hung `baysor` must not
    # stall CI forever (review COV-6).
    set_tests_properties(baysor_tests PROPERTIES TIMEOUT 900)
endif()

# --- Coverage report ---
# 'cmake --build <dir> --target coverage' deletes stale .gcda counters, reruns
# the tests, and writes text/JSON/HTML reports under <build-dir>/coverage/.
if(BAYSOR_WITH_COVERAGE)
    find_program(BAYSOR_GCOVR gcovr)
    if(BAYSOR_GCOVR)
        # Use the gcov that matches the compiler. For conda GCC that is the
        # triple-prefixed gcov sitting next to the compiler, e.g.
        # x86_64-conda-linux-gnu-g++ -> x86_64-conda-linux-gnu-gcov.
        get_filename_component(_baysor_cxx_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
        get_filename_component(_baysor_cxx_name "${CMAKE_CXX_COMPILER}" NAME)
        string(REGEX REPLACE "g\\+\\+$" "gcov" _baysor_gcov_name "${_baysor_cxx_name}")
        if(NOT _baysor_gcov_name STREQUAL _baysor_cxx_name
                AND EXISTS "${_baysor_cxx_dir}/${_baysor_gcov_name}")
            set(BAYSOR_GCOV "${_baysor_cxx_dir}/${_baysor_gcov_name}")
        else()
            # Fall back to a gcov on PATH; if none is configured gcovr
            # resolves plain 'gcov' itself.
            find_program(BAYSOR_GCOV NAMES gcov)
            if(NOT BAYSOR_GCOV)
                set(BAYSOR_GCOV gcov)
            endif()
        endif()

        add_custom_target(coverage
            COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_coverage.sh"
                    "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}"
                    "${CMAKE_CTEST_COMMAND}"
                    "${BAYSOR_GCOVR}"
                    "${BAYSOR_GCOV}"
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            COMMENT "Running tests and generating gcov coverage reports"
            VERBATIM
        )
        if(TARGET baysor_tests)
            add_dependencies(coverage baysor_tests)
        endif()
    else()
        message(STATUS "gcovr not found: the 'coverage' target will not be available")
    endif()
endif()
