function(baysor_find_package package_name)
    set(options)
    set(one_value_args VERSION TARGET INSTALL_HINT)
    set(multi_value_args COMPONENTS)
    cmake_parse_arguments(BAYSOR_DEP "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    set(_find_args)
    if(BAYSOR_DEP_VERSION)
        list(APPEND _find_args "${BAYSOR_DEP_VERSION}")
    endif()
    list(APPEND _find_args QUIET)
    if(BAYSOR_DEP_COMPONENTS)
        list(APPEND _find_args COMPONENTS ${BAYSOR_DEP_COMPONENTS})
    endif()

    find_package(${package_name} ${_find_args})

    if(package_name STREQUAL "Eigen3" AND BAYSOR_DEP_TARGET AND NOT TARGET "${BAYSOR_DEP_TARGET}")
        set(_eigen_hints)
        foreach(_prefix IN LISTS CMAKE_PREFIX_PATH)
            list(APPEND _eigen_hints
                "${_prefix}"
                "${_prefix}/include"
                "${_prefix}/include/eigen3"
            )
        endforeach()
        if(DEFINED VCPKG_INSTALLED_DIR AND DEFINED VCPKG_TARGET_TRIPLET)
            list(APPEND _eigen_hints
                "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}"
                "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include"
                "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include/eigen3"
            )
        endif()
        list(APPEND _eigen_hints
            /opt/homebrew
            /opt/homebrew/include
            /opt/homebrew/include/eigen3
            /usr/local
            /usr/local/include
            /usr/local/include/eigen3
            /usr
            /usr/include
            /usr/include/eigen3
        )

        find_path(EIGEN3_INCLUDE_DIR
            NAMES Eigen/Core
            HINTS ${_eigen_hints}
            PATH_SUFFIXES eigen3 include/eigen3
        )
        if(EIGEN3_INCLUDE_DIR)
            add_library("${BAYSOR_DEP_TARGET}" INTERFACE IMPORTED)
            target_include_directories("${BAYSOR_DEP_TARGET}" INTERFACE "${EIGEN3_INCLUDE_DIR}")
            set(EIGEN3_FOUND TRUE)
            set(EIGEN3_VERSION_STRING "unknown")
            message(STATUS "Found Eigen3 headers: ${EIGEN3_INCLUDE_DIR}")
        endif()
    endif()

    set(_found FALSE)
    if(DEFINED ${package_name}_FOUND AND ${${package_name}_FOUND})
        set(_found TRUE)
    elseif(package_name STREQUAL "Eigen3" AND DEFINED EIGEN3_FOUND AND EIGEN3_FOUND)
        set(_found TRUE)
    elseif(BAYSOR_DEP_TARGET AND TARGET "${BAYSOR_DEP_TARGET}")
        set(_found TRUE)
    endif()

    if(NOT _found)
        message(FATAL_ERROR
            "Required dependency '${package_name}' was not found.\n"
            "Install it with your system package manager, or use the vcpkg preset documented in docs/installation.md.\n"
            "${BAYSOR_DEP_INSTALL_HINT}\n"
            "If it is already installed in a custom prefix, reconfigure with:\n"
            "  cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/prefix\n"
            "or set ${package_name}_DIR to the directory containing ${package_name}Config.cmake."
        )
    endif()

    if(BAYSOR_DEP_TARGET AND NOT TARGET "${BAYSOR_DEP_TARGET}")
        message(FATAL_ERROR
            "Dependency '${package_name}' was found, but CMake target '${BAYSOR_DEP_TARGET}' is missing.\n"
            "This usually means the installed package is too old or was built without CMake package metadata.\n"
            "${BAYSOR_DEP_INSTALL_HINT}"
        )
    endif()

    set(_version "unknown")
    set(_version_vars
        ${package_name}_VERSION
        ${package_name}_VERSION_STRING
        ${package_name}_VERSION_MAJOR)
    if(package_name STREQUAL "Eigen3")
        list(APPEND _version_vars EIGEN3_VERSION_STRING EIGEN3_VERSION)
    endif()

    foreach(_version_var ${_version_vars})
        if(DEFINED ${_version_var} AND NOT "${${_version_var}}" STREQUAL "")
            set(_version "${${_version_var}}")
            break()
        endif()
    endforeach()

    if(BAYSOR_DEP_TARGET)
        message(STATUS "Found ${package_name}: version ${_version}, target ${BAYSOR_DEP_TARGET}")
    else()
        message(STATUS "Found ${package_name}: version ${_version}")
    endif()
endfunction()

function(baysor_select_target out_var)
    foreach(_target IN LISTS ARGN)
        if(TARGET "${_target}")
            set(${out_var} "${_target}" PARENT_SCOPE)
            return()
        endif()
    endforeach()

    message(FATAL_ERROR
        "None of the expected CMake targets exists: ${ARGN}\n"
        "The dependency was found, but its CMake package does not expose a supported target."
    )
endfunction()
