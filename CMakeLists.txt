# set up project
project("dune-grid-multiscale" CXX)

cmake_minimum_required(VERSION 3.1)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

# dune-common
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH
     "${PROJECT_SOURCE_DIR}/cmake/modules"
     "${dune-common_MODULE_PATH}")
include(DuneMacros)

find_package(dune-xt-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH
     "${dune-xt-common_MODULE_PATH}")
include(DuneUtils)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages(INCLUDE_DIRS ${dune-grid-multiscale_SOURCE_DIR}/dune
                         VERBOSE)


add_header_listing()
make_dependent_modules_sys_included()

add_subdirectory(dune)
add_subdirectory(cmake/modules)

# enable headercheck
add_definitions("-DENABLE_HEADERCHECK=1")
add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})

# finalize the dune project, e.g., generate config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
