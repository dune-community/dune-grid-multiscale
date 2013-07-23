# This file is part of the dune-grid-multiscale project:
#   http://users.dune-project.org/projects/dune-grid-multiscale
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

project( dune-grid-multiscale )

# cmake specific
cmake_minimum_required( VERSION 2.8 )
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
set( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
set( CMAKE_COLOR_MAKEFILE ON)
set( CMAKE_VERBOSE_MAKEFILE ON )

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

# includes from dune-stuff
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)
include(GridUtils)

# needs to come first to include correct config.h
include_dir(${dune-grid-multiscale_BINARY_DIR})
include_dir(${dune-grid-multiscale_SOURCE_DIR})

# options, that can be changed by the user in order to customise the module
set( ENABLE_MPI
     "0" CACHE STRING
     "set to 1 to enable mpi" )

SET( CUSTOM_FLAGS
     "-Wall -Wextra -O0 -DDEBUG -g -ggdb -fdiagnostics-show-option" CACHE LIST
     "CUSTOM_FLAGS" )

add_definitions( ${CUSTOM_FLAGS} )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

add_dune_modules(
  grid
  stuff
  common
  geometry
)

# header
file( GLOB_RECURSE grid-multiscale "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )
set( COMMON_HEADER ${grid-multiscale} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

#libs
set( COMMON_LIBS
     ${DUNE_LIBS}
     ${ParaLIBS}
     "boost_filesystem" "boost_system" 
     ${GRIDLIBS}
)

add_subdirectory( doc )
add_subdirectory( examples )
add_subdirectory( dune )

