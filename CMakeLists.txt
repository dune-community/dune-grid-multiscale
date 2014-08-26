# This file is part of the dune-grid-multiscale project:
#   http://users.dune-project.org/projects/dune-grid-multiscale
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# set up project
project("dune-grid-multiscale" C CXX)

if(CMAKE_WITH_AUTOTOOLS)

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
  set(CMAKE_MODULE_PATH
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake/modules)
  include(DuneUtils)
  include(GridUtils)

  # needs to come first to include correct config.h
  include(${dune-grid-multiscale_BINARY_DIR})
  include(${dune-grid-multiscale_SOURCE_DIR})

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
          grid-multiscale
          fem
          pdelab
          typetree
          localfunctions
          istl
    grid
    stuff
    common
    geometry
  )

  # header
  file( GLOB_RECURSE grid-multiscale "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" )
  set( COMMON_HEADER ${grid-multiscale} ${DUNE_HEADERS} )
  set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

  #libs
  set( COMMON_LIBS
       ${DUNE_LIBS}
       ${ParaLIBS}
       "boost_filesystem" "boost_system"
       ${GRIDLIBS}
  )

  add_subdirectory(test)

else(CMAKE_WITH_AUTOTOOLS)

  #circumvent not building docs
  set(BUILD_DOCS 1)

  # general stuff
  cmake_minimum_required(VERSION 2.8.6)

  #find dune-common and set the module path
  find_package(dune-common)
  list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
       "${CMAKE_SOURCE_DIR}/cmake/modules")# make sure our own modules are found

  #include the dune macros
  include(DuneMacros)

  # start a dune project with information from dune.module
  dune_project()

  include(DuneUtils)

  # pass most important options when "make distcheck" is used
  set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT}  CXX="${CXX}" CC="${CC}")
  # include not needed for CMake
  # include $(top_srcdir)/am/top-rules
  # include not needed for CMake
  # include $(top_srcdir)/am/global-rules
  #foreach(i dune doc)
  #  if(${i} STREQUAL "test")
  #    set(opt EXCLUDE_FROM_ALL)
  #  endif(${i} STREQUAL "test")
  #  add_subdirectory(${i} ${opt})
  #  unset(opt)
  #endforeach(i ${SUBDIRS})

  # finalize the dune project, e.g., generate config.h etc.
  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

endif(CMAKE_WITH_AUTOTOOLS)
