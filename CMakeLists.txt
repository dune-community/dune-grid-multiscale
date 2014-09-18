# This file is part of the dune-grid-multiscale project:
#   http://users.dune-project.org/projects/dune-grid-multiscale
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

if(NOT CMAKE_BUILD_TYPE)
  message( STATUS "No build type selected, setting default build type 'Debug'" )
  set( CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "DEBUG/RELEASE")
endif(NOT CMAKE_BUILD_TYPE)

# set up project
project("dune-grid-multiscale" C CXX)

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

# cmake specific
cmake_minimum_required( VERSION 2.8 )

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

  set(CMAKE_COLOR_MAKEFILE ON)

  # guess dune-common build dir
  if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR ${PROJECT_BINARY_DIR})
  endif()

  #find dune-common and set the module path
  find_package(dune-common REQUIRED)
  list(APPEND CMAKE_MODULE_PATH
       ${dune-common_MODULE_PATH}
       "${CMAKE_SOURCE_DIR}/../dune-stuff/cmake")

  #include the dune macros
  include(DuneMacros)

  # start a dune project with information from dune.module
  dune_project()

  include(DuneUtils)
  
  link_directories("/usr/lib" ${CCGNU_LIBRARY_DIRS})
  
    # header
  file(GLOB_RECURSE grid-multiscale "${CMAKE_CURRENT_SOURCE_DIR}/*.hh")
  set(COMMON_HEADER ${grid-multiscale} ${DUNE_HEADERS})
  set_source_files_properties(${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1)

  set(GRIDDIM "2" CACHE STRING "GRIDDIM")
  set(GRIDTYPE "ALUGRID_CONFORM" CACHE STRING "GRIDTYPE")
  add_definitions("-DGRIDDIM=${GRIDDIM}" "-DGRIDTYPE=${GRIDTYPE}")
  if (ALUGRID_FOUND)
    add_definitions("-DENABLE_ALUGRID=1")
    include_directories(${ALUGRID_INCLUDES})
  endif (ALUGRID_FOUND)

  # libs
  set(COMMON_LIBS ${DUNE_DEFAULT_LIBS})

  add_subdirectory(test EXCLUDE_FROM_ALL)
  add_subdirectory(cmake/modules)

  #disable most warnings from dependent modules
  foreach(_mod ${ALL_DEPENDENCIES})
      dune_module_to_uppercase(_upper_case "${_mod}")
      if(${_mod}_INCLUDE_DIRS)
          foreach( _idir ${${_mod}_INCLUDE_DIRS} )
              add_definitions("-isystem ${_idir}")
          endforeach( _idir )
      endif(${_mod}_INCLUDE_DIRS)
  endforeach(_mod DEPENDENCIES)

  # finalize the dune project, e.g., generate config.h etc.
  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

endif(CMAKE_WITH_AUTOTOOLS)
