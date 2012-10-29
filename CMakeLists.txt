PROJECT(dune-grid-multiscale)

SET( CMAKE_CXX_COMPILER
  "g++-4.6" CACHE FORCE
    "" )

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0009 NEW)
# tell cmake where to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

# stuff needed for dune
include(DuneUtils)
include(UseDoxygen)
include (FindPkgConfig)

# General Settings
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a" )
SET( CMAKE_COLOR_MAKEFILE ON )

# Options, that can be changed by the User in order to customise the module
SET( BLAS_LIB
  "blas" CACHE STRING
  "BLAS LIB" )

SET( GRIDTYPE
  "SGRID" CACHE STRING
  "GRIDTYPE" )

SET( GRIDDIM
  "2" CACHE STRING
  "GRIDDIM" )

SET( ALUGRID_BASE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../local" CACHE STRING
  "ALUGRID_BASE_PATH" )

SET( CUSTOM_FLAGS
  "-Wall -O0 -DDEBUG -g -fdiagnostics-show-option" CACHE LIST
  "CUSTOM FLAGS" )

set( BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../local" )
find_package( Boost 1.48.0 REQUIRED )
include_directories(${Boost_INCLUDE_DIR})
link_directories( ${Boost_LIBRARY_DIRS} )

# needs to come first to include correct config.h
include_dir( ${dune-grid-multiscale_SOURCE_DIR} )
add_cxx_flags( -DHAVE_CONFIG_H )
add_cxx_flags( -D${GRIDTYPE} )
add_cxx_flags( -DGRIDDIM=${GRIDDIM} )
add_cxx_flags( -I${ALUGRID_BASE_PATH}/include
  -I${ALUGRID_BASE_PATH}/include/serial
  -I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )

add_cxx_flags(-DBOOST_DATE_TIME_NO_LOCALE)
link_directories( "${ALUGRID_BASE_PATH}/lib" )

add_cxx_flags(-std=c++0x)
add_cxx_flags(${CUSTOM_FLAGS})

# the other dune modules
ADD_DUNE_MODULES(common grid geometry stuff fem)

# Source listing, relative to the position of this file
FILE( GLOB_RECURSE common "../dune-common/*.hh" )
FILE( GLOB_RECURSE common_cc "../dune-common/*.cc" )
FILE( GLOB_RECURSE grid "../dune-grid/*.hh" )
FILE( GLOB_RECURSE grid_cc "../dune-grid/*.cc" )
FILE( GLOB_RECURSE stuff "../dune-stuff/*.hh" )
FILE( GLOB_RECURSE stuff_cc "../dune-stuff/*.cc" )
FILE( GLOB_RECURSE grid-multiscale "../dune-grid-multiscale/*.hh" )

set_source_files_properties( ${common} ${grid} ${stuff} ${subgrid} ${grid-multiscale}
${common_cc} ${grid_cc} ${stuff_cc} ${subgrid_cc} PROPERTIES HEADER_FILE_ONLY 1 )

set( LIST_LIBS ${common} ${grid} ${stuff} ${subgrid} ${grid-multiscale})

set( COMMON_LIBS "grid" "dunegeometry" "common" ${BLAS_LIB} "alugrid" "boost_filesystem" "boost_system" )

ADD_EXECUTABLE( examples_provider_cube "examples/provider/cube.cc" ${LIST_LIBS} )
TARGET_LINK_LIBRARIES( examples_provider_cube ${COMMON_LIBS} )