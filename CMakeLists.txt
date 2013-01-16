PROJECT(dune-grid-multiscale)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )

SET(ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
SET(ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
SET(ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

include(DuneUtils)
include(UseDoxygen)

#needs to come first to include correct config.h
INCLUDE_DIR(${dune-grid-multiscale_BINARY_DIR})
INCLUDE_DIR(${dune-grid-multiscale_SOURCE_DIR})

# General Settings
SET( CMAKE_COLOR_MAKEFILE ON )

# Options, that can be changed by the User in order to customise the module
set(ENABLE_MPI
    "0" CACHE STRING
    "set to 1 to enable mpi")

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
  "ALUGRID BASE PATH" )

SET( ALBERTA_GRID_BASE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../local" CACHE STRING
  "ALBERTA GRID BASE_PATH" )

SET( UGGRID_BASE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../local" CACHE STRING
  "UGGRID BASE PATH" )

SET( CUSTOM_FLAGS
  "-Wall -O0 -DDEBUG -g -ggdb -fdiagnostics-show-option" CACHE LIST
  "CUSTOM FLAGS" )

SET( ENABLE_ALUGRID OFF )
IF( ${GRIDTYPE} MATCHES "ALUGRID" )
  INCLUDE_SYS_DIR(  ${GRAPE_PATH}
    ${ALUGRID_BASE_PATH}/include
    ${ALUGRID_BASE_PATH}/include/serial
    ${ALUGRID_BASE_PATH}/include/duneinterface )
  LINK_DIRECTORIES( "${ALUGRID_BASE_PATH}/lib" )
  SET( GRIDLIBS "alugrid" )
  SET( ENABLE_ALUGRID ON )
  IF ( ENABLE_PARALELL )
    INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
  ENDIF ( ENABLE_PARALELL )
  ADD_DEFINITIONS( -DENABLE_ALUGRID )
ENDIF( ${GRIDTYPE} MATCHES "ALUGRID" )

SET( ENABLE_ALBERTA OFF )
IF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )
  INCLUDE_SYS_DIR(  ${ALBERTA_GRID_BASE_PATH}/include
    ${ALBERTA_GRID_BASE_PATH}/include/alberta )
  LINK_DIRECTORIES( "${ALBERTA_GRID_BASE_PATH}/lib" )
  SET( GRIDLIBS alberta_utilities alberta_2d dunealbertagrid_2d dunegrid)
  SET( ENABLE_ALBERTA ON )
  ADD_DEFINITIONS( -DENABLE_ALBERTA=1 -DDEBUG )
ENDIF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )
ADD_DEFINITIONS(${CUSTOM_FLAGS})

# the other dune modules
ADD_DUNE_MODULES(common grid geometry stuff)
set( COMMON_LIBS "dunegrid" "dunecommon" "dunegeometry"
    ${BLAS_LIB} ${GRIDLIBS} ${ParaLIBS}
        "boost_filesystem" "boost_system" )
FILE( GLOB_RECURSE header "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )
set( COMMON_HEADER ${header} ${DUNE_HEADERS} )

set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

add_subdirectory(doc)
add_subdirectory(examples)

#HEADERCHECK( ${header} )
#DEPENDENCYCHECK( ${header} )
