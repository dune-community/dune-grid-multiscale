# This file is part of the dune-grid-multiscale project:
#   http://users.dune-project.org/projects/dune-grid-multiscale
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

set(lib_dune_grid_multiscale_sources
    exists.cc
    default_2d_aluconform.cc)

dune_add_library("dunegridmultiscale" ${lib_dune_grid_multiscale_sources}
                 ADD_LIBS ${DUNE_LIBS})

target_link_dune_default_libraries(dunegridmultiscale)

if(dune-grid_FOUND)
  add_dune_alugrid_flags(dunegridmultiscale)
endif(dune-grid_FOUND)

