# This file is part of the dune-grid-multiscale project:
#   http://users.dune-project.org/projects/dune-grid-multiscale
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

ADD_EXECUTABLE( examples_provider "provider.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_provider ${COMMON_LIBS} )
