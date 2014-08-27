// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include "default.hh"

namespace Dune {
namespace grid {
namespace Multiscale {

//#if HAVE_ALUGRID


template class Default< ALUGrid< 2, 2, simplex, conforming > >;


//#endif // HAVE_ALUGRID

} // namespace Multiscale
} // namespace grid
} // namespace Dune

