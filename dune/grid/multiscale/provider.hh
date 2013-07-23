// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/exceptions.hh>

#include <dune/grid/sgrid.hh>

#include <dune/stuff/common/color.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"
#include "provider/functionbased.hh"

namespace Dune {
namespace grid {
namespace Multiscale {

#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridType = Dune::GridSelector::GridType >
#else
template< class GridType = Dune::SGrid< 2, 2 > >
#endif
class Providers
{
public:
  static std::vector< std::string > available()
  {
    return {
        "grid.multiscale.provider.cube",
        "grid.multiscale.provider.functionbased"
    };
  } // ... available(...)

  static Dune::ParameterTree createSampleDescription(const std::string type, const std::string subname = "")
  {
  if (type == "grid.multiscale.provider.cube")
    return ProviderCube< GridType >::createSampleDescription(subname);
  else if (type == "grid.multiscale.provider.functionbased")
    return ProviderFunctionbased< GridType >::createSampleDescription(subname);
  else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << " unknown multiscale gridprovider '" << type << "' requested!");
  } // ... createSampleDescription(...)

  static ProviderInterface< GridType >* create(const std::string& type = available()[0],
                                               const Dune::ParameterTree description = Dune::ParameterTree())
  {
    // choose provider
    if (type == "grid.multiscale.provider.cube")
      return Dune::grid::Multiscale::ProviderCube< GridType >::create(description);
    else if (type == "grid.multiscale.provider.functionbased")
      return Dune::grid::Multiscale::ProviderFunctionbased< GridType >::create(description);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " unknown multiscale gridprovider '" << type << "' requested!");
  } // ... create(...)
}; // class Providers


} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_HH
