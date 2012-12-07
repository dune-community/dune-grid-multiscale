#ifndef DUNE_GRID_MULTISCALE_PROVIDER_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include "provider/interface.hh"
#include "provider/cube.hh"

namespace Dune {
namespace grid {
namespace Multiscale {
namespace Provider {

#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridType = Dune::GridSelector::GridType >
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridType = Dune::SGrid< 2, 2 > >
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
Interface< GridType >* create(const std::string& type = "grid.multiscale.provider.cube",
                              const Dune::ParameterTree paramTree = Dune::ParameterTree())
{
  // choose provider
  if (type == "grid.multiscale.provider.cube") {
    typedef Dune::grid::Multiscale::Provider::Cube< GridType > CubeProviderType;
    CubeProviderType* cubeProvider = new CubeProviderType(CubeProviderType::createFromParamTree(paramTree));
    return cubeProvider;
  } else
    DUNE_THROW(Dune::RangeError,
               "\nERROR: unknown grid provider '" << type << "' requested!");
} // Interface< GridImp >* create(const std::string& type, const Dune::ParameterTree paramTree = Dune::ParameterTree())

} // namespace Provider
} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_HH
