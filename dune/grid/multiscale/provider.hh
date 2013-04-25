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
        "gridprovider.multiscale.cube"
    };
  } // ... available(...)

  static Dune::ParameterTree createSampleDescription(const std::string type, const std::string subname = "")
  {
  if (type == "gridprovider.multiscale.cube")
    return ProviderCube< GridType >::createSampleDescription(subname);
  else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
               << " unknown multiscale gridprovider '" << type << "' requested!");
  } // ... createSampleDescription(...)

  static ProviderInterface< GridType >* create(const std::string& type = available()[0],
                                               const Dune::ParameterTree description = Dune::ParameterTree())
  {
    // choose provider
    if (type == "gridprovider.multiscale.cube")
      return Dune::grid::Multiscale::ProviderCube< GridType >::create(description);
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
