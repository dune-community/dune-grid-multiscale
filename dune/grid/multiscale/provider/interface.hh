#ifndef DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#else
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/shared_ptr.hh>

#include <dune/grid/multiscale/default.hh>

namespace Dune {
namespace grid {
namespace Multiscale {
namespace Provider {

#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridImp = Dune::GridSelector::GridType >
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridImp = Dune::SGrid< 2, 2 > >
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class Interface
{
public:
  typedef GridImp GridType;

  typedef Interface< GridType > ThisType;

  typedef Dune::grid::Multiscale::Default< GridType > MsGridType;

  static const std::string id()
  {
    return "grid.multiscale.provider";
  }

  virtual const Dune::shared_ptr< const GridType > grid() const = 0;

  virtual const Dune::shared_ptr< const MsGridType > msGrid() const = 0;

//  virtual void visualize(const std::string filename = id()) const;
}; // class Interface

} // namespace Provider
} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
