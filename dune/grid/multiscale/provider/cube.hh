#ifndef DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH

//#ifdef HAVE_DUNE_STUFF

// dune-common
#include <dune/common/parametertree.hh>

// dune-stuff
#include <dune/stuff/grid/provider/cube.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace Provider {

template< class GridImp >
class Cube
  : public Dune::Stuff::Grid::Provider::Cube< GridImp >
{
public:
  typedef Dune::Stuff::Grid::Provider::Cube< GridImp > BaseType;

  typedef typename BaseType::GridType GridType;

  typedef Cube< GridType > ThisType;

  Cube(Dune::ParameterTree paramTree)
    : BaseType(paramTree)
  {
  }

}; // class Cube

} // namespace Provider

} // namespace Multiscale

} // namespace Grid

} // namespace Dune

//#endif // HAVE_DUNE_STUFF

#endif // DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
