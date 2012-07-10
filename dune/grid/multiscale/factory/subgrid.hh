#ifndef DUNE_GRID_MULTISCALE_FACTORY_SUBGRID_HH
#define DUNE_GRID_MULTISCALE_FACTORY_SUBGRID_HH

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/subgrid.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace Factory {

namespace Subgrid {

template< class GridImp >
class FromFineGrid
{
public:
  typedef GridImp HostGridType;

  typedef FromFineGrid< HostGridType > ThisType;

  typedef Dune::grid::Multiscale::Subgrid< HostGridType > MsGridType;

  typedef typename MsGridType::GridPartType GridPartType;

  typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity EntityType;

  FromFineGrid(HostGridType& hostGrid)
    : hostGrid_(hostGrid),
      msGrid_(Dune::shared_ptr< MsGridType >(new MsGridType(hostGrid_)))
  {
  }

  const HostGridType& hostGrid() const
  {
    return hostGrid_;
  }

  void prepare()
  {}

  void add(const EntityType& entity, const unsigned int subdomain)
  {
    // create local grid if needed
    if (!msGrid_.localGridExists(subdomain))
      msGrid_.createLocalGrid(subdomain);

    // add entity to subdomain
    msGrid_.localGrid(subdomain).insertPartial(entity);
  } // void add(const EntityType& entity, const unsigned int subdomain)


private:
  HostGridType& hostGrid_;
  Dune::shared_ptr< MsGridType > msGrid_;
}; // class FromFineGrid

} // namespace Subgrid

} // namespace Factory

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_FACTORY_SUBGRID_HH
