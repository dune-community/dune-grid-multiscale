#ifndef DUNE_GRID_MULTISCALE_FACTORY_SUBGRID_HH
#define DUNE_GRID_MULTISCALE_FACTORY_SUBGRID_HH

// system
#include <sstream>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

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

  static const std::string id;

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
    if (!msGrid_->localGridExists(subdomain))
      msGrid_->createLocalGrid(subdomain);

    // add entity to subdomain
    msGrid_->localGrid(subdomain).insertPartial(entity);
  } // void add(const EntityType& entity, const unsigned int subdomain)

  void finalize()
  {
    // test for consecutive numbering of subdomains and finalize subgrids
    bool consecutive = true;
    for (unsigned int subdomain = 0; subdomain < msGrid_->numSubdomains(); ++subdomain) {
      if (msGrid_->localGridExists(subdomain))
        msGrid_->localGrid(subdomain).createEnd();
      else
        consecutive = false;
    }
    if (consecutive)
      msGrid_->finalize();
    else {
      std::stringstream msg;
      msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
  } // void finalize()

private:
  HostGridType& hostGrid_;
  Dune::shared_ptr< MsGridType > msGrid_;
}; // class FromFineGrid

template< class GridType >
const std::string FromFineGrid< GridType >::id = "grid.multiscale.factory.subgrid.fromfinegrid";

template< class GridImp >
class FromCoarseGrid
{
public:
  typedef GridImp HostGridType;

  typedef FromCoarseGrid< HostGridType > ThisType;

  static const std::string id;

  typedef Dune::grid::Multiscale::Subgrid< HostGridType > MsGridType;

  typedef typename MsGridType::GlobalEntityType GlobalEntityType;

  FromCoarseGrid(HostGridType& hostGrid)
    : hostGrid_(hostGrid),
      msGrid_(Dune::shared_ptr< MsGridType >(new MsGridType(hostGrid_)))
  {
  }

  const HostGridType& hostGrid() const
  {
    return hostGrid_;
  }

  void prepare()
  {
    return;
  }

  void add(const GlobalEntityType& entity, const unsigned int subdomain, const unsigned int refineLevel = 0)
  {
    // create local grid if needed
    if (!msGrid_->localGridExists(subdomain))
      msGrid_->createLocalGrid(subdomain);

    // add entity to subdomain
    msGrid_->localGrid(subdomain).insert(entity);

    return;
  } // void add(const EntityType& entity, const unsigned int subdomain)

  Dune::shared_ptr< MsGridType > finalize()
  {
    // test for consecutive numbering of subdomains and finalize subgrids
    bool consecutive = true;
    for (unsigned int subdomain = 0; subdomain < msGrid_->numSubdomains(); ++subdomain) {
      if (msGrid_->localGridExists(subdomain))
        msGrid_->localGrid(subdomain).createEnd();
      else
        consecutive = false;
    }
    if (consecutive)
      msGrid_->finalize();
    else {
      std::stringstream msg;
      msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return msGrid_;
  } // void finalize()

private:
  HostGridType& hostGrid_;
  Dune::shared_ptr< MsGridType > msGrid_;
}; // class FromCoarseGrid

template< class GridType >
const std::string FromCoarseGrid< GridType >::id = "grid.multiscale.factory.subgrid.fromfinegrid";

} // namespace Subgrid

} // namespace Factory

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_FACTORY_SUBGRID_HH