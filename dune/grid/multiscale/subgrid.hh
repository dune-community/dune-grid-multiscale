#ifndef DUNE_GRID_MULTISCALE_SUBGRID_HH
#define DUNE_GRID_MULTISCALE_SUBGRID_HH

// system
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-fem
#include <dune/fem/gridpart/gridpart.hh>

// dune-subgrid
#include <dune/subgrid/subgrid.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

template< class GridImp >
class Subgrid
{
public:
  typedef GridImp HostGridType;

  typedef Subgrid< HostGridType > ThisType;

  static const unsigned int dim = HostGridType::dimension;

  typedef Dune::SubGrid< dim, HostGridType > LocalGridType;

private:
  typedef std::map< unsigned int, Dune::shared_ptr< LocalGridType > > LocalGridMapType;

public:
  typedef Dune::LeafGridPart< HostGridType > GridPartType;

  Subgrid(HostGridType& hostGrid)
    : hostGrid_(hostGrid),
      finalized_(false)
  {}

  bool localGridExists(const unsigned int subdomain) const
  {
    std::cout << "subdomain: " << subdomain << std::endl;
    std::cout << "map: ";
    for (typename LocalGridMapType::const_iterator it = localGridMap_.begin(); it != localGridMap_.end(); ++it) {
        std::cout << it->first << " ";
    }
    std::cout << std::endl;
    return localGridMap_.find(subdomain) != localGridMap_.end();
  }

  void createLocalGrid(const unsigned int subdomain)
  {
    assert(!finalized_);
    if (localGridMap_.find(subdomain) != localGridMap_.end()) {
        localGridMap_[subdomain] = Dune::shared_ptr< LocalGridType >(new LocalGridType(hostGrid_));
        localGridMap_[subdomain]->createBegin();
    }
  } // void createLocalGrid(const unsigned int subdomain)

  LocalGridType& localGrid(const unsigned int subdomain)
  {
    assert(!finalized_);
    assert(localGridExists(subdomain));
    return *(localGridMap_[subdomain]);
  }

  const LocalGridType& localGrid(const unsigned int subdomain) const
  {
    assert(localGridExists(subdomain));
    return *(localGridMap_[subdomain]);
  }

private:
  HostGridType& hostGrid_;
  bool finalized_;
  LocalGridMapType localGridMap_;
}; // class Subgrid

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_SUBGRID_HH
