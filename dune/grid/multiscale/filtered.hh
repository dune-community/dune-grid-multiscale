#ifndef DUNE_GRID_MULTISCALE_FILTERED_HH
#define DUNE_GRID_MULTISCALE_FILTERED_HH

// system
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-fem
#include <dune/fem/gridpart/filteredgrid.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

template< class GridImp >
class Filtered
{
public:
  typedef GridImp HostGridType;

  typedef Filtered< HostGridType > ThisType;

  static const unsigned int dim = HostGridType::dimension;

//  typedef Dune::SubGrid< dim, HostGridType > LocalGridType;

//private:
//  typedef std::map< unsigned int, Dune::shared_ptr< LocalGridType > > LocalGridMapType;

//public:
//  typedef Dune::LeafGridPart< HostGridType > GlobalGridPartType;

//  typedef typename GlobalGridPartType::template Codim< 0 >::IteratorType::Entity GlobalEntityType;

//  typedef Dune::LeafGridPart< LocalGridType > LocalGridPartType;

//  typedef typename LocalGridPartType::template Codim< 0 >::IteratorType::Entity LocalEntityType;

  Filtered(HostGridType& hostGrid)
    : hostGrid_(hostGrid)/*,
      finalized_(false),
      numSubdomains_(0)*/
  {}

//  bool localGridExists(const unsigned int subdomain) const
//  {
//    return localGridMap_.find(subdomain) != localGridMap_.end();
//  }

//  void createLocalGrid(const unsigned int subdomain)
//  {
//    assert(!finalized_);
//    if (localGridMap_.find(subdomain) == localGridMap_.end()) {
//        localGridMap_[subdomain] = Dune::shared_ptr< LocalGridType >(new LocalGridType(hostGrid_));
//        localGridMap_[subdomain]->createBegin();
//        ++numSubdomains_;
//    }
//  } // void createLocalGrid(const unsigned int subdomain)

//  unsigned int numSubdomains() const
//  {
//    return numSubdomains_;
//  }

//  void finalize()
//  {
//    finalized_ = true;
//  }

//  LocalGridType& localGrid(const unsigned int subdomain)
//  {
//    assert(localGridExists(subdomain));
//    return *(localGridMap_[subdomain]);
//  }

//  const LocalGridType& localGrid(const unsigned int subdomain) const
//  {
//    assert(localGridExists(subdomain));
//    return *(localGridMap_[subdomain]);
//  }

//  GlobalGridPartType globalGridPart()
//  {
//    return GlobalGridPartType(hostGrid_);
//  }

//  LocalGridPartType localGridPart(const unsigned int subdomain)
//  {
//    assert(localGridExists(subdomain));
//    return LocalGridPartType(*(localGridMap_[subdomain]));
//  }

//  const GlobalEntityType& globalEntity(const unsigned int subdomain, const LocalEntityType& localEntity) const
//  {
//    assert(localGridExists(subdomain));
//    return *(localGrid(subdomain).template getHostEntity< 0 >(localEntity));
//  }

private:
  HostGridType& hostGrid_;
//  bool finalized_;
//  unsigned int numSubdomains_;
//  LocalGridMapType localGridMap_;
}; // class Filtered

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_FILTERED_HH
