#ifndef DUNE_GRID_MULTISCALE_FILTERED_HH
#define DUNE_GRID_MULTISCALE_FILTERED_HH

// system
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-fem
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/filteredgrid.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

template< class GridImp, class FilterImp >
class Filtered
{
public:
  typedef GridImp HostGridType;

  typedef FilterImp FilterType;

  typedef Filtered< HostGridType, FilterType > ThisType;

  static const std::string id;

  static const unsigned int dim = HostGridType::dimension;

//  typedef Dune::LeafGridPart< HostGridType > GlobalGridPartType;

  typedef Dune::FilteredGridPart< Dune::LeafGridPart< HostGridType >, FilterType > LocalGridPartType;

  Filtered(HostGridType& hostGrid)
    : hostGrid_(hostGrid),
      finalized_(false),
      size_(0)
  {}

  const unsigned int size() const
  {
    return size_;
  }

  bool exists(const unsigned int subdomain) const
  {
    return localGridPartMap_.find(subdomain) != localGridPartMap_.end();
  }

  const Dune::shared_ptr< const LocalGridPartType > localGridPart(const unsigned int subdomain) const
  {
    assert(finalized_);
    if (!exists(subdomain)) {
      std::stringstream msg;
      msg << "Error in " << id << ": map broken!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return localGridPartMap_.find(subdomain)->second;
  }

  void create(Dune::shared_ptr< FilterType > filter, const unsigned int subdomain)
  {
    assert(!finalized_);
    assert(!exists(subdomain));
    filterMap_[subdomain] = filter;
    localGridPartMap_[subdomain] = Dune::shared_ptr< LocalGridPartType >(new LocalGridPartType(hostGrid_, *(filterMap_[subdomain])));
    ++size_;
  }

  void finalize()
  {
    // test for consecutive numbering of subdomains and finalize subgrids
    bool consecutive = true;
    for (unsigned int subdomain = 0; subdomain < size_; ++subdomain) {
      if (!exists(subdomain))
        consecutive = false;
    }
    if (consecutive)
      finalized_ = true;
    else {
      std::stringstream msg;
      msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
  } // void finalize()

private:
  HostGridType& hostGrid_;
  bool finalized_;
  unsigned int size_;
  std::map< unsigned int, Dune::shared_ptr< FilterType > > filterMap_;
  std::map< unsigned int, Dune::shared_ptr< LocalGridPartType > > localGridPartMap_;
}; // class Filtered

template< class GridType, class FilterType >
const std::string Filtered< GridType, FilterType >::id = "grid.multiscale.filtered";

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_FILTERED_HH
