
#ifndef DUNE_GRID_MULTISCALE_GRIDPART_INDEXSET_LOCAL_HH
#define DUNE_GRID_MULTISCALE_GRIDPART_INDEXSET_LOCAL_HH

// system
#include <vector>

// dune-common
#include <dune/common/exceptions.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/gridpart/indexset/default.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace GridPart {

namespace IndexSet {

namespace Local {

template< class GridPartImp >
class IndexBased
{
public:
  typedef GridPartImp GridPartType;

  typedef IndexBased< GridPartType > ThisType;

  typedef typename GridPartType::GridType GridType;

  typedef typename GridPartType::EntityType EntityType;

private:
  typedef Dune::grid::Multiscale::GridPart::IndexSet::Default::Leaf< GridType > HostIndexSetType;

public:
  static const unsigned int dimension = HostIndexSetType::dimension;

  typedef typename HostIndexSetType::IndexType IndexType;

  IndexBased(const GridPartType& gridPart)
    : gridPart_(gridPart),
      hostIndexSetType_(gridPart_.globalGridPart_.indexSet())
  {}

  IndexType index(const EntityType& entity) const
  {
    assert(contains(entity));
    const IndexType globalIndex = hostIndexSetType_.index(entity);
    return gridPart_.globalToLocaIndexMap_->find(globalIndex)->second;
  }

  IndexType subIndex(const EntityType& /*entity*/, int /*i*/, unsigned int /*codim*/) const
  {
    DUNE_THROW(Dune::NotImplemented, "Will be implemented, as soon as I know what it does.");
    return -1;
  }

  const std::vector< GeometryType >& geomTypes(int /*codim*/) const
  {
    DUNE_THROW(Dune::NotImplemented, "Will be implemented, as soon as I know what it does.");
    return -1;
  }

  IndexType size(GeometryType /*type*/) const
  {
    DUNE_THROW(Dune::NotImplemented, "Will be implemented, as soon as I know what it does.");
    return -1;
  }

  IndexType size(int /*codim*/) const
  {
    DUNE_THROW(Dune::NotImplemented, "Will be implemented, as soon as I know what it does.");
    return -1;
  }

  bool contains(const EntityType& entity) const
  {
    const IndexType globalIndex = hostIndexSetType_.index(entity);
    return gridPart_.globalToLocaIndexMap_->find(globalIndex) != gridPart_.globalToLocaIndexMap_->end();
  }

private:
  const GridPartType& gridPart_;
  const HostIndexSetType& hostIndexSetType_;
}; // class IndexBased

} // namespace IndexSet

} // namespace IndexSet

} // namespace GridPart

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_GRIDPART_INDEXSET_LOCAL_HH
