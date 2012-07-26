
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
  {
    // the most stupid way of checking that we have only one geometry
    assert(hostIndexSetType_.geomTypes(0).size() == 1);
  }

  IndexType index(const EntityType& entity) const
  {
//    assert(contains(entity));
    const IndexType globalIndex = hostIndexSetType_.index(entity);
    return globalIndex;
//    return gridPart_.globalToLocaIndexMap_->find(globalIndex)->second;
  }

  IndexType subIndex(const EntityType& entity, int i, unsigned int codim) const
  {
    const IndexType ret = hostIndexSetType_.subIndex(entity, i, codim);
    std::cout << "    IndexBased::subIndex(entity, " << i << ", " << codim << "): " << ret << std::endl;
    return ret;
  }

  //! \attention Not thought about this yet!
  //! \todo Think about this!
  const std::vector< GeometryType >& geomTypes(int codim) const
  {
    const std::vector< GeometryType >& ret = hostIndexSetType_.geomTypes(codim);
    std::cout << "  IndexBased::geomTypes(" << codim << "): ";
    for (unsigned int i = 0; i < ret.size(); ++i)
      std::cout << ret[i] << " ";
    std::cout << std::endl;
    return ret;
  }

  //! \attention Not thought about this yet!
  //! \todo Think about this!
  IndexType size(GeometryType type) const
  {
    const IndexType ret = hostIndexSetType_.size(type);
    std::cout << "    IndexBased::size(" << type << "): " << ret << std::endl;
//    return gridPart_.globalToLocaIndexMap_->size();
    return ret;
  }

  //! \attention Not thought about this yet!
  //! \todo Think about this!
  IndexType size(int codim) const
  {
    const IndexType ret = hostIndexSetType_.size(codim);
    std::cout << "  IndexBased::size(" << codim << "): " << ret << std::endl;
//    assert(codim == 0);
//    return gridPart_.globalToLocaIndexMap_->size();
    return ret;
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
