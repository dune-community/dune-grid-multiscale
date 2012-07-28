
#ifndef DUNE_GRID_PART_ITERATOR_CODIM0_HH
#define DUNE_GRID_PART_ITERATOR_CODIM0_HH

// system
#include <map>

// dune-common
#include <dune/common/exceptions.hh>

// dune-geometry
#include <dune/geometry/type.hh>

// dune-grid
#include <dune/grid/common/grid.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Iterator {

namespace Local {

namespace Codim0 {

/**
 *  \brief  Iterates over those entities of a grid part, the indices of which match predefined ones.
 *  \todo   Replace GridPartImp with Interface< GridPartTraitsImp >!
 *  \todo   Document!
 */
template< class GridPartImp >
class IndexBased
{
public:
  typedef GridPartImp GridPartType;

  typedef IndexBased< GridPartType > ThisType;

  typedef typename GridPartType::IndexSetType::IndexType IndexType;

private:
  typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;

  typedef Dune::GeometryType GeometryType;

  typedef std::map< IndexType, IndexType > IndexMapType;

public:
  typedef std::map< GeometryType, std::map< IndexType, IndexType > > IndexContainerType;

  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

  IndexBased(const GridPartType& globalGridPart,
             const Dune::shared_ptr< const IndexContainerType > indexContainer,
             const bool end = false)
    : globalGridPart_(globalGridPart),
      indexContainer_(indexContainer),
      workAtAll_(0)
  {
    if (end)
      iterator_ = new IteratorType(globalGridPart_.template end< 0 >());
    else {
      for (typename IndexContainerType::const_iterator iterator = indexContainer_->begin();
           iterator != indexContainer_->end();
           ++iterator) {
        if (iterator->first.dim() == GridPartType::GridType::dimension) {
          ++workAtAll_;
          last_.insert(std::pair< GeometryType, IndexType >(iterator->first, iterator->second.rbegin()->first));
          end_.insert(std::pair< GeometryType, typename IndexMapType::const_iterator >(iterator->first, iterator->second.end()));
        }
      }
      iterator_ = new IteratorType(globalGridPart_.template begin< 0 >());
      forward();
    }
  }

  ThisType& operator++()
  {
    if (workAtAll_ > 0) {
      iterator_->operator++();
      forward();
    } else
      iterator_ = new IteratorType(globalGridPart_.template end< 0 >());
    return *this;
  }

  bool operator==(const ThisType& other) const
  {
    return *iterator_ == *(other.iterator_);
  }

  bool operator!=(const ThisType& other) const
  {
    return *iterator_ != *(other.iterator_);
  }

  EntityType& operator*() const
  {
    return iterator_->operator*();
  }

private:
  //! iterates forward until we find the next entity that belongs to the local grid part
  void forward()
  {
    bool found = false;
    while (!found && (workAtAll_ > 0)) {
      const EntityType& entity = iterator_->operator*();
      const IndexType index = globalGridPart_.indexSet().index(entity);
      const GeometryType geometryType = entity.type();
      typename IndexContainerType::const_iterator indexMap = indexContainer_->find(geometryType);
      if (indexMap != indexContainer_->end()) {
        const typename IndexMapType::const_iterator result = indexMap->second.find(index);
        if ((result != end_.find(geometryType)->second)) {
          found = true;
          if (result->first == last_.find(geometryType)->second)
            --workAtAll_;
        } else
          iterator_->operator++();
      } else
        iterator_->operator++();
    }
  } // void forward()

  const GridPartType& globalGridPart_;
  const Dune::shared_ptr< const IndexContainerType > indexContainer_;
  unsigned int workAtAll_;
  std::map< GeometryType, IndexType > last_;
  std::map< GeometryType, typename IndexMapType::const_iterator > end_;
  IteratorType* iterator_;
}; // class IndexBased

} // namespace Codim0

} // namespace Local

} // namespace Iterator

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_ITERATOR_CODIM0_HH
