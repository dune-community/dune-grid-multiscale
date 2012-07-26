
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

template< class LocalGridPartImp >
class IndexBased
{
public:
  typedef LocalGridPartImp LocalGridPartType;

  typedef IndexBased< LocalGridPartType > ThisType;

private:
  typedef typename LocalGridPartType::GlobalGridPartType::template Codim< 0 >::IteratorType HostIteratorType;

  typedef Dune::GeometryType GeometryType;

  typedef typename LocalGridPartType::IndexSetType::IndexType IndexType;

  typedef std::map< IndexType, IndexType > IndexMapType;

  typedef typename IndexMapType::iterator IndexPairType;

  typedef std::map< GeometryType, IndexMapType > GeometryMapType;
public:
  typedef typename LocalGridPartType::template Codim< 0 >::EntityType EntityType;

  IndexBased(const LocalGridPartType& localGridPart, const bool end = false)
    : localGridPart_(localGridPart),
      workAtAll_(0)
  {
    if (end)
      hostIterator_ = new HostIteratorType(localGridPart_.globalGridPart_.template end< 0 >());
    else {
      for (typename GeometryMapType::const_iterator iterator = localGridPart_.geometryMap_->begin();
           iterator != localGridPart_.geometryMap_->end();
           ++iterator) {
        if (iterator->first.dim() == LocalGridPartType::GridType::dimension) {
          ++workAtAll_;
          last_.insert(std::pair< GeometryType, IndexType >(iterator->first, iterator->second.rbegin()->first));
          end_.insert(std::pair< GeometryType, typename IndexMapType::const_iterator >(iterator->first, iterator->second.end()));
        }
      }
      hostIterator_ = new HostIteratorType(localGridPart_.globalGridPart_.template begin< 0 >());
      forward();
    }
  }

  ThisType& operator++()
  {
    if (workAtAll_ > 0) {
      hostIterator_->operator++();
      forward();
    } else
      hostIterator_ = new HostIteratorType(localGridPart_.globalGridPart_.template end< 0 >());
    return *this;
  }

  bool operator==(const ThisType& other) const
  {
    return *hostIterator_ == *(other.hostIterator_);
  }

  bool operator!=(const ThisType& other) const
  {
    return *hostIterator_ != *(other.hostIterator_);
  }

  EntityType& operator*() const
  {
    return hostIterator_->operator*();
  }

private:
  void forward()
  {
    // iterate forward until we find the next entity that belongs to the local grid part
    bool found = false;
    while (!found && (workAtAll_ > 0)) {
      const EntityType& entity = hostIterator_->operator*();
      const IndexType index = localGridPart_.globalGridPart_.indexSet().index(entity);
      const GeometryType geometryType = entity.type();
      typename GeometryMapType::const_iterator indexMap = localGridPart_.geometryMap_->find(geometryType);
      if (indexMap != localGridPart_.geometryMap_->end()) {
        const typename IndexMapType::const_iterator result = indexMap->second.find(index);
        if ((result != end_.find(geometryType)->second)) {
          found = true;
          if (result->first == last_.find(geometryType)->second)
            --workAtAll_;
        } else
          hostIterator_->operator++();
      } else
        hostIterator_->operator++();
    }
  } // void forward()

  const LocalGridPartType& localGridPart_;
  unsigned int workAtAll_;
  std::map< GeometryType, IndexType > last_;
  std::map< GeometryType, typename IndexMapType::const_iterator > end_;
  HostIteratorType* hostIterator_;
}; // class IndexBased

} // namespace Codim0

} // namespace Local

} // namespace Iterator

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_ITERATOR_CODIM0_HH
