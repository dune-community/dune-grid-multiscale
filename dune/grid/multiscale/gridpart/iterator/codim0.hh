
#ifndef DUNE_GRID_MULTISCALE_GRIDPART_ITERATOR_CODIM0_HH
#define DUNE_GRID_MULTISCALE_GRIDPART_ITERATOR_CODIM0_HH

// system
#include <map>

// dune-common
#include <dune/common/exceptions.hh>

// dune-grid
#include <dune/grid/common/grid.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace GridPart {

namespace Iterator {

namespace Codim0 {

template< class GridPartImp, int PartitionIteratorImp >
class IndexBased
{
public:
  IndexBased()
  {
    DUNE_THROW(Dune::NotImplemented, "Error: only implemented for PartitionIteratorType = InteriorBorder_Partition!");
  }
}; // class IndexBased

template< class GridPartImp >
class IndexBased< GridPartImp,  InteriorBorder_Partition >
//  : public GridPartImp::GridType::template Codim< 0 >::template Partition< InteriorBorder_Partition >::LeafIterator
{
public:
  typedef GridPartImp GridPartType;

  typedef typename GridPartType::GridType GridType;

  typedef IndexBased< GridPartType, InteriorBorder_Partition > ThisType;

private:
  typedef typename GridType::template Codim< 0 >::template Partition< InteriorBorder_Partition >::LeafIterator HostIteratorType;

  typedef typename GridPartType::IndexSetType::IndexType IndexType;

  typedef typename std::map< IndexType, IndexType >::iterator IndexPairType;

public:
  typedef typename GridType::template Codim< 0 >::Entity Entity;

  IndexBased(const GridPartType& gridPart, const bool end = false)
    : gridPart_(gridPart),
//      latest_(new IndexPairType(gridPart_.globalToLocaIndexMap_->begin())),
      last_(gridPart_.globalToLocaIndexMap_->rbegin()->first),
      end_(gridPart_.globalToLocaIndexMap_->end()),
      workAtAll_(true)
  {
    if (end)
      hostIterator_ = new HostIteratorType(gridPart_.globalGridPart_.template end< 0 >());
    else {
      hostIterator_ = new HostIteratorType(gridPart_.globalGridPart_.template begin< 0 >());
      forward();
    }
  }

  ThisType& operator++()
  {
    if (workAtAll_) {
      hostIterator_->operator++();
      forward();
    } else
      hostIterator_ = new HostIteratorType(gridPart_.globalGridPart_.template end< 0 >());
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

  Entity& operator*() const
  {
    return hostIterator_->operator*();
  }

private:
  void forward()
  {
    // iterate forward until we find the next entity that belongs to the local grid part
    bool found = false;
    while (!found && workAtAll_) {
      IndexType index = gridPart_.globalGridPart_.indexSet().index(hostIterator_->operator*());
      IndexPairType result = gridPart_.globalToLocaIndexMap_->find(index);
      if ((result != end_)) {
        found = true;
        if (result->first == last_)
          workAtAll_ = false;
//          latest_ = result;
      } else
        hostIterator_->operator++();
    }
  } // void forward()

  const GridPartType& gridPart_;
//  IndexPairType* latest_;
  const IndexType last_;
  const IndexPairType end_;
  bool workAtAll_;
  HostIteratorType* hostIterator_;
}; // class IndexBased

} // namespace Codim0

} // namespace Iterator

} // namespace GridPart

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_GRIDPART_ITERATOR_CODIM0_HH
