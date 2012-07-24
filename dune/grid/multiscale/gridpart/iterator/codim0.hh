
#ifndef DUNE_GRID_MULTISCALE_GRIDPART_ITERATOR_CODIM0_HH
#define DUNE_GRID_MULTISCALE_GRIDPART_ITERATOR_CODIM0_HH

// system
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-grid
#include <dune/grid/common/grid.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace GridPart {

namespace Iterator {

namespace Codim0 {

template< class GridPartImp, class GlobalIndexImp, int PartitionIteratorImp >
class IndexBased
{
public:
  IndexBased()
  {
    DUNE_THROW(Dune::NotImplemented, "Error: only implemented for PartitionIteratorType = InteriorBorder_Partition!");
  }
}; // class IndexBased

template< class GridPartImp, class GlobalIndexImp >
class IndexBased< GridPartImp, GlobalIndexImp, InteriorBorder_Partition >
//  : public GridPartImp::GridType::template Codim< 0 >::template Partition< InteriorBorder_Partition >::LeafIterator
{
public:
  typedef GridPartImp GridPartType;

  typedef typename GridPartType::GridType GridType;

  typedef GlobalIndexImp GlobalIndexType;

  typedef IndexBased< GridPartType, GlobalIndexType, InteriorBorder_Partition > ThisType;

  typedef typename GridType::template Codim< 0 >::template Partition< InteriorBorder_Partition >::LeafIterator HostIteratorType;

  typedef typename std::set< GlobalIndexType >::iterator IndexIteratorType;

  typedef typename GridType::template Codim< 0 >::Entity Entity;

  IndexBased(const GridPartType& gridPart, const bool end = false)
    : gridPart_(gridPart),
      latest_(new IndexIteratorType(gridPart_.globalIndicesSet_->begin())),
      last_(*(gridPart_.globalIndicesSet_->rbegin())),
      end_(gridPart_.globalIndicesSet_->end()),
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
      GlobalIndexType index = gridPart_.indexSet().index(hostIterator_->operator*());
      IndexIteratorType* result = new IndexIteratorType(std::find(*latest_, end_, index));
      if ((*result != end_)) {
        found = true;
        if (**result == last_)
          workAtAll_ = false;
        else
          latest_ = result;
      } else
        hostIterator_->operator++();
    }
  } // void forward()

  const GridPartType& gridPart_;
  IndexIteratorType* latest_;
  const GlobalIndexType last_;
  const IndexIteratorType end_;
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
