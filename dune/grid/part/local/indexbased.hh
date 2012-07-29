
#ifndef DUNE_GRID_PART_LOCAL_INDEXBASED_HH
#define DUNE_GRID_PART_LOCAL_INDEXBASED_HH

// system
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-geometry
#include <dune/geometry/type.hh>

// dune-grid-multiscale
#include <dune/grid/part/interface.hh>
#include <dune/grid/part/iterator/local/indexbased.hh>
#include <dune/grid/part/indexset/local.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Local {

namespace IndexBased {

template< class GlobalGridPartImp >
class Const;

template< class GlobalGridPartImp >
struct ConstTraits
{
  typedef Dune::grid::Part::Interface< typename GlobalGridPartImp::Traits > GlobalGridPartType;

  typedef Dune::grid::Part::Local::IndexBased::Const< GlobalGridPartImp > GridPartType;

  typedef typename GlobalGridPartType::GridType GridType;

  typedef typename Dune::grid::Part::IndexSet::Local::IndexBased< GridType, typename GlobalGridPartType::IndexSetType > IndexSetType;

  template< int codim >
  struct Codim
  {
    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Dune::grid::Part::Iterator::Local::IndexBased< GlobalGridPartType, codim, pitype > IteratorType;
    };
  };

  static const PartitionIteratorType indexSetPartitionType = GlobalGridPartType::indexSetPartitionType;

  static const bool conforming = GlobalGridPartType::conforming;

  typedef typename GlobalGridPartType::IntersectionIteratorType IntersectionIteratorType;

}; // class ConstTraits

template< class GlobalGridPartImp >
class Const
  : public Dune::grid::Part::Interface< Dune::grid::Part::Local::IndexBased::ConstTraits< GlobalGridPartImp > >
{
public:
  typedef Const< GlobalGridPartImp > ThisType;

  typedef Dune::grid::Part::Local::IndexBased::ConstTraits< GlobalGridPartImp > Traits;

  typedef Dune::grid::Part::Interface< Traits > BaseType;

  typedef typename Traits::GridType GridType;

  typedef typename Traits::GlobalGridPartType GlobalGridPartType;

  typedef typename Traits::IndexSetType IndexSetType;

  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef typename IndexSetType::IndexType IndexType;

  typedef std::map< IndexType, IndexType > IndexMapType;

  typedef Dune::GeometryType GeometryType;

  //! container type for the indices
  typedef std::map< GeometryType, std::map< IndexType, IndexType > > IndexContainerType;

  Const(const GlobalGridPartType& globalGridPart,
        const Dune::shared_ptr< const IndexContainerType > indexContainer)
    : globalGridPart_(globalGridPart),
      indexContainer_(indexContainer),
      indexSet_(globalGridPart_.indexSet(), indexContainer_)
  {}

  const IndexSetType& indexSet() const
  {
    return indexSet_;
  }

  const GridType& grid() const
  {
    return globalGridPart_.grid();
  }

  const GlobalGridPartType& globalGridPart() const
  {
    return globalGridPart_;
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType begin() const
  {
    return typename BaseType::template Codim< codim >::IteratorType(globalGridPart_, indexContainer_);
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType begin() const
  {
    return typename BaseType::template Codim< codim >::template Partition< pitype >::IteratorType(globalGridPart_, indexContainer_);
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType end() const
  {
    return typename BaseType::template Codim< codim >::IteratorType(globalGridPart_, indexContainer_, true);
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType end() const
  {
    return typename BaseType::template Codim< codim >::template Partition< pitype >::IteratorType(globalGridPart_, indexContainer_, true);
  }

  //! \todo Think about this, is obviously needed!
  IntersectionIteratorType ibegin(const EntityType& en) const
  {
//    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return en.ileafbegin();
  }

  //! \todo Think about this, is obviously needed!
  IntersectionIteratorType iend(const EntityType& en) const
  {
//    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return en.ileafend();
  }

  int boundaryId(const IntersectionType& intersection) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return intersection.boundaryId();
  }

  int level() const
  {
    return globalGridPart_.level();
  }

  template< class DataHandleImp ,class DataType >
  void communicate(CommDataHandleIF< DataHandleImp, DataType > & data, InterfaceType iftype, CommunicationDirection dir) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    globalGridPart_.communicate(data,iftype,dir);
  }

private:
  const GlobalGridPartType& globalGridPart_;
  const Dune::shared_ptr< const IndexContainerType > indexContainer_;
  const IndexSetType indexSet_;
}; // class Const

} // namespace IndexBased

} // namespace Local

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_LOCAL_INDEXBASED_HH
