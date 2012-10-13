
#ifndef DUNE_GRID_PART_LOCAL_INDEXBASED_HH
#define DUNE_GRID_PART_LOCAL_INDEXBASED_HH

// system
#include <map>
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-geometry
#include <dune/geometry/type.hh>

// dune-grid-multiscale
#include <dune/grid/part/interface.hh>
#include <dune/grid/part/iterator/local/indexbased.hh>
#include <dune/grid/part/iterator/intersection/local.hh>
#include <dune/grid/part/iterator/intersection/wrapper.hh>
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

  typedef typename Dune::grid::Part::IndexSet::Local::IndexBased< GlobalGridPartType > IndexSetType;

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

  typedef Dune::grid::Part::Iterator::Intersection::Wrapper::FakeDomainBoundary< GlobalGridPartType > IntersectionIteratorType;
}; // class ConstTraits

/**
 *  \todo Wrap entity, so that entity.ileaf{begin,end}() returns i{begin,end}(entity)
 *  \todo Implement boundaryId(intersection) by adding a std::map< intersectionIndex, boundaryId >!
 */
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

  //! container type for the boundary information
  typedef std::map< IndexType, std::map< int, int > > BoundaryInfoContainerType;

  Const(const Dune::shared_ptr< const GlobalGridPartType > globalGridPart,
        const Dune::shared_ptr< const IndexContainerType > indexContainer,
        const Dune::shared_ptr< const BoundaryInfoContainerType > boundaryInfoContainer)
    : globalGridPart_(globalGridPart),
      indexContainer_(indexContainer),
      boundaryInfoContainer_(boundaryInfoContainer),
      indexSet_(*globalGridPart_, indexContainer_)
  {}

  const IndexSetType& indexSet() const
  {
    return indexSet_;
  }

  const GridType& grid() const
  {
    return globalGridPart_->grid();
  }

  const GlobalGridPartType& globalGridPart() const
  {
    return *globalGridPart_;
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType begin() const
  {
    return typename BaseType::template Codim< codim >::IteratorType(*globalGridPart_, indexContainer_);
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType begin() const
  {
    return typename BaseType::template Codim< codim >::template Partition< pitype >::IteratorType(*globalGridPart_, indexContainer_);
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType end() const
  {
    return typename BaseType::template Codim< codim >::IteratorType(*globalGridPart_, indexContainer_, true);
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType end() const
  {
    return typename BaseType::template Codim< codim >::template Partition< pitype >::IteratorType(*globalGridPart_, indexContainer_, true);
  }

  IntersectionIteratorType ibegin(const EntityType& entity) const
  {
    const IndexType& globalIndex = globalGridPart_->indexSet().index(entity);
    const typename BoundaryInfoContainerType::const_iterator result = boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map< int, int >& info = result->second;
      // return wrapped iterator
      return IntersectionIteratorType(*globalGridPart_, entity, info);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIteratorType(*globalGridPart_, entity);
    } // if this is an entity at the boundary
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& entity) const
  {
    const IndexType& globalIndex = globalGridPart_->indexSet().index(entity);
    const typename BoundaryInfoContainerType::const_iterator result = boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map< int, int >& info = result->second;
      // return wrapped iterator
      return IntersectionIteratorType(*globalGridPart_, entity, info, true);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIteratorType(*globalGridPart_, entity, true);
    } // if this is an entity at the boundary
  }

  int boundaryId(const IntersectionType& intersection) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return intersection.boundaryId();
  }

  int level() const
  {
    return globalGridPart_->level();
  }

  template< class DataHandleImp ,class DataType >
  void communicate(CommDataHandleIF< DataHandleImp, DataType > & data, InterfaceType iftype, CommunicationDirection dir) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    globalGridPart_->communicate(data,iftype,dir);
  }

private:
  const Dune::shared_ptr< const GlobalGridPartType > globalGridPart_;
  const Dune::shared_ptr< const IndexContainerType > indexContainer_;
  const Dune::shared_ptr< const BoundaryInfoContainerType > boundaryInfoContainer_;
  const IndexSetType indexSet_;
}; // class Const

template< class GlobalGridPartImp >
class ConstCoupling;

template< class GlobalGridPartImp >
struct ConstCouplingTraits
  : public ConstTraits< GlobalGridPartImp >
{
  typedef Dune::grid::Part::Interface< typename GlobalGridPartImp::Traits > GlobalGridPartType;

  typedef Dune::grid::Part::Local::IndexBased::ConstCoupling< GlobalGridPartImp > GridPartType;

  //! localized intersection iterator
  typedef Dune::grid::Part::Iterator::Intersection::Local< GlobalGridPartType > IntersectionIteratorType;
}; // class ConstCouplingTraits

template< class GlobalGridPartImp >
class ConstCoupling
  : public Const< GlobalGridPartImp >
{
public:
  typedef ConstCoupling< GlobalGridPartImp > ThisType;

  typedef Dune::grid::Part::Local::IndexBased::ConstCouplingTraits< GlobalGridPartImp > Traits;

  typedef Const< GlobalGridPartImp > BaseType;

  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::GlobalGridPartType GlobalGridPartType;

  typedef typename BaseType::IndexType IndexType;

  typedef typename BaseType::IndexContainerType IndexContainerType;

  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;

  typedef BaseType InsideType;

  typedef BaseType OutsideType;

  //! container type for the intersection information
  typedef std::map< IndexType, std::set< int > > IntersectionInfoContainerType;

  ConstCoupling(const Dune::shared_ptr< const GlobalGridPartType > globalGridPart,
                const Dune::shared_ptr< const IndexContainerType > indexContainer,
                const Dune::shared_ptr< const IntersectionInfoContainerType > intersectionContainer,
                const Dune::shared_ptr< const InsideType > inside,
                const Dune::shared_ptr< const OutsideType > outside)
    : BaseType(globalGridPart, indexContainer, Dune::shared_ptr< const BoundaryInfoContainerType >(new BoundaryInfoContainerType())),
      intersectionContainer_(intersectionContainer),
      inside_(inside),
      outside_(outside)
  {}

  IntersectionIteratorType ibegin(const EntityType& entity) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const std::set< int >& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& entity) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const std::set< int >& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info, true);
  } // IntersectionIteratorType iend(const EntityType& entity) const

  Dune::shared_ptr< const InsideType > inside() const
  {
    return inside_;
  }

  Dune::shared_ptr< const InsideType > outside() const
  {
    return outside_;
  }

private:
  const Dune::shared_ptr< const IntersectionInfoContainerType > intersectionContainer_;
  const Dune::shared_ptr< const InsideType > inside_;
  const Dune::shared_ptr< const OutsideType > outside_;
}; // class ConstCoupling

template< class GlobalGridPartImp >
class ConstBoundary;

template< class GlobalGridPartImp >
struct ConstBoundaryTraits
  : public ConstTraits< GlobalGridPartImp >
{
  typedef Dune::grid::Part::Interface< typename GlobalGridPartImp::Traits > GlobalGridPartType;

  typedef Dune::grid::Part::Local::IndexBased::ConstBoundary< GlobalGridPartImp > GridPartType;

  //! localized intersection iterator
  typedef Dune::grid::Part::Iterator::Intersection::Local< GlobalGridPartType > IntersectionIteratorType;
}; // class ConstBoundaryTraits

template< class GlobalGridPartImp >
class ConstBoundary
  : public Const< GlobalGridPartImp >
{
public:
  typedef ConstBoundary< GlobalGridPartImp > ThisType;

  typedef Dune::grid::Part::Local::IndexBased::ConstBoundaryTraits< GlobalGridPartImp > Traits;

  typedef Const< GlobalGridPartImp > BaseType;

  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::GlobalGridPartType GlobalGridPartType;

  typedef typename BaseType::IndexType IndexType;

  typedef typename BaseType::IndexContainerType IndexContainerType;

  typedef typename BaseType::BoundaryInfoContainerType BoundaryInfoContainerType;

  typedef BaseType InsideType;

  typedef BaseType OutsideType;

  //! container type for the intersection information
  typedef std::map< IndexType, std::set< int > > IntersectionInfoContainerType;

  ConstBoundary(const Dune::shared_ptr< const GlobalGridPartType > globalGridPart,
                const Dune::shared_ptr< const IndexContainerType > indexContainer,
                const Dune::shared_ptr< const IntersectionInfoContainerType > intersectionContainer,
                const Dune::shared_ptr< const InsideType > inside)
    : BaseType(globalGridPart, indexContainer, Dune::shared_ptr< const BoundaryInfoContainerType >(new BoundaryInfoContainerType())),
      intersectionContainer_(intersectionContainer),
      inside_(inside)
  {}

  IntersectionIteratorType ibegin(const EntityType& entity) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const std::set< int >& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& entity) const
  {
    const IndexType& globalIndex = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const std::set< int >& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info, true);
  } // IntersectionIteratorType iend(const EntityType& entity) const

  Dune::shared_ptr< const InsideType > inside() const
  {
    return inside_;
  }

private:
  const Dune::shared_ptr< const IntersectionInfoContainerType > intersectionContainer_;
  const Dune::shared_ptr< const InsideType > inside_;
}; // class ConstBoundary

} // namespace IndexBased

} // namespace Local

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_LOCAL_INDEXBASED_HH
