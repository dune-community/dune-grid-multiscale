// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_PART_LOCAL_INDEXBASED_HH
#define DUNE_GRID_PART_LOCAL_INDEXBASED_HH

#include <map>
#include <set>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/common/gridpart.hh>
#endif

#include <dune/grid/part/iterator/local/indexbased.hh>
#include <dune/grid/part/iterator/intersection/local.hh>
#include <dune/grid/part/iterator/intersection/wrapper.hh>
#include <dune/grid/part/indexset/local.hh>

namespace Dune {
namespace grid {
namespace Part {
namespace Local {
namespace IndexBased {

template <class GlobalGridPartImp>
class Const;

template <class GlobalGridPartImp>
class ConstTraits
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef Const<GlobalGridPartImp> GridPartType;
  typedef typename GlobalGridPartType::GridType GridType;
  typedef typename IndexSet::Local::IndexBased<GlobalGridPartType> IndexSetType;
  typedef typename GlobalGridPartType::CollectiveCommunicationType CollectiveCommunicationType;
  typedef typename GlobalGridPartType::TwistUtilityType TwistUtilityType;

  static const PartitionIteratorType indexSetPartitionType = GlobalGridPartType::indexSetPartitionType;
  static const InterfaceType indexSetInterfaceType         = GlobalGridPartType::indexSetInterfaceType;

  typedef Iterator::Intersection::Wrapper::FakeDomainBoundary<GlobalGridPartType> IntersectionIteratorType;

  template <int codim>
  struct Codim : public GlobalGridPartType::template Codim<codim>
  {
    template <PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename Iterator::Local::IndexBased<GlobalGridPartType, codim, pitype> IteratorType;
    };
  };

  static const bool conforming = GlobalGridPartType::Traits::conforming;
}; // class ConstTraits

/**
 *  \todo Wrap entity, so that entity.ileaf{begin,end}() returns i{begin,end}(entity)
 *  \todo Implement boundaryId(intersection) by adding a std::map< intersectionIndex, boundaryId >!
 */
template <class GlobalGridPartImp>
class Const
#if HAVE_DUNE_FEM
    : public Fem::GridPartInterface<ConstTraits<GlobalGridPartImp>>
#endif
{
public:
  typedef Const<GlobalGridPartImp> ThisType;
  typedef ConstTraits<GlobalGridPartImp> Traits;

private:
#if HAVE_DUNE_FEM
  typedef Fem::GridPartInterface<ConstTraits<GlobalGridPartImp>> BaseType;
  typedef BaseType BaseTraits;
#else
  typedef Traits BaseTraits;
#endif
public:
  typedef typename Traits::GridType GridType;
  typedef typename Traits::CollectiveCommunicationType CollectiveCommunicationType;
  typedef typename Traits::GlobalGridPartType GlobalGridPartType;
  typedef typename Traits::IndexSetType IndexSetType;
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef typename IndexSetType::IndexType IndexType;
  typedef std::map<IndexType, IndexType> IndexMapType;
  typedef Dune::GeometryType GeometryType;
  //! container type for the indices
  typedef std::map<GeometryType, IndexMapType> IndexContainerType;
  //! container type for the boundary information
  typedef std::map<IndexType, std::map<int, int>> BoundaryInfoContainerType;

  Const(const std::shared_ptr<const GlobalGridPartType> globalGridPart,
        const std::shared_ptr<const IndexContainerType> indexContainer,
        const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer)
    : globalGridPart_(globalGridPart)
    , indexContainer_(indexContainer)
    , boundaryInfoContainer_(boundaryInfoContainer)
    , indexSet_(*globalGridPart_, indexContainer_)
  {
  }

  Const(const ThisType& other) = default;

  Const(ThisType&& source) = default;

  const IndexSetType& indexSet() const { return indexSet_; }

  const GridType& grid() const { return globalGridPart_->grid(); }

  const GlobalGridPartType& globalGridPart() const { return *globalGridPart_; }

  template <int codim>
  typename BaseTraits::template Codim<codim>::IteratorType begin() const
  {
    return typename BaseTraits::template Codim<codim>::IteratorType(*globalGridPart_, indexContainer_);
  }

  template <int codim, PartitionIteratorType pitype>
  typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType begin() const
  {
    return typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType(*globalGridPart_,
                                                                                                indexContainer_);
  }

  template <int codim>
  typename BaseTraits::template Codim<codim>::IteratorType end() const
  {
    return typename BaseTraits::template Codim<codim>::IteratorType(*globalGridPart_, indexContainer_, true);
  }

  template <int codim, PartitionIteratorType pitype>
  typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType end() const
  {
    return typename BaseTraits::template Codim<codim>::template Partition<pitype>::IteratorType(
        *globalGridPart_, indexContainer_, true);
  }

  IntersectionIteratorType ibegin(const EntityType& entity) const
  {
    const IndexType& globalIndex                                    = globalGridPart_->indexSet().index(entity);
    const typename BoundaryInfoContainerType::const_iterator result = boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map<int, int>& info = result->second;
      // return wrapped iterator
      return IntersectionIteratorType(*globalGridPart_, entity, info);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIteratorType(*globalGridPart_, entity);
    } // if this is an entity at the boundary
  }   // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& entity) const
  {
    const IndexType& globalIndex                                    = globalGridPart_->indexSet().index(entity);
    const typename BoundaryInfoContainerType::const_iterator result = boundaryInfoContainer_->find(globalIndex);
    // if this is an entity at the boundary
    if (result != boundaryInfoContainer_->end()) {
      // get the information for this entity
      const std::map<int, int>& info = result->second;
      // return wrapped iterator
      return IntersectionIteratorType(*globalGridPart_, entity, info, true);
    } else {
      // return iterator which just passes everything thrugh
      return IntersectionIteratorType(*globalGridPart_, entity, true);
    } // if this is an entity at the boundary
  }

  int boundaryId(const IntersectionType& intersection) const
  {
    DUNE_THROW(Dune::NotImplemented, "Call intersection.boundaryId() instead!");
    return intersection.boundaryId();
  }

  int level() const { return globalGridPart_->level(); }

  template <class DataHandleImp, class DataType>
  void communicate(CommDataHandleIF<DataHandleImp, DataType>& /*data*/, InterfaceType /*iftype*/,
                   CommunicationDirection /*dir*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "As long as I am not sure what this does or is used for I will not implement this!");
    //    globalGridPart_->communicate(data,iftype,dir);
  }

  const CollectiveCommunicationType& comm() const { return grid().comm(); }

private:
  const std::shared_ptr<const GlobalGridPartType> globalGridPart_;
  const std::shared_ptr<const IndexContainerType> indexContainer_;
  const std::shared_ptr<const BoundaryInfoContainerType> boundaryInfoContainer_;
  const IndexSetType indexSet_;
}; // class Const

template <class GlobalGridPartImp>
class ConstCoupling;

template <class GlobalGridPartImp>
struct ConstCouplingTraits : public ConstTraits<GlobalGridPartImp>
{
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef Dune::grid::Part::Local::IndexBased::ConstCoupling<GlobalGridPartImp> GridPartType;

  //! localized intersection iterator
  typedef Dune::grid::Part::Iterator::Intersection::Local<GlobalGridPartType> IntersectionIteratorType;
}; // class ConstCouplingTraits

template <class GlobalGridPartImp>
class ConstCoupling : public Const<GlobalGridPartImp>
{
public:
  typedef ConstCoupling<GlobalGridPartImp> ThisType;

  typedef Dune::grid::Part::Local::IndexBased::ConstCouplingTraits<GlobalGridPartImp> Traits;

  typedef Const<GlobalGridPartImp> BaseType;

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
  typedef std::map<IndexType, std::vector<int>> IntersectionInfoContainerType;

  ConstCoupling(const std::shared_ptr<const GlobalGridPartType> globalGridPart,
                const std::shared_ptr<const IndexContainerType> indexContainer,
                const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer,
                const std::shared_ptr<const InsideType> inside, const std::shared_ptr<const OutsideType> outside)
    : BaseType(globalGridPart, indexContainer,
               std::shared_ptr<const BoundaryInfoContainerType>(new BoundaryInfoContainerType()))
    , intersectionContainer_(intersectionContainer)
    , inside_(inside)
    , outside_(outside)
  {
  }

  ConstCoupling(const ThisType& other) = default;

  ConstCoupling(ThisType&& source) = default;

  IntersectionIteratorType ibegin(const EntityType& entity) const
  {
    const IndexType& globalIndex                                        = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& entity) const
  {
    const IndexType& globalIndex                                        = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info, true);
  } // IntersectionIteratorType iend(const EntityType& entity) const

  std::shared_ptr<const InsideType> inside() const { return inside_; }

  std::shared_ptr<const InsideType> outside() const { return outside_; }

private:
  const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer_;
  const std::shared_ptr<const InsideType> inside_;
  const std::shared_ptr<const OutsideType> outside_;
}; // class ConstCoupling

template <class GlobalGridPartImp>
class ConstBoundary;

template <class GlobalGridPartImp>
struct ConstBoundaryTraits : public ConstTraits<GlobalGridPartImp>
{
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef Dune::grid::Part::Local::IndexBased::ConstBoundary<GlobalGridPartImp> GridPartType;

  //! localized intersection iterator
  typedef Dune::grid::Part::Iterator::Intersection::Local<GlobalGridPartType> IntersectionIteratorType;
}; // class ConstBoundaryTraits

template <class GlobalGridPartImp>
class ConstBoundary : public Const<GlobalGridPartImp>
{
public:
  typedef ConstBoundary<GlobalGridPartImp> ThisType;

  typedef Dune::grid::Part::Local::IndexBased::ConstBoundaryTraits<GlobalGridPartImp> Traits;

  typedef Const<GlobalGridPartImp> BaseType;

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
  typedef std::map<IndexType, std::vector<int>> IntersectionInfoContainerType;

  ConstBoundary(const std::shared_ptr<const GlobalGridPartType> globalGridPart,
                const std::shared_ptr<const IndexContainerType> indexContainer,
                const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer,
                const std::shared_ptr<const InsideType> inside)
    : BaseType(globalGridPart, indexContainer,
               std::shared_ptr<const BoundaryInfoContainerType>(new BoundaryInfoContainerType()))
    , intersectionContainer_(intersectionContainer)
    , inside_(inside)
  {
  }

  ConstBoundary(const ThisType& other) = default;

  ConstBoundary(ThisType&& source) = default;

  IntersectionIteratorType ibegin(const EntityType& entity) const
  {
    const IndexType& globalIndex                                        = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info);
  } // IntersectionIteratorType ibegin(const EntityType& entity) const

  IntersectionIteratorType iend(const EntityType& entity) const
  {
    const IndexType& globalIndex                                        = BaseType::globalGridPart().indexSet().index(entity);
    const typename IntersectionInfoContainerType::const_iterator result = intersectionContainer_->find(globalIndex);
    assert(result != intersectionContainer_->end());
    // get the information for this entity
    const auto& info = result->second;
    // return localized iterator
    return IntersectionIteratorType(BaseType::globalGridPart(), entity, info, true);
  } // IntersectionIteratorType iend(const EntityType& entity) const

  std::shared_ptr<const InsideType> inside() const { return inside_; }

private:
  const std::shared_ptr<const IntersectionInfoContainerType> intersectionContainer_;
  const std::shared_ptr<const InsideType> inside_;
}; // class ConstBoundary

} // namespace IndexBased
} // namespace Local
} // namespace Part
} // namespace grid

#if HAVE_DUNE_FEM

namespace Fem {
namespace GridPartCapabilities {

template <class GridPartType>
struct hasGrid<grid::Part::Local::IndexBased::Const<GridPartType>>
{
  static const bool v = hasGrid<GridPartType>::v;
};

template <class GridPartType>
struct hasSingleGeometryType<grid::Part::Local::IndexBased::Const<GridPartType>>
{
  static const bool v                  = hasSingleGeometryType<GridPartType>::v;
  static const unsigned int topologyId = hasSingleGeometryType<GridPartType>::topologyId;
};

template <class GridPartType>
struct isCartesian<grid::Part::Local::IndexBased::Const<GridPartType>>
{
  static const bool v = isCartesian<GridPartType>::v;
};

template <class GridPartType, int codim>
struct hasEntity<grid::Part::Local::IndexBased::Const<GridPartType>, codim>
{
  static const bool v = hasEntity<GridPartType, codim>::v;
};

template <class GridPartType>
struct isParallel<grid::Part::Local::IndexBased::Const<GridPartType>>
{
  static const bool v = false;
};

template <class GridPartType, int codim>
struct canCommunicate<grid::Part::Local::IndexBased::Const<GridPartType>, codim>
{
  static const bool v = false;
};

template <class GridPartType>
struct isConforming<grid::Part::Local::IndexBased::Const<GridPartType>>
{
  static const bool v = false;
};

} // namespace GridPartCapabilities
} // namespace Fem

#endif // HAVE_DUNE_FEM

} // namespace Dune

#endif // DUNE_GRID_PART_LOCAL_INDEXBASED_HH
