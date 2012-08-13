
#ifndef DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH
#define DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH

// system
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-grid-multiscale
#include <dune/grid/part/intersection/wrapper.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Iterator {

namespace Intersection {

namespace Wrapper {

template< class GlobalGridPartImp >
class FakeDomainBoundary
  : public GlobalGridPartImp::IntersectionIteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef FakeDomainBoundary< GlobalGridPartType > ThisType;

  typedef typename GlobalGridPartType::IntersectionIteratorType BaseType;

  typedef typename GlobalGridPartType::template Codim< 0 >::EntityType EntityType;

  typedef std::map< int, int > InfoContainerType;

private:
  typedef typename BaseType::Intersection BaseIntersectionType;

public:
  typedef Dune::grid::Part::Intersection::Wrapper::FakeDomainBoundary< ThisType, BaseIntersectionType > Intersection;

  FakeDomainBoundary(const GlobalGridPartType& globalGridPart,
        const EntityType& entity,
        bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity)),
      passThrough_(true),
      intersection_(*this)
  {}

  FakeDomainBoundary(const GlobalGridPartType& globalGridPart,
        const EntityType& entity,
        const InfoContainerType infoContainer,
        bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity)),
      passThrough_(false),
      intersection_(*this),
      infoContainer_(infoContainer)
  {}

  const Intersection& operator*() const
  {
    setIntersectionState();
    return intersection_;
  }

  const Intersection* operator->() const
  {
    setIntersectionState();
    return &intersection_;
  }

private:
  friend class Dune::grid::Part::Intersection::Wrapper::FakeDomainBoundary< ThisType, BaseIntersectionType >;

  const BaseIntersectionType& getBaseIntersection() const
  {
    return BaseType::operator*();
  }

  void setIntersectionState() const
  {
    // if we are on an entity of interest
    if (passThrough_) {
      intersection_.setPassThrough(true);
    } else {
      const int intersectionIndex = getBaseIntersection().indexInInside();
      // if this intersection is special
      typename InfoContainerType::const_iterator result = infoContainer_.find(intersectionIndex);
      if (result != infoContainer_.end()) {
        intersection_.setPassThrough(false);
        intersection_.setBoundaryId(result->second);
      } else {
        intersection_.setPassThrough(true);
      } // if this intersection is special
    } // if we are not on an entity of interest
  } // void setIntersectionState() const

  bool passThrough_;
  mutable Intersection intersection_;
  const InfoContainerType infoContainer_;
}; // class FakeDomainBoundary

} // namespace Wrapper

} // namespace Intersection

} // namespace Iterator

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH
