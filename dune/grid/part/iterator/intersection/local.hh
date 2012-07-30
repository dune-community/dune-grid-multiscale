
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

template< class GlobalGridPartImp >
class Local
  : public GlobalGridPartImp::IntersectionIteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef Local< GlobalGridPartType > ThisType;

  typedef typename GlobalGridPartType::IntersectionIteratorType BaseType;

  typedef typename GlobalGridPartType::template Codim< 0 >::EntityType EntityType;

  typedef std::map< int, int > InfoContainerType;

private:
  typedef typename BaseType::Intersection BaseIntersectionType;

public:
  typedef Dune::grid::Part::Intersection::Wrapper::FakeDomainBoundary< BaseIntersectionType > Intersection;

  Local(const GlobalGridPartType& globalGridPart,
        const EntityType& entity,
        bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity)),
      passThrough_(true)
  {}

  Local(const GlobalGridPartType& globalGridPart,
        const EntityType& entity,
        const InfoContainerType infoContainer,
        bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity)),
      passThrough_(false),
      infoContainer_(infoContainer)
  {}

//  ~Local()
//  {
//    delete intersection_;
//  }

  const Intersection& operator*() const
  {
    // if we are not on an entity of interest
    if (passThrough_) {
      // create wrapper
      baseIntersection_ = Dune::shared_ptr< const BaseIntersectionType >(BaseType::operator->());
      intersection_.bind(baseIntersection_);
      return intersection_;
    }
  } // const Intersection& operator*() const

//  const Intersection* operator->() const;

private:
  bool passThrough_;
  mutable Dune::shared_ptr< const BaseIntersectionType > baseIntersection_;
  mutable Intersection intersection_;
  const InfoContainerType infoContainer_;
}; // class Local

} // namespace Intersection

} // namespace Iterator

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH
