
#ifndef DUNE_GRID_PART_INTERSECTION_WRAPPER_HH
#define DUNE_GRID_PART_INTERSECTION_WRAPPER_HH

// dune-geometry
#include <dune/geometry/type.hh>

// dune-grid
#include <dune/grid/common/intersection.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Intersection {

namespace Wrapper {

template< class IntersectionIteratorImp, class WrappedIntersectionImp >
class FakeDomainBoundary
{
public:
  typedef IntersectionIteratorImp IntersectionIteratorType;

  typedef WrappedIntersectionImp WrappedIntersectionType;

  typedef FakeDomainBoundary< IntersectionIteratorType, WrappedIntersectionType > ThisType;

  static const int dimension = WrappedIntersectionType::dimension;

  static const int dimensionworld = WrappedIntersectionType::dimensionworld;

  typedef typename WrappedIntersectionType::Entity Entity;

  typedef typename WrappedIntersectionType::EntityPointer EntityPointer;

  typedef typename WrappedIntersectionType::Geometry Geometry;

  typedef typename WrappedIntersectionType::LocalCoordinate LocalCoordinate;

  typedef typename WrappedIntersectionType::GlobalCoordinate GlobalCoordinate;

  typedef typename WrappedIntersectionType::LocalGeometry LocalGeometry;

//  typedef Dune::Intersection<const GridImp, Dune::SIntersectionIterator> Intersection;

  typedef typename WrappedIntersectionType::ctype ctype;

  FakeDomainBoundary(const IntersectionIteratorType& intersectionIterator)
    : intersectionIterator_(intersectionIterator),
      passThrough_(true)
  {
  }

//  void bind(Dune::shared_ptr< const WrappedIntersectionType > wrappedIntersection)
//  {
//    wrappedIntersection_ = wrappedIntersection;
//    bound_ = true;
//  }

  bool boundary() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().boundary();
  }

  int boundaryId() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().boundaryId();
  }

  size_t boundarySegmentIndex() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().boundarySegmentIndex();
  }

  bool neighbor() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().neighbor();
  }

  EntityPointer inside() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().inside();
  }

  EntityPointer outside() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().outside();
  }

  bool conforming() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().conforming();
  }

  LocalGeometry geometryInInside() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().geometryInInside();
  }

  LocalGeometry geometryInOutside() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().geometryInOutside();
  }

  Geometry geometry() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().geometry();
  }

  Dune::GeometryType type() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().type();
  }

  int indexInInside() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().indexInInside();
  }

  int indexInOutside() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().indexInOutside();
  }

  GlobalCoordinate outerNormal(const LocalCoordinate& local) const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().outerNormal(local);
  }

  GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().integrationOuterNormal(local);
  }

  GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().unitOuterNormal(local);
  }

  GlobalCoordinate centerUnitOuterNormal() const
  {
    //if (passThrough_)
      return intersectionIterator_.getBaseIntersection().centerUnitOuterNormal();
  }

private:
  const IntersectionIteratorType& intersectionIterator_;
  bool passThrough_;
}; // class FakeDomainBoundary

} // namespace Wrapper

} // namespace Intersection

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_INTERSECTION_WRAPPER_HH
