
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

template< class WrappedIntersectionImp >
class FakeDomainBoundary
{
public:
  typedef WrappedIntersectionImp WrappedIntersectionType;

  typedef FakeDomainBoundary< WrappedIntersectionType > ThisType;

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

  FakeDomainBoundary()
    : bound_(false),
      passThrough_(true)
  {
  }

  void bind(Dune::shared_ptr< const WrappedIntersectionType > wrappedIntersection)
  {
    wrappedIntersection_ = wrappedIntersection;
    bound_ = true;
  }

  bool boundary() const
  {
    //if (passThrough_)
      return wrappedIntersection_->boundary();
  }

  int boundaryId() const
  {
    //if (passThrough_)
      return wrappedIntersection_->boundaryId();
  }

  size_t boundarySegmentIndex() const
  {
    //if (passThrough_)
      return wrappedIntersection_->boundarySegmentIndex();
  }

  bool neighbor() const
  {
    //if (passThrough_)
      return wrappedIntersection_->neighbor();
  }

  EntityPointer inside() const
  {
    //if (passThrough_)
      return wrappedIntersection_->inside();
  }

  EntityPointer outside() const
  {
    //if (passThrough_)
      return wrappedIntersection_->outside();
  }

  bool conforming() const
  {
    //if (passThrough_)
      return wrappedIntersection_->conforming();
  }

  LocalGeometry geometryInInside() const
  {
    //if (passThrough_)
      return wrappedIntersection_->geometryInInside();
  }

  LocalGeometry geometryInOutside() const
  {
    //if (passThrough_)
      return wrappedIntersection_->geometryInOutside();
  }

  Geometry geometry() const
  {
    //if (passThrough_)
      return wrappedIntersection_->geometry();
  }

  Dune::GeometryType type() const
  {
    //if (passThrough_)
      return wrappedIntersection_->type();
  }

  int indexInInside() const
  {
    std::cout << "=====================" << std::endl;
    //if (passThrough_)
      return wrappedIntersection_->indexInInside();
  }

  int indexInOutside() const
  {
    //if (passThrough_)
      return wrappedIntersection_->indexInOutside();
  }

  GlobalCoordinate outerNormal(const LocalCoordinate& local) const
  {
    //if (passThrough_)
      return wrappedIntersection_->outerNormal(local);
  }

  GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
  {
    //if (passThrough_)
      return wrappedIntersection_->integrationOuterNormal(local);
  }

  GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
  {
    //if (passThrough_)
      return wrappedIntersection_->unitOuterNormal(local);
  }

  GlobalCoordinate centerUnitOuterNormal() const
  {
    //if (passThrough_)
      return wrappedIntersection_->centerUnitOuterNormal();
  }

private:
  mutable Dune::shared_ptr< const WrappedIntersectionType > wrappedIntersection_;
  bool bound_;
  bool passThrough_;
}; // class FakeDomainBoundary

} // namespace Wrapper

} // namespace Intersection

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_INTERSECTION_WRAPPER_HH
