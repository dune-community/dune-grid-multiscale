// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_PART_INTERSECTION_WRAPPER_HH
#define DUNE_GRID_PART_INTERSECTION_WRAPPER_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// dune-geometry
#include <dune/geometry/type.hh>

// dune-grid
#include <dune/grid/common/intersection.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Intersection {

namespace Wrapper {

template <class IntersectionIteratorImp, class WrappedIntersectionImp>
class FakeDomainBoundary
{
public:
  typedef IntersectionIteratorImp IntersectionIteratorType;

  typedef WrappedIntersectionImp WrappedIntersectionType;

  typedef FakeDomainBoundary<IntersectionIteratorType, WrappedIntersectionType> ThisType;

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
    : intersectionIterator_(intersectionIterator)
    , passThrough_(true)
    , boundaryId_(-1)
  {
  }

  void setPassThrough(const bool passThrough)
  {
    passThrough_ = passThrough;
  }

  void setBoundaryId(const int boundaryId)
  {
    boundaryId_ = boundaryId;
  }

  bool neighbor() const
  {
    if (passThrough_)
      return intersectionIterator_.getBaseIntersection().neighbor();
    else
      return false;
  }

  bool boundary() const
  {
    if (passThrough_)
      return intersectionIterator_.getBaseIntersection().boundary();
    else
      return true;
  }

  int boundaryId() const
  {
    if (passThrough_)
      return intersectionIterator_.getBaseIntersection().boundaryId();
    else
      return boundaryId_;
  }

  size_t boundarySegmentIndex() const
  {
    return intersectionIterator_.getBaseIntersection().boundarySegmentIndex();
  }

  EntityPointer inside() const
  {
    return intersectionIterator_.getBaseIntersection().inside();
  }

  EntityPointer outside() const
  {
    return intersectionIterator_.getBaseIntersection().outside();
  }

  bool conforming() const
  {
    return intersectionIterator_.getBaseIntersection().conforming();
  }

  LocalGeometry geometryInInside() const
  {
    return intersectionIterator_.getBaseIntersection().geometryInInside();
  }

  LocalGeometry geometryInOutside() const
  {
    return intersectionIterator_.getBaseIntersection().geometryInOutside();
  }

  Geometry geometry() const
  {
    return intersectionIterator_.getBaseIntersection().geometry();
  }

  Dune::GeometryType type() const
  {
    return intersectionIterator_.getBaseIntersection().type();
  }

  int indexInInside() const
  {
    return intersectionIterator_.getBaseIntersection().indexInInside();
  }

  int indexInOutside() const
  {
    return intersectionIterator_.getBaseIntersection().indexInOutside();
  }

  GlobalCoordinate outerNormal(const LocalCoordinate& local) const
  {
    return intersectionIterator_.getBaseIntersection().outerNormal(local);
  }

  GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
  {
    return intersectionIterator_.getBaseIntersection().integrationOuterNormal(local);
  }

  GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
  {
    return intersectionIterator_.getBaseIntersection().unitOuterNormal(local);
  }

  GlobalCoordinate centerUnitOuterNormal() const
  {
    return intersectionIterator_.getBaseIntersection().centerUnitOuterNormal();
  }

  const typename IntersectionIteratorType::BaseIntersectionType& asBase() const
  {
    return intersectionIterator_.getBaseIntersection();
  }

private:
  const IntersectionIteratorType& intersectionIterator_;
  bool passThrough_;
  int boundaryId_;
}; // class FakeDomainBoundary

} // namespace Wrapper

} // namespace Intersection

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_INTERSECTION_WRAPPER_HH
