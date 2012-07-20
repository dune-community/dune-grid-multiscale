#ifndef DUNE_GRID_MULTISCALE_FILTERED_CUBE_HH
#define DUNE_GRID_MULTISCALE_FILTERED_CUBE_HH

// system
#include <sstream>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-fem
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/filteredgrid.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/filtered.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace Filter {

template< class GridImp >
class Cube
  : public FilterDefaultImplementation< DefaultFilterTraits< Cube< GridImp >, Dune::LeafGridPart< GridImp > > >
{
public:
  typedef GridImp GridType;

  typedef Cube< GridType > ThisType;

  typedef Dune::LeafGridPart< GridType > GridPartType;

  typedef DefaultFilterTraits< ThisType, GridPartType > Traits;

  typedef FilterDefaultImplementation< Traits > BaseType;

  typedef typename BaseType::FilterType FilterType;

  typedef typename BaseType::EntityCodim0Type EntityCodim0Type;

  typedef typename BaseType::EntityPointerCodim0Type EntityPointerCodim0Type;

private:
  typedef typename GridType::Traits::template Codim< 0 >::Geometry GeometryType;

  enum{ dim = GeometryType::dimension };

  typedef typename Dune::FieldVector< typename GridType::ctype, dim > FieldVectorType;

public:
  Cube(const FieldVectorType& lowerLeft, const FieldVectorType& upperRight)
    : lowerLeft_(lowerLeft),
      upperRight_(upperRight)
  {}

  Cube(const ThisType& other)
    : lowerLeft_(other.lowerLeft_),
      upperRight_(other.upperRight_)
  {}

  inline bool has0Entity(const EntityPointerCodim0Type& entity) const
  {
    return has0Entity(*entity);
  }

  inline bool has0Entity(const EntityCodim0Type& entity) const
  {
    const FieldVectorType center = entity.geometry().global(entity.geometry().center());
    if (dim == 1) {
      return (lowerLeft_[0] <= center[0]) && (center[0] < upperRight_[0]);
    } else if (dim == 2) {
      return (lowerLeft_[0] <= center[0]) && (center[0] < upperRight_[0])
          && (lowerLeft_[1] <= center[1]) && (center[1] < upperRight_[1]);
    } else if (dim == 3) {
      return (lowerLeft_[0] <= center[0]) && (center[0] < upperRight_[0])
          && (lowerLeft_[1] <= center[1]) && (center[1] < upperRight_[1])
          && (lowerLeft_[2] <= center[2]) && (center[2] < upperRight_[2]);
    } else {
        DUNE_THROW(Dune::InvalidStateException, "Error: only implemented for dimension 1, 2 and 3!");
    }
  }

  //! return what boundary id we have in case of boundary intersection
  //! which is either it.boundary == true or has0Entity (it.ouside()) == false
  //! so here true is a good choice
  template <class IntersectionIteratorType>
  inline bool intersectionBoundary(const IntersectionIteratorType & it) const
  {
    return true;
  }
  //! return what boundary id we have in case of boundary intersection
  //! which is either it.boundary == true or has0Entity (it.ouside()) == false
  template <class IntersectionIteratorType>
  inline int intersectionBoundaryId(const IntersectionIteratorType & it) const
  {
    return 1;
  }

  //! if has0Entity is true then we have an interior entity
  template <class IntersectionIteratorType>
  inline bool intersectionNeighbor(const IntersectionIteratorType & it) const
  {
    return true;
  }

  static ThisType createObject(const GridPartType& gridPart)
  {
    return ThisType(FieldVectorType(0.0), FieldVectorType(1.0));
  }

private:
  const FieldVectorType lowerLeft_;
  const FieldVectorType upperRight_;
}; // end Cube

} // namespace Filter

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_FILTERED_CUBE_HH
