
/**
 *  \attention  This file was copied from dune-fem, please observe the corresponding license.
 *  \todo       Add licensing information.
 **/

#ifndef DUNE_GRID_MULTISCALE_GRIDPART_HH
#define DUNE_GRID_MULTISCALE_GRIDPART_HH

// dune-common
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/deprecated.hh>

// dune-grid
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>

// dune-grid-multiscale
#include <dune/grid/part/interface.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Leaf {

template< class GridImp >
class Const;

//! Type definitions for the GridPart::Leaf class
template< class GridImp >
struct ConstTraits
{
  /** \brief The type of the grid */
  typedef GridImp GridType;

  /** \brief The type of the corresponding grid part class */
  typedef Const< GridImp > GridPartType;

  /** \brief The appropriate index set */
  typedef Dune::grid::Part::IndexSet::Default::Leaf< GridType > IndexSetType;

  static const PartitionIteratorType indexSetPartitionType = All_Partition;

  //! the leaf grid view from the grid
  typedef typename GridType::template Partition< All_Partition >::LeafGridView GridViewType;


  /** \brief The appropriate intersection iterator */
  typedef typename GridType::template Codim< 0 >::Entity::LeafIntersectionIterator IntersectionIteratorType;

  /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
  template< int codim >
  struct Codim
  {
    template< PartitionIteratorType pitype >
    struct Partition
    {
      /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
      typedef typename GridType::template Codim< codim >::template Partition< pitype >::LeafIterator IteratorType;
    };
  };

  //! \brief is true if grid on this view only has conforming intersections
  static const bool conforming = Capabilities::isLeafwiseConforming<GridType>::v;
}; // struct LeafTraits

//! \brief Selects the leaf level of a grid
template< class GridImp >
class Const
: public Dune::grid::Part::Interface< ConstTraits< GridImp > >
{
public:
  typedef Const< GridImp > ThisType;

  //! Type definitions
  typedef ConstTraits< GridImp > Traits;

  typedef Interface< Traits > BaseType;

  //! Grid implementation type
  typedef typename Traits::GridType GridType;

  //! The leaf index set of the grid implementation
  typedef typename Traits::IndexSetType IndexSetType;

  //! The corresponding IntersectionIterator
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

  //! Constructor
  Const(const GridType& grid)
  : grid_(grid),
    leafView_(grid_.leafView()),
    isetWrapper_(grid_)
  {}

  //! copy constructor
  Const(const ThisType &other)
  : grid_(other.grid()),
    leafView_(grid_),
    isetWrapper_(grid_)
  {}

  //! Returns reference to index set of the underlying grid
  const IndexSetType &indexSet() const
  {
    return isetWrapper_;
  }

  //! Begin iterator on the leaf level
  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType begin() const
  {
    return this->template begin< codim, InteriorBorder_Partition >();
  }

  //! Begin iterator on the leaf level
  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType begin() const
  {
    return leafView_.template begin< codim, pitype >();
  }

  //! Begin iterator on the leaf level
  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType end() const
  {
    return this->template end< codim, InteriorBorder_Partition >();
  }

  //! End iterator on the leaf level
  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType end() const
  {
    return leafView_.template end< codim, pitype >();
  }

  //! ibegin of corresponding intersection iterator for given entity
  IntersectionIteratorType ibegin(const EntityType& entity) const
  {
    return entity.ileafbegin();
  }

  //! iend of corresponding intersection iterator for given entity
  IntersectionIteratorType iend(const EntityType& entity) const
  {
    return entity.ileafend();
  }

  const GridType& grid() const
  {
    return grid_;
  }

  int boundaryId(const IntersectionType& intersection) const
  {
    return intersection.boundaryId();
  }

  //! Returns maxlevel of the grid
  int level() const
  {
    return grid_.maxLevel();
  }

  //! corresponding communication method for this grid part
  template< class DataHandleImp, class DataType >
  void communicate(CommDataHandleIF< DataHandleImp, DataType >& data, InterfaceType iftype, CommunicationDirection dir) const
  {
    grid_.communicate(data,iftype,dir);
  }

private:
  const GridType& grid_;
  const typename GridType::LeafGridView leafView_ ;
  const IndexSetType isetWrapper_;
}; // class Const

} // namespace Leaf

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_GRIDPART_HH
