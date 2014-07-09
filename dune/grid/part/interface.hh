// This file was copied from dune-fem and slightly modified.
// See http://dune.mathematik.uni-freiburg.de/ for the orginal
// file, authors, licence and readme!

#ifndef DUNE_GRID_PART_INTERFACE_HH
#define DUNE_GRID_PART_INTERFACE_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// dune-common
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/deprecated.hh>

// dune-grid
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>

// dune-grid-multiscale
#include <dune/grid/part/view.hh>
#include <dune/grid/part/indexset/default.hh>

namespace Dune {

namespace grid {

namespace Part {

template< class TraitsImp >
class Interface
{
public:
  //! \brief Type of the Traits
  typedef TraitsImp Traits;

  typedef Interface< Traits > ThisType;

  //! \brief Type of the implementation
  typedef typename Traits::GridPartType GridPartType;

  //! \brief type of Grid implementation
  typedef typename Traits::GridType GridType;

  typedef typename Traits::CollectiveCommunicationType CollectiveCommunicationType;

  //! \brief Index set implementation
  typedef typename Traits::IndexSetType IndexSetType;

  //! \brief Maximum Partition type, the index set provides indices for
  static const PartitionIteratorType indexSetPartitionType = Traits::indexSetPartitionType;

  //! \brief type of IntersectionIterator
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType;

  //! \brief type of Intersection
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  //! \brief is true if grid on this view only has conforming intersections
  static const bool conforming = Traits::conforming;

  typedef GridView< ViewTraits< GridPartType > > GridViewType;

  typedef typename GridType::ctype ctype;

  static const int dimension = GridType::dimension;

  template< int codim >
  struct Codim
  {
    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType
        IteratorType;
    };

    typedef typename Partition< InteriorBorder_Partition >::IteratorType IteratorType;

    typedef typename GridType::template Codim< codim >::Entity EntityType;
  };

  //! \brief Returns const reference to the underlying grid
  const GridType& grid() const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
    return asImp().grid();
  }

//  //! \brief Returns reference to the underlying grid
//  GridType& grid()
//  {
//    CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
//    return asImp().grid();
//  }

  GridViewType gridView() const
  {
    typedef typename GridViewType::GridViewImp Impl;
    return GridViewType(Impl(asImp()));
  }

  //! \brief Returns reference to index set of the underlying grid
  const IndexSetType& indexSet() const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().indexSet()));
    return asImp().indexSet();
  }

  /** \brief obtain begin iterator for the interior-border partition
   *
   *  \tparam  codim  codimension for which the iterator is requested
   */
  template< int codim >
  typename Codim< codim >::IteratorType begin() const
  {
    CHECK_INTERFACE_IMPLEMENTATION( (asImp().template begin< codim >()) );
    return asImp().template begin< codim >();
  }

  /** \brief obtain begin iterator for the given partition
   *
   *  \tparam  codim   codimension for which the iterator is requested
   *  \tparam  pitype  requested partition iterator type
   */
  template< int codim, PartitionIteratorType pitype >
  typename Codim< codim >::template Partition< pitype >::IteratorType begin() const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().template begin< codim, pitype >()));
    return asImp().template begin< codim, pitype >();
  }

  /** \brief obtain end iterator for the interior-border partition
   *
   *  \tparam  codim  codimension for which the iterator is requested
   */
  template< int codim >
  typename Codim< codim >::IteratorType end() const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().template end< codim >()));
    return asImp().template end< codim >();
  }

  /** \brief obtain end iterator for the given partition
   *
   *  \tparam  codim   codimension for which the iterator is requested
   *  \tparam  pitype  requested partition iterator type
   */
  template< int codim, PartitionIteratorType pitype >
  typename Codim< codim >::template Partition< pitype >::IteratorType end() const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().template end< codim, pitype >()));
    return asImp().template end< codim, pitype >();
  }

  //! \brief Level of the grid part
  int level() const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().level()));
    return asImp().level();
  }

  //! \brief ibegin of corresponding intersection iterator for given entity
  IntersectionIteratorType ibegin(const typename Codim< 0 >::EntityType& entity) const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().ibegin(entity)));
    return asImp().ibegin(entity);
  }

  //! \brief iend of corresponding intersection iterator for given entity
  IntersectionIteratorType iend(const typename Codim< 0 >::EntityType& entity) const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().iend(entity)));
    return asImp().iend(entity);
  }

  int boundaryId(const IntersectionType &intersection) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().boundaryId(intersection));
    return asImp().boundaryId(intersection);
  }

  const CollectiveCommunicationType& comm() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().comm());
    return asImp().comm();
  }

  //! \brief corresponding communication method for grid part
  template< class DataHandleImp, class DataType >
  void communicate(CommDataHandleIF< DataHandleImp, DataType >& data, InterfaceType iftype, CommunicationDirection dir) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().communicate(data,iftype,dir)));
  }

protected:
  Interface()
  {}

private:
  // Barton-Nackman
  GridPartType& asImp()
  {
    return static_cast< GridPartType& >(*this);
  }

  // const Barton-Nackman
  const GridPartType& asImp() const
  {
    return static_cast< const GridPartType& >(*this);
  }
}; // class Interface

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_INTERFACE_HH
