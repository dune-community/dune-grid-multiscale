// This file was copied from dune-fem and slightly modified.
// See http://dune.mathematik.uni-freiburg.de/ for the orginal
// file, authors, licence and readme!

#ifndef DUNE_GRID_PART_VIEW_HH
#define DUNE_GRID_PART_VIEW_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// dune-grid
#include <dune/grid/common/gridview.hh>

namespace Dune {

namespace grid {

namespace Part {

template< class GridPart >
class ViewImp;

template< class GridPart >
struct ViewTraits
{
  typedef ViewImp< GridPart > GridViewImp;

  typedef typename GridPart::GridType Grid;
  typedef typename GridPart::IndexSetType IndexSet;
  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;

  typedef typename IntersectionIterator :: Intersection Intersection;

  typedef typename Grid :: Traits :: CollectiveCommunication CollectiveCommunication;

  template< int codim >
  struct Codim
  : public Grid :: Traits :: template Codim< codim >
  {
    typedef typename GridPart :: template Codim< codim > :: IteratorType Iterator;

    typedef typename Grid :: template Codim< codim > :: Entity Entity;
    typedef typename Grid :: template Codim< codim > :: EntityPointer
      EntityPointer;

    typedef typename Grid :: template Codim< codim > :: Geometry Geometry;
    typedef typename Grid :: template Codim< codim > :: LocalGeometry
      LocalGeometry;

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename GridPart :: template Codim< codim >
          :: template Partition< pitype > :: IteratorType
          Iterator;
      };
  };

  static const bool conforming = GridPart :: conforming;
}; // struct ViewTraits

template< class GridPart >
class ViewImp
{
  typedef ViewImp< GridPart > ThisType;

public:
  typedef GridPart GridPartType;

  typedef ViewTraits< GridPartType > Traits;

  /** \brief type of the grid */
  typedef typename Traits :: Grid Grid;

  /** \brief type of the index set */
  typedef typename Traits :: IndexSet IndexSet;

  /** \brief type of the intersection */
  typedef typename Traits :: Intersection Intersection;

  /** \brief type of the intersection iterator */
  typedef typename Traits :: IntersectionIterator IntersectionIterator;

  /** \brief type of the collective communication */
  typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

  /** \brief Codim Structure */
  template< int codim >
  struct Codim
  : public Traits :: template Codim< codim >
  {};

  enum { conforming = Traits :: conforming };

  enum { dimension = Grid :: dimension };
  enum { dimensionworld = Grid :: dimensionworld };

private:
  const GridPartType &gridPart_;

public:
  explicit ViewImp ( const GridPartType &gridPart )
  : gridPart_( gridPart )
  {}

  ViewImp ( const ThisType &other )
  : gridPart_( other.gridPart_ )
  {}

private:
  ThisType &operator= ( const ThisType & );

public:
  const Grid &grid () const
  {
    return gridPart_.grid();
  }

  const IndexSet &indexSet () const
  {
    return gridPart_.indexSet();
  }

  int size ( int codim ) const
  {
    return indexSet().size( codim );
  }

  int size ( const GeometryType &type ) const
  {
    return indexSet().size( type );
  }

  template< int codim >
  typename Codim< codim > :: Iterator begin () const
  {
    return gridPart_.template begin< codim >();
  }

  template< int codim, PartitionIteratorType pitype >
  typename Codim< codim > :: template Partition< pitype > :: Iterator begin () const
  {
    return gridPart_.template begin< codim, pitype >();
  }

  template< int codim >
  typename Codim< codim > :: Iterator end () const
  {
    return gridPart_.template end< codim >();
  }

  template< int codim, PartitionIteratorType pitype >
  typename Codim< codim > :: template Partition< pitype > :: Iterator end () const
  {
    return gridPart_.template end< codim, pitype >();
  }

  IntersectionIterator ibegin ( const typename Codim< 0 > :: Entity &entity ) const
  {
    return gridPart_.ibegin( entity );
  }

  IntersectionIterator iend ( const typename Codim< 0 > :: Entity &entity ) const
  {
    return gridPart_.iend( entity );
  }

  const CollectiveCommunication &comm () const
  {
    return grid().comm();
  }

  template< class DataHandleImp, class DataType >
  void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                     InterfaceType iftype,
                     CommunicationDirection dir ) const
  {
    gridPart_.communicate( data, iftype, dir );
  }
};

template< class GridPart >
class View
: public GridView< ViewTraits< GridPart > >
{
  typedef View< GridPart > ThisType;
  typedef GridView< ViewTraits< GridPart > > BaseType;

  typedef typename BaseType :: GridViewImp GridViewImp;

public:
  explicit View ( const GridPart &gridPart )
  : BaseType( GridViewImp( gridPart ) )
  {}

  View ( const ThisType &other )
  : BaseType( other )
  {}
};

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_VIEW_HH
