#ifndef DUNE_GRID_MULTISCALE_GRIDPART_HH
#define DUNE_GRID_MULTISCALE_GRIDPART_HH

// dune-common
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/deprecated.hh>

// dune-grid
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/datahandleif.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/gridpart/view.hh>
#include <dune/grid/multiscale/gridpart/indexset/default.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace GridPart {

template< class GridImp >
class LeafTraits;

template< class TraitsImp >
class Interface
{
public:
  typedef Interface< TraitsImp > ThisType;

  //! \brief Type of the Traits
  typedef TraitsImp Traits;

  //! \brief Type of the implementation
  typedef typename Traits::GridPartType GridPartType;

  //! \brief type of Grid implementation
  typedef typename Traits::GridType GridType;

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
  const GridType & grid () const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
    return asImp().grid();
  }

  //! \brief Returns reference to the underlying grid
  GridType & grid ()
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().grid()));
    return asImp().grid();
  }

  GridViewType gridView () const
  {
    typedef typename GridViewType :: GridViewImp Impl;
    return GridViewType( Impl( asImp() ) );
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
  typename Codim< codim > :: IteratorType
  begin () const
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
  typename Codim< codim > :: template Partition< pitype > :: IteratorType
  begin () const
  {
    CHECK_INTERFACE_IMPLEMENTATION( (asImp().template begin< codim, pitype >()) );
    return asImp().template begin< codim, pitype >();
  }

  /** \brief obtain end iterator for the interior-border partition
   *
   *  \tparam  codim  codimension for which the iterator is requested
   */
  template< int codim >
  typename Codim< codim > :: IteratorType
  end () const
  {
    CHECK_INTERFACE_IMPLEMENTATION( (asImp().template end< codim >()) );
    return asImp().template end< codim >();
  }

  /** \brief obtain end iterator for the given partition
   *
   *  \tparam  codim   codimension for which the iterator is requested
   *  \tparam  pitype  requested partition iterator type
   */
  template< int codim, PartitionIteratorType pitype >
  typename Codim< codim > :: template Partition< pitype > :: IteratorType
  end () const
  {
    CHECK_INTERFACE_IMPLEMENTATION( (asImp().template end< codim, pitype >()) );
    return asImp().template end< codim, pitype >();
  }

  //! \brief Level of the grid part
  int level() const
  {
    CHECK_INTERFACE_IMPLEMENTATION((asImp().level()));
    return asImp().level();
  }

  //! \brief ibegin of corresponding intersection iterator for given entity
  IntersectionIteratorType
  ibegin ( const typename Codim< 0 >::EntityType &entity ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( (asImp().ibegin( entity )) );
    return asImp().ibegin( entity );
  }

  //! \brief iend of corresponding intersection iterator for given entity
  IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( (asImp().iend( entity )) );
    return asImp().iend( entity );
  }

  int boundaryId ( const IntersectionType &intersection ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().boundaryId( intersection ) );
    return asImp().boundaryId( intersection );
  }

  //! \brief corresponding communication method for grid part
  template <class DataHandleImp,class DataType>
  void communicate(CommDataHandleIF<DataHandleImp,DataType> & data,
                   InterfaceType iftype, CommunicationDirection dir) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().communicate(data,iftype,dir)));
  }

protected:
  //! do not create explicit instances of this class
  Interface () {}

private:
  // Barton-Nackman
  GridPartType& asImp() {
    return static_cast<GridPartType&>(*this);
  }

  // const Barton-Nackman
  const GridPartType& asImp() const {
    return static_cast<const GridPartType&>(*this);
  }
}; // class Interface

//! \brief Default implementation for the GridPart classes
template< class TraitsImp >
class Default
: public Interface< TraitsImp >
{
  typedef Default< TraitsImp > ThisType;
  typedef Interface< TraitsImp > BaseType;

public:
  //! Grid implementation
  typedef typename BaseType :: GridType GridType;
  //! Index set implementation
  typedef typename BaseType :: IndexSetType IndexSetType;

protected:
  GridType &grid_;

protected:
  //! constructor
  Default ( GridType &grid )
  : BaseType(),
    grid_( grid )
  {}

  Default ( const ThisType &other )
  : BaseType(),
    grid_( other.grid_ )
  {}

  ~Default ()
  {}

public:
  //! Returns const reference to the underlying grid
  const GridType &grid () const { return grid_; }

  //! Returns reference to the underlying grid
  GridType &grid () { return grid_; }

  /** \brief obtain begin iterator for the interior-border partition
   *
   *  \tparam  codim  codimension for which the iterator is requested
   */
  template< int codim >
  typename BaseType :: template Codim< codim > :: IteratorType
  begin () const
  {
    return BaseType :: template begin< codim, InteriorBorder_Partition >();
  }

  /** \brief obtain end iterator for the interior-border partition
   *
   *  \tparam  codim  codimension for which the iterator is requested
   */
  template< int codim >
  typename BaseType :: template Codim< codim > :: IteratorType
  end () const
  {
    return BaseType :: template end< codim, InteriorBorder_Partition >();
  }

private:
  template< int codim, PartitionIteratorType pitype >
  typename BaseType :: template Codim< codim > :: template Partition< pitype > :: IteratorType
  begin () const;

  template< int codim, PartitionIteratorType pitype >
  typename BaseType :: template Codim< codim > :: template Partition< pitype > :: IteratorType
  end () const;
}; // class Default

//! \brief Selects the leaf level of a grid
template< class GridImp >
class Leaf
: public Default< LeafTraits< GridImp > >
{
  typedef Leaf< GridImp > ThisType;
  typedef Default< LeafTraits< GridImp > > BaseType;

public:
  //- Public typedefs and enums
  //! Type definitions
  typedef LeafTraits< GridImp > Traits;

  //! Grid implementation type
  typedef typename Traits::GridType GridType;
  //! The leaf index set of the grid implementation
  typedef typename Traits::IndexSetType IndexSetType;

  //! The corresponding IntersectionIterator
  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  //! the leaf grid view from the grid
  typedef typename GridType :: template Partition < All_Partition > :: LeafGridView LeafGridView;

private:
  typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

public:
  //- Public methods
  //! Constructor
  explicit Leaf( GridType &grid )
  : BaseType( grid ),
    leafView_( grid.leafView() ),
    isetWrapper_( grid )
  {}

  //! copy constructor
  Leaf( const ThisType &other )
  : BaseType( other ),
    leafView_( other.leafView_ ),
    isetWrapper_( other.grid() )
  {}

  //! Returns reference to index set of the underlying grid
  const IndexSetType &indexSet () const
  {
    return isetWrapper_;
  }

  //! Begin iterator on the leaf level
  template< int codim >
  typename BaseType :: template Codim< codim > :: IteratorType
  begin () const
  {
    return BaseType :: template begin< codim >();
  }

  //! Begin iterator on the leaf level
  template< int codim, PartitionIteratorType pitype >
  typename Traits :: template Codim< codim > :: template Partition< pitype > :: IteratorType
  begin () const
  {
    return leafView_.template begin< codim, pitype >();
  }

  //! Begin iterator on the leaf level
  template< int codim >
  typename BaseType :: template Codim< codim > :: IteratorType
  end () const
  {
    return BaseType :: template end< codim >();
  }

  //! End iterator on the leaf level
  template< int codim, PartitionIteratorType pitype >
  typename Traits :: template Codim< codim > :: template Partition< pitype > :: IteratorType
  end () const
  {
    return leafView_.template end< codim, pitype >();
  }

  //! ibegin of corresponding intersection iterator for given entity
  IntersectionIteratorType ibegin(const EntityCodim0Type & en) const
  {
    return en.ileafbegin();
  }

  //! iend of corresponding intersection iterator for given entity
  IntersectionIteratorType iend(const EntityCodim0Type & en) const
  {
    return en.ileafend();
  }

  int boundaryId ( const IntersectionType &intersection ) const
  {
    return intersection.boundaryId();
  }

  //! Returns maxlevel of the grid
  int level() const { return this->grid().maxLevel(); }

  //! corresponding communication method for this grid part
  template <class DataHandleImp,class DataType>
  void communicate(CommDataHandleIF<DataHandleImp,DataType> & data,
                   InterfaceType iftype, CommunicationDirection dir) const
  {
    this->grid().communicate(data,iftype,dir);
  }

private:
  //! leaf grid view
  LeafGridView leafView_ ;
  //! GridDefaultIndexSet Wrapper
  IndexSetType isetWrapper_;
}; // class Leaf

//! Type definitions for the LeafGridPart class
template< class GridImp >
struct LeafTraits
{
  /** \brief The type of the grid */
  typedef GridImp GridType;

  /** \brief The type of the corresponding grid part class */
  typedef Leaf< GridImp > GridPartType;

  /** \brief The appropriate index set */
  typedef Dune::grid::Multiscale::GridPart::IndexSet::Default::Leaf<GridType> IndexSetType;

  static const PartitionIteratorType indexSetPartitionType = All_Partition;

  /** \brief The appropriate intersection iterator */
  typedef typename GridType::template Codim<0>::Entity::LeafIntersectionIterator
    IntersectionIteratorType;

  /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
  template< int codim >
  struct Codim
  {
    template< PartitionIteratorType pitype >
    struct Partition
    {
      /** \brief Iterators over the entities of codimension <tt>cd</tt> of this grid part */
      typedef typename GridType :: template Codim< codim >
        :: template Partition< pitype > :: LeafIterator
        IteratorType;
    };
  };

  //! \brief is true if grid on this view only has conforming intersections
  static const bool conforming = Capabilities::isLeafwiseConforming<GridType>::v;
}; // struct LeafTraits

} // namespace GridPart

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_GRIDPART_HH
