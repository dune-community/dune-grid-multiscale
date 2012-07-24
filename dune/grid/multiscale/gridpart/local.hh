
#ifndef DUNE_GRID_MULTISCALE_GRIDPART_LOCAL_HH
#define DUNE_GRID_MULTISCALE_GRIDPART_LOCAL_HH

// system
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/gridpart/leaf.hh>
#include <dune/grid/multiscale/gridpart/iterator/codim0.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace GridPart {

template< class GridImp >
class Local;

template< class GridImp >
struct LocalTraits
  : public Dune::grid::Multiscale::GridPart::LeafTraits< GridImp >
{
  typedef Dune::grid::Multiscale::GridPart::LeafTraits< GridImp > BaseType;

  typedef typename BaseType::GridType GridType;

  typedef Local< GridType > GridPartType;

  typedef typename BaseType::IndexSetType IndexSetType;

  static const PartitionIteratorType indexSetPartitionType = BaseType::indexSetPartitionType;

  typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;

  //! rene fragen, wie das hier vernünftig geht (für codim = 0 spezialisieren, für rest boom)
  template< int codim >
  struct Codim
  {
    template< PartitionIteratorType pitype >
    struct Partition
    {
//      typedef typename BaseType::template Codim< codim >::template Partition< pitype >::IteratorType IteratorType;
      typedef typename Dune::grid::Multiscale::GridPart::Iterator::Codim0::IndexBased< GridPartType, pitype > IteratorType;
    };
  };

  static const bool conforming = BaseType::conforming;
}; // class LocalTraits

template< class GridImp >
class Local
  : public Dune::grid::Multiscale::GridPart::Default< LocalTraits< GridImp > >
{
public:
  typedef Local< GridImp > ThisType;

  typedef Dune::grid::Multiscale::GridPart::Default< LocalTraits< GridImp > > BaseType;

  typedef LocalTraits< GridImp > Traits;

  typedef typename Traits::GridType GridType;

  typedef typename Traits::IndexSetType IndexSetType;

  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename GridType::template Partition< All_Partition >::LeafGridView LeafGridView;

  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef Dune::grid::Multiscale::GridPart::Leaf< GridType > GlobalGridPartType;

  typedef typename GlobalGridPartType::IndexSetType::IndexType GlobalIndexType;

//  typedef typename Dune::grid::Multiscale::GridPart::Iterator::Codim0::IndexBased< ThisType, GlobalIndexType, InteriorBorder_Partition > LocalIteratorType;

  explicit Local(GlobalGridPartType& globalGridPart, Dune::shared_ptr< std::set< GlobalIndexType > > globalIndicesSet)
    : BaseType(globalGridPart.grid()),
      globalGridPart_(globalGridPart),
      globalIndicesSet_(globalIndicesSet)
  {}

  const IndexSetType& indexSet() const
  {
    return globalGridPart_.indexSet();
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType begin() const
  {
    return typename BaseType::template Codim< codim >::IteratorType(*this);
//    return globalGridPart_.template begin< codim >();
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType begin() const
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return globalGridPart_.template begin< codim, pitype >();
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType end() const
  {
    return typename BaseType::template Codim< codim >::IteratorType(*this, true);
//    return globalGridPart_.template end< codim >();
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType end() const
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return globalGridPart_.template end< codim, pitype >();
  }

  IntersectionIteratorType ibegin(const EntityType& en) const
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return en.ileafbegin();
  }

  IntersectionIteratorType iend(const EntityType& en) const
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return en.ileafend();
  }

  int boundaryId(const IntersectionType& intersection) const
  {
    DUNE_THROW(Dune::NotImplemented, "");
    return intersection.boundaryId();
  }

  int level() const
  {
    return globalGridPart_.level();
  }

  template< class DataHandleImp ,class DataType >
  void communicate(CommDataHandleIF< DataHandleImp, DataType > & data, InterfaceType iftype, CommunicationDirection dir) const
  {
    DUNE_THROW(Dune::NotImplemented, "");
    globalGridPart_.communicate(data,iftype,dir);
  }

private:
  friend class Dune::grid::Multiscale::GridPart::Iterator::Codim0::IndexBased< ThisType, InteriorBorder_Partition >;

  GlobalGridPartType& globalGridPart_;
  Dune::shared_ptr< std::set< GlobalIndexType > > globalIndicesSet_;
}; // class Local

} // namespace GridPart

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_GRIDPART_LOCAL_HH
