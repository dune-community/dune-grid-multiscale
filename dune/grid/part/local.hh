
#ifndef DUNE_GRID_MULTISCALE_GRIDPART_LOCAL_HH
#define DUNE_GRID_MULTISCALE_GRIDPART_LOCAL_HH

// system
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/gridpart/leaf.hh>
#include <dune/grid/multiscale/gridpart/iterator/codim0.hh>
#include <dune/grid/multiscale/gridpart/indexset/local.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace GridPart {

template< class GlobalGridPartImp >
class Local;

template< class GlobalGridPartImp >
struct LocalTraits
  : public Dune::grid::Multiscale::GridPart::LeafTraits< typename GlobalGridPartImp::GridType >
{
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef Dune::grid::Multiscale::GridPart::Local< GlobalGridPartType > GridPartType;

  typedef typename GlobalGridPartType::GridType GridType;

  typedef Dune::grid::Multiscale::GridPart::IndexSet::Local::IndexBased< GridPartType > IndexSetType;

  typedef GlobalGridPartType::GridViewType GridViewType;

  static const PartitionIteratorType indexSetPartitionType = GlobalGridPartType::indexSetPartitionType;

  static const bool conforming = GlobalGridPartType::conforming;

  typedef typename GlobalGridPartType::IntersectionIteratorType IntersectionIteratorType;

  //! rene fragen, wie das hier vernünftig geht (für codim = 0 spezialisieren, für rest boom)
  template< int codim >
  struct Codim
  {
    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename GlobalGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType IteratorType;
//      typedef typename Dune::grid::Multiscale::GridPart::Iterator::Codim0::IndexBased< GridPartType, pitype > IteratorType;
    };
  };
}; // class LocalTraits

//! \attention Only grids with one GeometryType supported!
template< class GlobalGridPartImp >
class Local
  : public Dune::grid::Multiscale::GridPart::Default< LocalTraits< GlobalGridPartImp > >
{
public:
  typedef Local< GlobalGridPartImp > ThisType;

  typedef Dune::grid::Multiscale::GridPart::Default< LocalTraits< GlobalGridPartImp > > BaseType;

  typedef Dune::grid::Multiscale::GridPart::LocalTraits< GlobalGridPartImp > Traits;

  typedef typename Traits::GridType GridType;

  typedef typename Traits::GlobalGridPartType GlobalGridPartType;

  typedef typename Traits::IndexSetType IndexSetType;

  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename Traits::GridViewType GridViewType;

  typedef typename GridType::template Codim<0>::Entity EntityType;

private:
  typedef typename IndexSetType::IndexType IndexType;

public:
  explicit Local(GlobalGridPartType& globalGridPart, Dune::shared_ptr< std::map< IndexType, IndexType > > globalToLocaIndexMap)
    : BaseType(globalGridPart.grid()),
      globalGridPart_(globalGridPart),
      globalToLocaIndexMap_(globalToLocaIndexMap),
      indexSet_(*this)
  {}

  const IndexSetType& indexSet() const
  {
    return indexSet_;
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
  friend class Dune::grid::Multiscale::GridPart::IndexSet::Local::IndexBased< ThisType >;

  GlobalGridPartType& globalGridPart_;
  Dune::shared_ptr< std::map< IndexType, IndexType > > globalToLocaIndexMap_;
  const IndexSetType indexSet_;
}; // class Local

} // namespace GridPart

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_GRIDPART_LOCAL_HH
