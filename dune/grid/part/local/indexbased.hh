
#ifndef DUNE_GRID_PART_LOCAL_INDEXBASED_HH
#define DUNE_GRID_PART_LOCAL_INDEXBASED_HH

// system
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// dune-grid-multiscale
#include <dune/grid/part/interface.hh>
//#include <dune/grid/part/iterator/codim0.hh>
//#include <dune/grid/part/indexset/local.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Local {

namespace IndexBased {

template< class GlobalGridPartImp >
class Const;

template< class GlobalGridPartImp >
struct ConstTraits
{
  typedef Dune::grid::Part::Interface< typename GlobalGridPartImp::Traits > GlobalGridPartType;

  typedef Dune::grid::Part::Local::IndexBased::Const< GlobalGridPartImp > GridPartType;

  typedef typename GlobalGridPartType::GridType GridType;

  typedef typename GlobalGridPartType::IndexSetType IndexSetType;
//  typedef Dune::grid::Multiscale::GridPart::IndexSet::Local::IndexBased< GridPartType > IndexSetType;

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

template< class GlobalGridPartImp >
class Const
  : public Dune::grid::Part::Interface< Dune::grid::Part::Local::IndexBased::ConstTraits< typename GlobalGridPartImp::Traits > >
{
public:
  typedef Const< GlobalGridPartImp > ThisType;

  typedef Dune::grid::Part::Local::IndexBased::ConstTraits< typename GlobalGridPartImp::Traits > Traits;

  typedef Dune::grid::Part::Interface< Traits > BaseType;

  typedef typename Traits::GridType GridType;

  typedef typename Traits::GlobalGridPartType GlobalGridPartType;

  typedef typename Traits::IndexSetType IndexSetType;

  typedef typename Traits::IntersectionIteratorType IntersectionIteratorType ;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename GridType::template Codim<0>::Entity EntityType;

private:
  typedef typename IndexSetType::IndexType IndexType;

public:
  Const(const GlobalGridPartType& globalGridPart/*, Dune::shared_ptr< std::map< IndexType, IndexType > > globalToLocaIndexMap*/)
    : globalGridPart_(globalGridPart)/*,
      globalToLocaIndexMap_(globalToLocaIndexMap),
      indexSet_(*this)*/
  {}

  const IndexSetType& indexSet() const
  {
    return globalGridPart_.indexSet();
  }

  const GridType& grid() const
  {
    return globalGridPart_.grid();
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType begin() const
  {
//    return typename BaseType::template Codim< codim >::IteratorType(*this);
    return globalGridPart_.template begin< codim >();
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType begin() const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return globalGridPart_.template begin< codim, pitype >();
  }

  template< int codim >
  typename BaseType::template Codim< codim >::IteratorType end() const
  {
//    return typename BaseType::template Codim< codim >::IteratorType(*this, true);
    return globalGridPart_.template end< codim >();
  }

  template< int codim, PartitionIteratorType pitype >
  typename Traits::template Codim< codim >::template Partition< pitype >::IteratorType end() const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return globalGridPart_.template end< codim, pitype >();
  }

  IntersectionIteratorType ibegin(const EntityType& en) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return en.ileafbegin();
  }

  IntersectionIteratorType iend(const EntityType& en) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return en.ileafend();
  }

  int boundaryId(const IntersectionType& intersection) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    return intersection.boundaryId();
  }

  int level() const
  {
    return globalGridPart_.level();
  }

  template< class DataHandleImp ,class DataType >
  void communicate(CommDataHandleIF< DataHandleImp, DataType > & data, InterfaceType iftype, CommunicationDirection dir) const
  {
    DUNE_THROW(Dune::NotImplemented, "As long as I am not sure what this does or is used for I will not implement this!");
    globalGridPart_.communicate(data,iftype,dir);
  }

private:
//  friend class Dune::grid::Multiscale::GridPart::Iterator::Codim0::IndexBased< ThisType, InteriorBorder_Partition >;
//  friend class Dune::grid::Multiscale::GridPart::IndexSet::Local::IndexBased< ThisType >;

  const GlobalGridPartType& globalGridPart_;
//  Dune::shared_ptr< std::map< IndexType, IndexType > > globalToLocaIndexMap_;
//  const IndexSetType indexSet_;
}; // class Const

} // namespace IndexBased

} // namespace Local

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_LOCAL_INDEXBASED_HH
