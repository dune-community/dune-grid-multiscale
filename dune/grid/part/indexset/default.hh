
/**
 *  \attention  This file was copied from dune-fem, please observe the corresponding license.
 *  \todo       Add licensing information.
 **/

#ifndef DUNE_GRID_PART_INDEXSET_DEFAULT_HH
#define DUNE_GRID_PART_INDEXSET_DEFAULT_HH

// system
#include <vector>
#include <rpc/rpc.h>

// dune-common
#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>

// dune-grid
#include <dune/grid/common/grid.hh>
//#include <dune/grid/common/adaptcallback.hh> // for compatibility only

// dune-grid-multiscale
#include <dune/grid/part/indexset/empty.hh>

//#include <dune/fem/misc/capabilities.hh>
//#include <dune/fem/misc/mpimanager.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace IndexSet {

namespace Default {

//! Wraps the interface methods of indexsets and adds the addiotnal needed
//! functions
template <class IndexSetImp >
class Wrapper : public Empty
{
public:
  typedef typename IndexSetImp :: IndexType  IndexType ;
  //! The types of the iterator
  template<int cd>
  struct Codim
  {
    template<PartitionIteratorType pitype>
    struct Partition
    {
      typedef typename IndexSetImp::template Codim<cd>::
        template Partition<pitype>::Iterator  Iterator;
    };
  };

  //! store const reference to set
  Wrapper(const IndexSetImp & set, bool adaptive = false)
    : Empty()
    , set_(set)
  {}

  //! store const reference to set
  Wrapper(const Wrapper<IndexSetImp> & s)
    : Empty(/*s.adaptive_*/),
      set_(s.set_)
  {}

  //! return persistent status
  bool persistent () const { return false; }

  //! return size of set for codim
  int size ( GeometryType type ) const
  {
    return set_.size(type);
  }

  //! return size of grid entities per level and codim
  int size ( int codim ) const
  {
    return set_.size(codim);
  }

  //! return index of en
  template <class EntityType>
  int index (const EntityType & en) const
  {
    return set_.index(en);
  }

  //! return sub index of given entities sub entity with codim and number
  template< int codim, class EntityType >
  int DUNE_DEPRECATED subIndex ( const EntityType &entity, int num ) const
  {
    return set_.template subIndex< codim >( entity, num );
  }

  template< class Entity >
  int subIndex ( const Entity &entity, int num, unsigned int codim ) const
  {
    return set_.subIndex( entity, num, codim );
  }

  //! wrap geomTypes method of set
  const std::vector< GeometryType > & geomTypes (int codim) const
  {
    return set_.geomTypes(codim);
  }

  //! returns true if this set provides an index for given entity
  template<class EntityType>
  bool contains (const EntityType& en) const
  {
    return set_.contains(en);
  }

  //! return index
  template< int codim, class Entity >
  int DUNE_DEPRECATED index ( const Entity &entity, int i ) const
  {
    static const int cc = Entity::codimension;
    return set_.subIndex( entity, i, codim );
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   */
  template<int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin () const
  {
    return set_.template begin<cd,pitype> ();
  }

  /** @brief Iterator to one past the last entity of given codim for partition type
   */
  template<int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end () const
  {
    return set_.template end<cd,pitype> ();
  }

private:
  const IndexSetImp & set_;
}; // class Wrapper

//! Wraps LeafIndexSet of Dune Grids for use with LagrangeFunctionSpace
template <class GridType>
class Leaf
:  public Wrapper<typename GridType :: Traits :: LeafIndexSet>
{
  // my type, to be revised
  enum { myType = 5 };

  // my index set type
  typedef typename GridType :: Traits :: LeafIndexSet IndexSetType;
  typedef Leaf<GridType> ThisType;
public:
  //! number of codimensions
  enum { ncodim = GridType::dimension + 1 };
  //! constructor
  Leaf ( const GridType & grid , const int level =-1 )
    : Wrapper < IndexSetType > (grid.leafIndexSet()) {}

  //! constructor taking grid part
  template <class GridPartType>
  Leaf ( const GridPartType& gridPart )
    : Wrapper < IndexSetType > (gridPart.grid().leafIndexSet()) {}

  //! return type (for Grape In/Output)
  static int type() { return myType; }
  //! returns reference to singleton
  static ThisType & instance (const GridType & grid)
  {
    static ThisType set(grid,grid.maxLevel());
    return set;
  }
}; // class Leaf

} // namespace Default

} // namespace IndexSet

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_INDEXSET_DEFAULT_HH
