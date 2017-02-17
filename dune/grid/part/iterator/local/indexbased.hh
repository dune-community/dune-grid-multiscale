// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_PART_ITERATOR_CODIM0_HH
#define DUNE_GRID_PART_ITERATOR_CODIM0_HH

// system
#include <map>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

// dune-geometry
#include <dune/geometry/type.hh>

// dune-grid
#include <dune/grid/common/grid.hh>

namespace Dune {
namespace grid {
namespace Part {
namespace Iterator {
namespace Local {

/**
 *  \brief  Iterates over those entities of a grid part, the indices of which match predefined ones.
 *  \todo   Replace GlobalGridPartImp with Interface< GlobalGridPartTraitsImp >!
 *  \todo   Document!
 */
template <class GlobalGridPartImp, int codim, Dune::PartitionIteratorType pitype>
class IndexBased : public GlobalGridPartImp::template Codim<codim>::template Partition<pitype>::IteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef IndexBased<GlobalGridPartType, codim, pitype> ThisType;

  typedef typename GlobalGridPartImp::template Codim<codim>::template Partition<pitype>::IteratorType BaseType;

  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef Dune::GeometryType GeometryType;

private:
  typedef std::map<IndexType, IndexType> IndexMapType;

public:
  typedef std::map<GeometryType, std::map<IndexType, IndexType>> IndexContainerType;

  typedef typename BaseType::Entity Entity;

  IndexBased(const GlobalGridPartType& globalGridPart,
             const Dune::shared_ptr<const IndexContainerType> indexContainer,
             const bool end = false)
    : BaseType(end ? globalGridPart.template end<codim, pitype>() : globalGridPart.template begin<codim, pitype>())
    , globalGridPart_(globalGridPart)
    , indexContainer_(indexContainer)
    , workAtAll_(0)
  {
    if (!end) {
      // loop over all GeometryTypes
      for (typename IndexContainerType::const_iterator iterator = indexContainer_->begin();
           iterator != indexContainer_->end();
           ++iterator) {
        // treat only the codim 0 ones
        if (iterator->first.dim() == (GlobalGridPartType::GridType::dimension - codim)) {
          ++workAtAll_;
          last_.insert(std::make_pair(iterator->first, iterator->second.rbegin()->first));
          end_.insert(std::make_pair(iterator->first, iterator->second.end()));
        } // treat only the codim 0 ones
      } // loop over all GeometryTypes
      forward();
    } // if (!end)
  } // IndexBased

  ThisType& operator++()
  {
    if (workAtAll_ > 0) {
      BaseType::operator++();
      forward();
    } else
      BaseType::operator=(globalGridPart_.template end<codim, pitype>());
    return *this;
  } // ThisType& operator++()

private:
  //! iterates forward until we find the next entity that belongs to the local grid part
  void forward()
  {
    bool found = false;
    while (!found && (workAtAll_ > 0)) {
      const Entity& entity = BaseType::operator*();
      const IndexType& index = globalGridPart_.indexSet().index(entity);
      const GeometryType& geometryType = entity.type();
      typename IndexContainerType::const_iterator indexMap = indexContainer_->find(geometryType);
      if (indexMap != indexContainer_->end()) {
        const typename IndexMapType::const_iterator result = indexMap->second.find(index);
        if ((result != end_.find(geometryType)->second)) {
          found = true;
          if (result->first == last_.find(geometryType)->second)
            --workAtAll_;
        } else
          BaseType::operator++();
      } else
        BaseType::operator++();
    }
  } // void forward()

  const GlobalGridPartType& globalGridPart_;
  const Dune::shared_ptr<const IndexContainerType> indexContainer_;
  unsigned int workAtAll_;
  std::map<GeometryType, IndexType> last_;
  std::map<GeometryType, typename IndexMapType::const_iterator> end_;
}; // class IndexBased

} // namespace Local
} // namespace Iterator
} // namespace Part
} // namespace grid
} // namespace Dune
namespace std {

template <class GlobalGridPartImp, int codim, Dune::PartitionIteratorType pitype>
struct iterator_traits<Dune::grid::Part::Iterator::Local::IndexBased<GlobalGridPartImp, codim, pitype>>
{
  typedef ptrdiff_t difference_type;
  typedef const typename Dune::grid::Part::Iterator::Local::IndexBased<GlobalGridPartImp, codim, pitype>::Entity
      value_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef forward_iterator_tag iterator_category;
};

} // namespace std

#endif // DUNE_GRID_PART_ITERATOR_CODIM0_HH
