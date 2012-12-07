#ifndef DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH
#define DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// system
#include <map>
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace Iterator {

namespace Intersection {

template< class GlobalGridPartImp >
class Local
  : public GlobalGridPartImp::IntersectionIteratorType
{
public:
  typedef GlobalGridPartImp GlobalGridPartType;

  typedef Local< GlobalGridPartType > ThisType;

  typedef typename GlobalGridPartType::IntersectionIteratorType BaseType;

  typedef typename BaseType::Intersection Intersection;

  typedef typename GlobalGridPartType::template Codim< 0 >::EntityType EntityType;

  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef std::set< int > IndexContainerType;

  Local(const GlobalGridPartType& globalGridPart,
        const EntityType& entity,
        const IndexContainerType& indexContainer,
        const bool end = false)
    : BaseType(end ? globalGridPart.iend(entity) : globalGridPart.ibegin(entity)),
      globalGridPart_(globalGridPart),
      entity_(entity),
      indexContainer_(indexContainer),
      workAtAll_(0)
  {
    if (!end) {
      if (indexContainer_.size() > 0) {
        last_ = *(indexContainer_.rbegin());
        end_ = indexContainer_.end();
        ++workAtAll_;
      }
      forward();
    } // if (!end)
  } // Local

  ThisType& operator++()
  {
    if (workAtAll_ > 0) {
      BaseType::operator++();
      forward();
    } else
      BaseType::operator=(globalGridPart_.iend(entity_));
    return *this;
  } // ThisType& operator++()

private:
  //! iterates forward until we find the next intersection of interest
  void forward()
  {
    bool found = false;
    while (!found && (workAtAll_ > 0)) {
      const Intersection& intersection = BaseType::operator*();
      const int index = intersection.indexInInside();
      typename IndexContainerType::const_iterator result = indexContainer_.find(index);
      if (result != end_) {
        found = true;
        if (*result == last_)
          --workAtAll_;
      } else
        BaseType::operator++();
    } // while (!found && (workAtAll_ > 0))
  } // void forward()

  const GlobalGridPartType& globalGridPart_;
  const EntityType& entity_;
  const IndexContainerType& indexContainer_;
  unsigned int workAtAll_;
  int last_;
  typename IndexContainerType::const_iterator end_;
}; // class Local

} // namespace Intersection

} // namespace Iterator

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_ITERATOR_INTERSECTION_LOCAL_HH
