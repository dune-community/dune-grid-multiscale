// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_DEFAULT_HH

#include <vector>
#include <set>
#include <map>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/leafgridpart.hh>
#endif

#include <dune/xt/common/color.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/grid/part/local/indexbased.hh>

namespace Dune {
namespace grid {
namespace Multiscale {

#if HAVE_DUNE_FEM

/**
 *  \brief      Default implementation of a Dune::grid::Multiscale::Interface
 *
 *              The basic assumption is, that each entity of the global grid part is added to at most one local grid
 *              part. Overlapping or the like is then handled afterwards by adding additional overlapping regions to
 *              an existing local grid part (not yet implemented).
 *  \attention  Works only for one GeometryType per Codim (I think)! Problem is, that indices are unique per
 * GeometryType, not per Codim.
 *  \todo       Resolve the above Issue (should be easy, compute the size as a sum over all GeometryTypes for a given
 * codim)!
 *
 *  \todo       Giving the local and global gridparts as shared pointers is quite misleading, since it is not
 *              guaranteed that the underlying grid parts and the grid will exist forever. So we should change those to
 *              reference imho.
 */
template <class GridImp>
class Default
{
public:
  typedef GridImp GridType;

  typedef Default<GridType> ThisType;

  static const unsigned int dim = GridType::dimension;
  static const unsigned int dimension = GridType::dimension;

  typedef typename GridType::ctype ctype;

  typedef Fem::LeafGridPart<GridType> GlobalGridPartType;

//  typedef typename GlobalGridPartType::GridViewType GlobalGridViewType;

  typedef Dune::grid::Part::Local::IndexBased::Const<GlobalGridPartType> LocalGridPartType;

//  typedef typename LocalGridPartType::GridViewType LocalGridViewType;

  typedef Dune::grid::Part::Local::IndexBased::ConstBoundary<GlobalGridPartType> BoundaryGridPartType;

//  typedef typename BoundaryGridPartType::GridViewType BoundaryGridViewType;

  typedef Dune::grid::Part::Local::IndexBased::ConstCoupling<GlobalGridPartType> CouplingGridPartType;

//  typedef typename CouplingGridPartType::GridViewType CouplingGridViewType;

  typedef typename GlobalGridPartType::template Codim<0>::EntityType EntityType;

  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef std::set<size_t> NeighborSetType;

  //! map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map<IndexType, size_t> EntityToSubdomainMapType;

  static const std::string id()
  {
    return "grid.multiscale.default";
  }

  Default(const std::shared_ptr<const GridType> grd,
          const std::shared_ptr<const GlobalGridPartType> globalGrdPrt,
          const size_t sz,
          const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSets,
          const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdMap,
          const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> localGridParts,
          const std::shared_ptr<const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>> boundaryGridParts,
          const std::shared_ptr<const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>>
              couplingGridPartsMaps)
    : grid_(grd)
    , globalGridPart_(globalGrdPrt)
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(false)
  //    , localGridViews_(new std::vector<std::shared_ptr<const LocalGridViewType>>(size_))
  //    , boundaryGridViews_(new std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>())
  //    , couplingGridViewsMaps_(new std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>(
  //          size_, std::map<size_t, std::shared_ptr<const CouplingGridViewType>>()))
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error in " << id() << ":" << std::endl;
    if (localGridParts_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << localGridParts_->size() << ", should be " << size_ << ")!"
          << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'couplingGridPartsMaps' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be "
          << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    //    // create grid views
    //    createGridViews();
  } // Default()

  Default(const std::shared_ptr<const GridType> grd,
          const std::shared_ptr<const GlobalGridPartType> globalGrdPrt,
          const size_t sz,
          const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSets,
          const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdMap,
          const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> localGridParts,
          const std::shared_ptr<const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>> boundaryGridParts,
          const std::shared_ptr<const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>>
              couplingGridPartsMaps,
          const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> oversampledLocalGridParts)
    : grid_(grd)
    , globalGridPart_(globalGrdPrt)
    , size_(sz)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(true)
    , oversampledLocalGridParts_(oversampledLocalGridParts)
  //    , localGridViews_(new std::vector<std::shared_ptr<const LocalGridViewType>>(size_))
  //    , boundaryGridViews_(new std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>())
  //    , couplingGridViewsMaps_(new std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>(
  //          size_, std::map<size_t, std::shared_ptr<const CouplingGridViewType>>()))
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error in " << id() << ":" << std::endl;
    if (localGridParts_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << localGridParts_->size() << ", should be " << size_ << ")!"
          << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'couplingGridPartsMaps' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be "
          << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    //    // create grid views
    //    createGridViews();
  } // Default()

  Default(const ThisType& other) = default;

  Default(ThisType&& source) = default;

  const std::shared_ptr<const GridType>& grid() const
  {
    return grid_;
  }

  const GlobalGridPartType& globalGridPart() const
  {
    return *globalGridPart_;
  }

  //  GlobalGridViewType globalGridView() const
  //  {
  //    return *globalGridView_;
  //  }

  size_t size() const
  {
    return size_;
  }

  bool oversampling() const
  {
    return oversampling_;
  }

  LocalGridPartType localGridPart(const size_t subdomain, const bool ovrsmplng = false) const
  {
    assert(subdomain < size_);
    if (!ovrsmplng) {
      const std::vector<std::shared_ptr<const LocalGridPartType>>& localGridParts = *localGridParts_;
      return *(localGridParts[subdomain]);
    } else {
      if (!oversampling_)
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::XT::Common::color_string_red("ERROR:")
                        << " oversampled local gridpart requested from a grid without oversampling!");
      const std::vector<std::shared_ptr<const LocalGridPartType>>& oversampledLocalGridParts =
          *oversampledLocalGridParts_;
      return *(oversampledLocalGridParts[subdomain]);
    }
  } // ... localGridPart(...)

  //  LocalGridViewType localGridView(const size_t subdomain) const
  //  {
  //    DUNE_THROW(NotImplemented, "The grid views are unsafe at the moment!");
  //    assert(subdomain < size_);
  //    const std::vector<std::shared_ptr<const LocalGridViewType>>& localGridViews = *localGridViews_;
  //    return *(localGridViews[subdomain]);
  //  }

  bool boundary(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>& boundaryGridParts = *boundaryGridParts_;
    return (boundaryGridParts.find(subdomain) != boundaryGridParts.end());
  }

  BoundaryGridPartType boundaryGridPart(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>& boundaryGridParts = *boundaryGridParts_;
    typename std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>::const_iterator result =
        boundaryGridParts.find(subdomain);
    assert(result != boundaryGridParts.end()
           && "Only call boundaryGridPart(subdomain), if boundary(subdomain) is true!");
    return *(result->second);
  }

  //  BoundaryGridViewType boundaryGridView(const size_t subdomain) const
  //  {
  //    DUNE_THROW(NotImplemented, "The grid views are unsafe at the moment!");
  //    assert(subdomain < size_);
  //    const std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>& boundaryGridViews = *boundaryGridViews_;
  //    typename std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>::const_iterator result =
  //        boundaryGridViews.find(subdomain);
  //    assert(result != boundaryGridViews.end()
  //           && "Only call boundaryGridView(subdomain), if boundary(subdomain) is true!");
  //    return *(result->second);
  //  }

  CouplingGridPartType couplingGridPart(const size_t subdomain, const size_t neighbor) const
  {
    assert(subdomain < size_);
    assert(neighbor < size_);
    const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>& couplingGridPartsMaps =
        *couplingGridPartsMaps_;
    const std::map<size_t, std::shared_ptr<const CouplingGridPartType>>& couplingGridPartsMap =
        couplingGridPartsMaps[subdomain];
    const typename std::map<size_t, std::shared_ptr<const CouplingGridPartType>>::const_iterator result =
        couplingGridPartsMap.find(neighbor);
    if (result == couplingGridPartsMap.end()) {
      std::stringstream msg;
      msg << "Error in " << id() << ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain
          << "!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return *(result->second);
  } // const std::shared_ptr< const CouplingGridPartType > couplingGridPart(const size_t subdomain, const size_t
  // neighbor) const

  //  CouplingGridViewType couplingGridView(const size_t subdomain, const size_t neighbor) const
  //  {
  //    DUNE_THROW(NotImplemented, "The grid views are unsafe at the moment!");
  //    assert(subdomain < size_);
  //    assert(neighbor < size_);
  //    const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>& couplingGridViewsMaps =
  //        *couplingGridViewsMaps_;
  //    const std::map<size_t, std::shared_ptr<const CouplingGridViewType>>& couplingGridViewsMap =
  //        couplingGridViewsMaps[subdomain];
  //    const typename std::map<size_t, std::shared_ptr<const CouplingGridViewType>>::const_iterator result =
  //        couplingGridViewsMap.find(neighbor);
  //    if (result == couplingGridViewsMap.end()) {
  //      std::stringstream msg;
  //      msg << "Error in " << id() << ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain
  //          << "!";
  //      DUNE_THROW(Dune::InvalidStateException, msg.str());
  //    }
  //    return *(result->second);
  //  } // const std::shared_ptr< const CouplingGridViewType > couplingGridView(const size_t subdomain, const size_t
  //  // neighbor) const

  const std::shared_ptr<const EntityToSubdomainMapType>& entityToSubdomainMap() const
  {
    return entityToSubdomainMap_;
  }

  const NeighborSetType& neighborsOf(const size_t subdomain) const
  {
    assert(subdomain < size_);
    const std::vector<NeighborSetType>& neighboringSets = *neighboringSetsPtr_;
    return neighboringSets[subdomain];
  } // const NeighborSetType& neighborsOf(const size_t subdomain) const

  size_t subdomainOf(const IndexType& globalIndex) const
  {
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "Error in " << id() << ": missing information for entity " << globalIndex << "in entityToSubdomainMap_!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // size_t getSubdomainOf(const IndexType& globalIndex) const

  size_t subdomainOf(const EntityType& entity) const
  {
    return subdomainOf(globalGridPart_->indexSet().index(entity));
  } // size_t subdomainOf(const EntityType& entity) const

private:
  //  void createGridViews()
  //  {
  //    // create global grid view
  //    globalGridView_ = std::shared_ptr<const GlobalGridViewType>(new
  //    GlobalGridViewType(globalGridPart_->gridView()));
  //    // walk the subdomains
  //    //   * to create the local grid views
  //    //   * to create the boundary grid views
  //    //   * to create the coupling grid views
  //    for (size_t subdomain = 0; subdomain < size_; ++subdomain) {
  //      // for the local grid view
  //      //   * get the local grid part
  //      const std::vector<std::shared_ptr<const LocalGridPartType>>& localGridParts = *localGridParts_;
  //      const std::shared_ptr<const LocalGridPartType>& localGridPartPtr = localGridParts[subdomain];
  //      const LocalGridPartType& localGrdPrt = *localGridPartPtr;
  //      //   * and create the local grid view
  //      std::vector<std::shared_ptr<const LocalGridViewType>>& localGridViews = *localGridViews_;
  //      localGridViews[subdomain] =
  //          std::shared_ptr<const LocalGridViewType>(new LocalGridViewType(localGrdPrt.gridView()));
  //      // for the coupling grid views
  //      //   * get the grid parts map for this subdomain
  //      const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>& couplingGridPartsMaps =
  //          *couplingGridPartsMaps_;
  //      const std::map<size_t, std::shared_ptr<const CouplingGridPartType>>& couplingGridPartsMap =
  //          couplingGridPartsMaps[subdomain];
  //      //   * get the target grid view map for this subdomain
  //      std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>& couplingGridViewsMaps =
  //          *couplingGridViewsMaps_;
  //      std::map<size_t, std::shared_ptr<const CouplingGridViewType>>& couplingGridViewsMap =
  //          couplingGridViewsMaps[subdomain];
  //      //   * walk the neighbors
  //      const std::vector<NeighborSetType>& neighboringSets = *neighboringSetsPtr_;
  //      const NeighborSetType& neighbors = neighboringSets[subdomain];
  //      for (typename NeighborSetType::const_iterator neighborIt = neighbors.begin(); neighborIt != neighbors.end();
  //           ++neighborIt) {
  //        const size_t neighbor = *neighborIt;
  //        // * get the coupling grid part
  //        typename std::map<size_t, std::shared_ptr<const CouplingGridPartType>>::const_iterator
  //        couplingGridPartsMapIt =
  //            couplingGridPartsMap.find(neighbor);
  //        assert(couplingGridPartsMapIt != couplingGridPartsMap.end()
  //               && "Error: missing coupling grid part in given 'couplingGridPartsMaps'!");
  //        const std::shared_ptr<const CouplingGridPartType>& couplingGridPartPtr = couplingGridPartsMapIt->second;
  //        const CouplingGridPartType& couplingGrdPrt = *couplingGridPartPtr;
  //        // * and create the coupling grid view
  //        couplingGridViewsMap.insert(std::pair<size_t, std::shared_ptr<const CouplingGridViewType>>(
  //            neighbor,
  //            std::shared_ptr<const CouplingGridViewType>(new CouplingGridViewType(couplingGrdPrt.gridView()))));
  //      } // walk the neighbors
  //    } // walk the subdomains
  //    // walk those subdomains that have a boundary grid part
  //    //   * to create the boundary grid views
  //    for (typename std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>::const_iterator boundaryGridPartIt =
  //             boundaryGridParts_->begin();
  //         boundaryGridPartIt != boundaryGridParts_->end();
  //         ++boundaryGridPartIt) {
  //      const size_t subdomain = boundaryGridPartIt->first;
  //      const std::shared_ptr<const BoundaryGridPartType> boundaryGrdPrt = boundaryGridPartIt->second;
  //      assert(boundaryGridViews_->find(subdomain) == boundaryGridViews_->end() && "This should not happen!");
  //      boundaryGridViews_->insert(std::pair<size_t, std::shared_ptr<const BoundaryGridViewType>>(
  //          subdomain,
  //          std::shared_ptr<const BoundaryGridViewType>(new BoundaryGridViewType(boundaryGrdPrt->gridView()))));
  //    } // walk those subdomains that have a boundary grid part
  //  } // void createGridViews()

  const std::shared_ptr<const GridType> grid_;
  const std::shared_ptr<const GlobalGridPartType> globalGridPart_;
  const size_t size_;
  const std::shared_ptr<const std::vector<NeighborSetType>> neighboringSetsPtr_;
  const std::shared_ptr<const EntityToSubdomainMapType> entityToSubdomainMap_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> localGridParts_;
  const std::shared_ptr<const std::map<size_t, std::shared_ptr<const BoundaryGridPartType>>> boundaryGridParts_;
  const std::shared_ptr<const std::vector<std::map<size_t, std::shared_ptr<const CouplingGridPartType>>>>
      couplingGridPartsMaps_;
  bool oversampling_;
  const std::shared_ptr<const std::vector<std::shared_ptr<const LocalGridPartType>>> oversampledLocalGridParts_;
  //  std::shared_ptr<std::vector<std::shared_ptr<const LocalGridViewType>>> localGridViews_;
  //  std::shared_ptr<std::map<size_t, std::shared_ptr<const BoundaryGridViewType>>> boundaryGridViews_;
  //  std::shared_ptr<std::vector<std::map<size_t, std::shared_ptr<const CouplingGridViewType>>>>
  //  couplingGridViewsMaps_;
  //  std::shared_ptr<const GlobalGridViewType> globalGridView_;
}; // class Default

#else // HAVE_DUNE_FEM

template <class GridImp>
class Default
{
  static_assert(AlwaysFalse<GridImp>::value, "You are missing dune-fem!");
};

#endif // HAVE_DUNE_FEM

} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_DEFAULT_HH
