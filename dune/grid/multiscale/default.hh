// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_DEFAULT_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>
#include <set>
#include <map>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/common/color.hh>

#include <dune/grid/part/leaf.hh>
#include <dune/grid/part/local/indexbased.hh>

namespace Dune {
namespace grid {
namespace Multiscale {

/**
 *  \brief      Default implementation of a Dune::grid::Multiscale::Interface
 *
 *              The basic assumption is, that each entity of the global grid part is added to at most one local grid
 *              part. Overlapping or the like is then handled afterwards by adding additional overlapping regions to
 *              an existing local grid part (not yet implemented).
 *  \attention  Works only for one GeometryType per Codim (I think)! Problem is, that indices are unique per GeometryType, not per Codim.
 *  \todo       Resolve the above Issue (should be easy, compute the size as a sum over all GeometryTypes for a given codim)!
 */
template< class GridImp >
class Default
{
public:
  typedef GridImp GridType;

  typedef Default< GridType > ThisType;

  static const unsigned int dim = GridType::dimension;

  typedef Dune::grid::Part::Leaf::Const< GridType > GlobalGridPartType;

  typedef typename GlobalGridPartType::GridViewType GlobalGridViewType;

  typedef Dune::grid::Part::Local::IndexBased::Const< GlobalGridPartType > LocalGridPartType;

  typedef typename LocalGridPartType::GridViewType LocalGridViewType;

  typedef Dune::grid::Part::Local::IndexBased::ConstBoundary< GlobalGridPartType > BoundaryGridPartType;

  typedef typename BoundaryGridPartType::GridViewType BoundaryGridViewType;

  typedef Dune::grid::Part::Local::IndexBased::ConstCoupling< GlobalGridPartType > CouplingGridPartType;

  typedef typename CouplingGridPartType::GridViewType CouplingGridViewType;

  typedef typename GlobalGridPartType::template Codim< 0 >::EntityType EntityType;

  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef std::set< unsigned int > NeighborSetType;

  //! map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map< IndexType, unsigned int > EntityToSubdomainMapType;

  static const std::string id()
  {
    return "grid.multiscale.default";
  }

  Default(const Dune::shared_ptr< const GridType > grid,
          const Dune::shared_ptr< const GlobalGridPartType > globalGridPart,
          const unsigned int size,
          const Dune::shared_ptr< const std::vector< NeighborSetType > > neighboringSets,
          const Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap,
          const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const LocalGridPartType > > > localGridParts,
          const Dune::shared_ptr< const std::map< unsigned int, Dune::shared_ptr< const BoundaryGridPartType > > > boundaryGridParts,
          const Dune::shared_ptr< const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > > > couplingGridPartsMaps)
    : grid_(grid)
    , globalGridPart_(globalGridPart)
    , size_(size)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdomainMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(false)
    , localGridViews_(new std::vector< Dune::shared_ptr< const LocalGridViewType > >(size_))
    , boundaryGridViews_(new std::map< unsigned int, Dune::shared_ptr< const BoundaryGridViewType > >())
    , couplingGridViewsMaps_(new std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > >(size_, std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >()))
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error in " << id() << ":" << std::endl;
    if (localGridParts_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << localGridParts_->size() << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'couplingGridPartsMaps' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    // create grid views
    createGridViews();
  } // Default()

  Default(const Dune::shared_ptr< const GridType > grid,
          const Dune::shared_ptr< const GlobalGridPartType > globalGridPart,
          const unsigned int size,
          const Dune::shared_ptr< const std::vector< NeighborSetType > > neighboringSets,
          const Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap,
          const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const LocalGridPartType > > > localGridParts,
          const Dune::shared_ptr< const std::map< unsigned int, Dune::shared_ptr< const BoundaryGridPartType > > > boundaryGridParts,
          const Dune::shared_ptr< const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > > > couplingGridPartsMaps,
          const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const LocalGridPartType > > > oversampledLocalGridParts)
    : grid_(grid)
    , globalGridPart_(globalGridPart)
    , size_(size)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdomainMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , oversampling_(true)
    , oversampledLocalGridParts_(oversampledLocalGridParts)
    , localGridViews_(new std::vector< Dune::shared_ptr< const LocalGridViewType > >(size_))
    , boundaryGridViews_(new std::map< unsigned int, Dune::shared_ptr< const BoundaryGridViewType > >())
    , couplingGridViewsMaps_(new std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > >(size_, std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >()))
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error in " << id() << ":" << std::endl;
    if (localGridParts_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << localGridParts_->size() << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'couplingGridPartsMaps' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    // create grid views
    createGridViews();
  } // Default()

  const Dune::shared_ptr< const GlobalGridPartType > globalGridPart() const
  {
    return globalGridPart_;
  }

  const Dune::shared_ptr< const GlobalGridViewType > globalGridView() const
  {
    return globalGridView_;
  }

  unsigned int size() const
  {
    return size_;
  }

  bool oversampling() const
  {
    return oversampling_;
  }

  const Dune::shared_ptr< const LocalGridPartType > localGridPart(const unsigned int subdomain,
                                                                  const bool oversampling = false) const
  {
    assert(subdomain < size_);
    if (!oversampling) {
      const std::vector< Dune::shared_ptr< const LocalGridPartType > >& localGridParts = *localGridParts_;
      return localGridParts[subdomain];
    } else {
      if (!oversampling_)
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                   << " oversampled local gridpart requested from a grid without oversampling!");
      const std::vector< Dune::shared_ptr< const LocalGridPartType > >&
          oversampledLocalGridParts = *oversampledLocalGridParts_;
      return oversampledLocalGridParts[subdomain];
    }
  } // ... localGridPart(...)

  const Dune::shared_ptr< const LocalGridViewType > localGridView(const unsigned int subdomain) const
  {
    assert(subdomain < size_);
    const std::vector< Dune::shared_ptr< const LocalGridViewType > >& localGridViews = *localGridViews_;
    return localGridViews[subdomain];
  }

  bool boundary(const unsigned int subdomain) const {
    assert(subdomain < size_);
    const std::map< unsigned int, Dune::shared_ptr< const BoundaryGridPartType > >& boundaryGridParts = *boundaryGridParts_;
    return (boundaryGridParts.find(subdomain) != boundaryGridParts.end());
  }

  const Dune::shared_ptr< const BoundaryGridPartType > boundaryGridPart(const unsigned int subdomain) const
  {
    assert(subdomain < size_);
    const std::map< unsigned int, Dune::shared_ptr< const BoundaryGridPartType > >& boundaryGridParts = *boundaryGridParts_;
    typename std::map< unsigned int, Dune::shared_ptr< const BoundaryGridPartType > >::const_iterator result = boundaryGridParts.find(subdomain);
    assert(result != boundaryGridParts.end() && "Only call boundaryGridPart(subdomain), if boundary(subdomain) is true!");
    return result->second;
  }

  const Dune::shared_ptr< const BoundaryGridViewType > boundaryGridView(const unsigned int subdomain) const
  {
    assert(subdomain < size_);
    const std::map< unsigned int, Dune::shared_ptr< const BoundaryGridViewType > >& boundaryGridViews = *boundaryGridViews_;
    typename std::map< unsigned int, Dune::shared_ptr< const BoundaryGridViewType > >::const_iterator result = boundaryGridViews.find(subdomain);
    assert(result != boundaryGridViews.end() && "Only call boundaryGridView(subdomain), if boundary(subdomain) is true!");
    return result->second;
  }

  const Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const
  {
    assert(subdomain < size_);
    assert(neighbor < size_);
    const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > >& couplingGridPartsMaps = *couplingGridPartsMaps_;
    const std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > >& couplingGridPartsMap = couplingGridPartsMaps[subdomain];
    const typename std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > >::const_iterator result = couplingGridPartsMap.find(neighbor);
    if (result == couplingGridPartsMap.end()) {
      std::stringstream msg;
      msg << "Error in " << id()<< ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // const Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const

  const Dune::shared_ptr< const CouplingGridViewType > couplingGridView(const unsigned int subdomain, const unsigned int neighbor) const
  {
    assert(subdomain < size_);
    assert(neighbor < size_);
    const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > >& couplingGridViewsMaps = *couplingGridViewsMaps_;
    const std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >& couplingGridViewsMap = couplingGridViewsMaps[subdomain];
    const typename std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >::const_iterator result = couplingGridViewsMap.find(neighbor);
    if (result == couplingGridViewsMap.end()) {
      std::stringstream msg;
      msg << "Error in " << id()<< ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // const Dune::shared_ptr< const CouplingGridViewType > couplingGridView(const unsigned int subdomain, const unsigned int neighbor) const

  const Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap() const
  {
    return entityToSubdomainMap_;
  }

  const NeighborSetType& neighborsOf(const unsigned int subdomain) const
  {
    assert(subdomain < size_);
    const std::vector< NeighborSetType >& neighboringSets = *neighboringSetsPtr_;
    return neighboringSets[subdomain];
  } // const NeighborSetType& neighborsOf(const unsigned int subdomain) const

  unsigned int subdomainOf(const IndexType& globalIndex) const
  {
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "Error in " << id()<< ": missing information for entity " << globalIndex << "in entityToSubdomainMap_!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // unsigned int getSubdomainOf(const IndexType& globalIndex) const

  unsigned int subdomainOf(const EntityType& entity) const
  {
    return subdomainOf(globalGridPart_->indexSet().index(entity));
  } // unsigned int subdomainOf(const EntityType& entity) const

private:
  void createGridViews()
  {
    // create global grid view
    globalGridView_ = Dune::shared_ptr< const GlobalGridViewType >(new GlobalGridViewType(globalGridPart_->gridView()));
    // walk the subdomains
    //   * to create the local grid views
    //   * to create the boundary grid views
    //   * to create the coupling grid views
    for (unsigned int subdomain = 0; subdomain < size_; ++subdomain) {
      // for the local grid view
      //   * get the local grid part
      const std::vector< Dune::shared_ptr< const LocalGridPartType > >& localGridParts = *localGridParts_;
      const Dune::shared_ptr< const LocalGridPartType >& localGridPartPtr = localGridParts[subdomain];
      const LocalGridPartType& localGridPart = *localGridPartPtr;
      //   * and create the local grid view
      std::vector< Dune::shared_ptr< const LocalGridViewType > >& localGridViews = *localGridViews_;
      localGridViews[subdomain] = Dune::shared_ptr< const LocalGridViewType >(new LocalGridViewType(localGridPart.gridView()));
      // for the coupling grid views
      //   * get the grid parts map for this subdomain
      const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > >& couplingGridPartsMaps = *couplingGridPartsMaps_;
      const std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > >& couplingGridPartsMap = couplingGridPartsMaps[subdomain];
      //   * get the target grid view map for this subdomain
      std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > >& couplingGridViewsMaps = *couplingGridViewsMaps_;
      std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >& couplingGridViewsMap = couplingGridViewsMaps[subdomain];
      //   * walk the neighbors
      const std::vector< NeighborSetType >& neighboringSets = *neighboringSetsPtr_;
      const NeighborSetType& neighbors = neighboringSets[subdomain];
      for (typename NeighborSetType::const_iterator neighborIt = neighbors.begin();
           neighborIt != neighbors.end();
           ++neighborIt) {
        const unsigned int neighbor = *neighborIt;
        // * get the coupling grid part
        typename std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > >::const_iterator
            couplingGridPartsMapIt = couplingGridPartsMap.find(neighbor);
        assert(couplingGridPartsMapIt != couplingGridPartsMap.end() && "Error: missing coupling grid part in given 'couplingGridPartsMaps'!");
        const Dune::shared_ptr< const CouplingGridPartType >& couplingGridPartPtr = couplingGridPartsMapIt->second;
        const CouplingGridPartType& couplingGridPart = *couplingGridPartPtr;
        // * and create the coupling grid view
        couplingGridViewsMap.insert(std::pair< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >(
                                      neighbor,
                                      Dune::shared_ptr< const CouplingGridViewType >(new CouplingGridViewType(couplingGridPart.gridView()))));
      } // walk the neighbors
    } // walk the subdomains
    // walk those subdomains that have a boundary grid part
    //   * to create the boundary grid views
    for (typename std::map< unsigned int, Dune::shared_ptr< const BoundaryGridPartType > >::const_iterator boundaryGridPartIt = boundaryGridParts_->begin();
         boundaryGridPartIt != boundaryGridParts_->end();
         ++boundaryGridPartIt) {
      const unsigned int subdomain = boundaryGridPartIt->first;
      const Dune::shared_ptr< const BoundaryGridPartType > boundaryGridPart = boundaryGridPartIt->second;
      assert(boundaryGridViews_->find(subdomain) == boundaryGridViews_->end() && "This should not happen!");
      boundaryGridViews_->insert(std::pair< unsigned int, Dune::shared_ptr< const BoundaryGridViewType > >(
                                   subdomain,
                                   Dune::shared_ptr< const BoundaryGridViewType >(new BoundaryGridViewType(boundaryGridPart->gridView()))));
    } // walk those subdomains that have a boundary grid part
  } // void createGridViews()

  const Dune::shared_ptr< const GridType > grid_;
  const Dune::shared_ptr< const GlobalGridPartType > globalGridPart_;
  const unsigned int size_;
  const Dune::shared_ptr< const std::vector< NeighborSetType > > neighboringSetsPtr_;
  const Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap_;
  const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const LocalGridPartType > > > localGridParts_;
  const Dune::shared_ptr< const std::map< unsigned int, Dune::shared_ptr< const BoundaryGridPartType > > > boundaryGridParts_;
  const Dune::shared_ptr< const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > > > couplingGridPartsMaps_;
  bool oversampling_;
  const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const LocalGridPartType > > > oversampledLocalGridParts_;
  Dune::shared_ptr< std::vector< Dune::shared_ptr< const LocalGridViewType > > > localGridViews_;
  Dune::shared_ptr< std::map< unsigned int, Dune::shared_ptr< const BoundaryGridViewType > > > boundaryGridViews_;
  Dune::shared_ptr< std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > > > couplingGridViewsMaps_;
  Dune::shared_ptr< const GlobalGridViewType > globalGridView_;
}; // class Default

} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_DEFAULT_HH
