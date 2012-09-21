
#ifndef DUNE_GRID_MULTISCALE_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_DEFAULT_HH

// system
#include <vector>
#include <set>
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>

// gune-grid
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// dune-grid-multiscale
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

  static const std::string id;

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

  Default(const Dune::shared_ptr< const GridType > grid,
          const Dune::shared_ptr< const GlobalGridPartType > globalGridPart,
          const unsigned int size,
          const Dune::shared_ptr< const std::vector< NeighborSetType > > neighboringSets,
          const Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap,
          const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const LocalGridPartType > > > localGridParts,
          const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const BoundaryGridPartType > > > boundaryGridParts,
          const Dune::shared_ptr< const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > > > couplingGridPartsMaps)
    : grid_(grid)
    , globalGridPart_(globalGridPart)
    , size_(size)
    , neighboringSetsPtr_(neighboringSets)
    , entityToSubdomainMap_(entityToSubdomainMap)
    , localGridParts_(localGridParts)
    , boundaryGridParts_(boundaryGridParts)
    , couplingGridPartsMaps_(couplingGridPartsMaps)
    , localGridViews_(new std::vector< Dune::shared_ptr< const LocalGridViewType > >(size_))
    , boundaryGridViews_(new std::vector< Dune::shared_ptr< const BoundaryGridViewType > >(size_))
    , couplingGridViewsMaps_(new std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > >(size_, std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >()))
  {
    // check for correct sizes
    std::stringstream msg;
    bool error = false;
    msg << "Error in " << id << ":" << std::endl;
    if (localGridParts_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << localGridParts_->size() << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (boundaryGridParts_->size() != size_) {
      msg << "  - 'boundaryGridParts' has wrong size (is " << boundaryGridParts_->size() << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (couplingGridPartsMaps_->size() != size_) {
      msg << "  - 'localGridParts' has wrong size (is " << couplingGridPartsMaps_->size() << ", should be " << size_ << ")!" << std::endl;
      error = true;
    }
    if (error)
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    // create grid views
    createGridViews();
  } // Default()

//  Dune::shared_ptr< const GlobalGridPartType > globalGridPart() const
//  {
//    return globalGridPart_;
//  }

//  Dune::shared_ptr< const GlobalGridViewType > globalGridView() const
//  {
//    return globalGridView_;
//  }

//  unsigned int size() const
//  {
//    assert(finalized_);
//    return size_;
//  } // unsigned int size() const

//  Dune::shared_ptr< const LocalGridPartType > localGridPart(const unsigned int subdomain) const
//  {
//    assert(finalized_);
//    assert(subdomain < size());
//    return localGridParts_[subdomain];
//  }

//  Dune::shared_ptr< const LocalGridViewType > localGridView(const unsigned int subdomain) const
//  {
//    assert(finalized_);
//    assert(subdomain < size());
//    return localGridViews_[subdomain];
//  }

//  Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const
//  {
//    assert(finalized_);
//    assert(subdomain < size());
//    assert(neighbor < size());
//    const CouplingGridPartMapType& couplingGridPartMap = couplingGridPartMaps_[subdomain];
//    const typename CouplingGridPartMapType::const_iterator result = couplingGridPartMap.find(neighbor);
//    if (result == couplingGridPartMap.end()) {
//      std::stringstream msg;
//      msg << "Error in " << id << ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!";
//      DUNE_THROW(Dune::InvalidStateException, msg.str());
//    }
//    return result->second;
//  } // Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const

//  Dune::shared_ptr< const CouplingGridViewType > couplingGridView(const unsigned int subdomain, const unsigned int neighbor) const
//  {
//    assert(finalized_);
//    assert(subdomain < size());
//    assert(neighbor < size());
//    const CouplingGridViewMapType& couplingGridViewMap = couplingGridViewMaps_[subdomain];
//    const typename CouplingGridViewMapType::const_iterator result = couplingGridViewMap.find(neighbor);
//    if (result == couplingGridViewMap.end()) {
//      std::stringstream msg;
//      msg << "Error in " << id << ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!";
//      DUNE_THROW(Dune::InvalidStateException, msg.str());
//    }
//    return result->second;
//  } // Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const

//  Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap() const
//  {
//    return entityToSubdomainMap_;
//  }

//  const std::set< unsigned int >& neighborsOf(const unsigned int subdomain) const
//  {
//    assert(finalized_);
//    assert(subdomain < size());
//    return neighboringSubdomainsInfoVector_[subdomain];
//  } //  Dune::shared_ptr< const NeighboringSubdomainsSetType > neighborsOf(const unsigned int subdomain) const

//  unsigned int getSubdomainOf(const IndexType& globalIndex) const
//  {
//    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
//    if (result == entityToSubdomainMap_->end()) {
//      std::stringstream msg;
//      msg << "Error in " << id << ": missing information for entity " << globalIndex << "in entityToSubdomainMap_!";
//      DUNE_THROW(Dune::InvalidStateException, msg.str());
//    }
//    return result->second;
//  } // unsigned int getSubdomainOf(const IndexType& globalIndex) const

//  void visualize(const std::string filename = "msGrid_visualization") const
//  {
//    // preparations
//    assert(finalized_);

//    // vtk writer
//    Dune::VTKWriter< GlobalGridViewType > vtkwriter(*globalGridView_);

//    // subdomain id
//    std::vector< double > subdomainId = generateSubdomainVisualization();
//    vtkwriter.addCellData(subdomainId, "subdomainId");

//    // boundary id
//    std::vector< double > boundaryId = generateBoundaryIdVisualization();
//    vtkwriter.addCellData(boundaryId, "boundaryId");

//    // codim 0 entity id
//    std::vector< double > entityId = generateEntityVisualization();
//    vtkwriter.addCellData(entityId, "entityId");

//    // write
//    vtkwriter.write(filename, Dune::VTK::ascii);
//  } // visualize()

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
      // for the boundry grid view
      //   * get the boundary grid part
      const std::vector< Dune::shared_ptr< const BoundaryGridPartType > >& boundaryGridParts = *boundaryGridParts_;
      const Dune::shared_ptr< const BoundaryGridPartType >& boundaryGridPartPtr = boundaryGridParts[subdomain];
      const BoundaryGridPartType& boundaryGridPart = *boundaryGridPartPtr;
      //   * and create the boundary grid view
      std::vector< Dune::shared_ptr< const BoundaryGridViewType > >& boundaryGridViews = *boundaryGridViews_;
      boundaryGridViews[subdomain] = Dune::shared_ptr< const BoundaryGridViewType >(new BoundaryGridViewType(boundaryGridPart.gridView()));
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
  } // void createGridViews()

//  std::vector< double > generateSubdomainVisualization() const
//  {
//    std::vector< double > data(globalGridView_->indexSet().size(0));
//    // walk the grid
//    for (typename GlobalGridViewType::template Codim< 0 >::Iterator it = globalGridView_->template begin< 0 >();
//         it != globalGridView_->template end< 0 >();
//         ++it)
//    {
//      const EntityType& entity = *it;
//      const IndexType& index = globalGridView_->indexSet().index(entity);
//      const unsigned int subdomain = getSubdomainOf(index);
//      data[index] = subdomain;
//    } // walk the grid
//    return data;
//  } // std::vector< double > generateSubdomainVisualization() const

//  std::vector< double > generateEntityVisualization() const
//  {
//    std::vector< double > data(globalGridView_->indexSet().size(0));
//    // walk the grid
//    for (typename GlobalGridViewType::template Codim< 0 >::Iterator it = globalGridView_->template begin< 0 >();
//         it != globalGridView_->template end< 0 >();
//         ++it)
//    {
//      const EntityType& entity = *it;
//      const IndexType& index = globalGridView_->indexSet().index(entity);
//      data[index] = double(index);
//    } // walk the grid
//    return data;
//  } // std::vector< double > generateEntityVisualization() const

//  std::vector< double > generateBoundaryIdVisualization() const
//  {
//    std::vector< double > data(globalGridView_->indexSet().size(0));
//    // walk the grid
//    for (typename GlobalGridViewType::template Codim< 0 >::Iterator it = globalGridView_->template begin< 0 >();
//         it != globalGridView_->template end< 0 >();
//         ++it)
//    {
//      const EntityType& entity = *it;
//      const IndexType& index = globalGridView_->indexSet().index(entity);
//      data[index] = 0.0;
//      int numberOfBoundarySegments = 0;
//      bool isOnBoundary = false;
//      for (typename GlobalGridViewType::IntersectionIterator intersectionIt = globalGridView_->ibegin(entity);
//           intersectionIt != globalGridView_->iend(entity);
//           ++intersectionIt) {
//        if (!intersectionIt->neighbor() && intersectionIt->boundary()){
//          isOnBoundary = true;
//          numberOfBoundarySegments += 1;
//          data[index] += double(intersectionIt->boundaryId());
//        }
//      }
//      if (isOnBoundary) {
//        data[index] /= double(numberOfBoundarySegments);
//      }
//    } // walk the grid
//    return data;
//  } // std::vector< double > generateBoundaryIdVisualization() const

  const Dune::shared_ptr< const GridType > grid_;
  const Dune::shared_ptr< const GlobalGridPartType > globalGridPart_;
  const unsigned int size_;
  const Dune::shared_ptr< const std::vector< NeighborSetType > > neighboringSetsPtr_;
  const Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap_;
  const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const LocalGridPartType > > > localGridParts_;
  const Dune::shared_ptr< const std::vector< Dune::shared_ptr< const BoundaryGridPartType > > > boundaryGridParts_;
  const Dune::shared_ptr< const std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > > > couplingGridPartsMaps_;
  Dune::shared_ptr< std::vector< Dune::shared_ptr< const LocalGridViewType > > > localGridViews_;
  Dune::shared_ptr< std::vector< Dune::shared_ptr< const BoundaryGridViewType > > > boundaryGridViews_;
  Dune::shared_ptr< std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > > > couplingGridViewsMaps_;
  Dune::shared_ptr< const GlobalGridViewType > globalGridView_;
}; // class Default

template< class GridType >
const std::string Default< GridType >::id = "grid.multiscale.default";

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_DEFAULT_HH
