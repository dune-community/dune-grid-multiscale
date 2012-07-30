
#ifndef DUNE_GRID_MULTISCALE_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_DEFAULT_HH

// system
#include <vector>
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

// gune-grid
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// dune-geometry
#include <dune/geometry/type.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>

// dune-grid-multiscale
#include <dune/grid/part/leaf.hh>
#include <dune/grid/part/local/indexbased.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

/**
 *  \attention  When generating neighbor information, the boundary id is taken as the neighbors (first) subdomain. This is not unique for some overlapping cases!
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

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

private:
  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef Dune::GeometryType GeometryType;

  // for the indices per subdomain
  typedef std::map< IndexType, IndexType > IndexMapType;

  typedef std::map< GeometryType, IndexMapType > GeometryMapType;

  typedef std::map< unsigned int, Dune::shared_ptr< GeometryMapType > > SubdomainMapType;

  // for entity index -> subdomain map
  typedef std::map< IndexType, std::set< unsigned int > > IndexToSubdomainMapType;

  // for the neighbor information
  typedef std::map< int, int > NeighborInfoContainerType;

  typedef std::map< IndexType, NeighborInfoContainerType > IndexToNeighborInfoMapType;

  typedef std::vector< Dune::shared_ptr< IndexToNeighborInfoMapType > > SubdomainNeighborInfoVectorType;

  typedef std::vector< Dune::shared_ptr< const LocalGridPartType > > LocalGridPartVectorType;

  typedef std::vector< Dune::shared_ptr< const LocalGridViewType > > LocalGridViewVectorType;

  template< int c, int d >
  struct Add
  {
    static void subEntities(ThisType& msGrid, const EntityType& entity, GeometryMapType& geometryMap/*, Dune::Stuff::Common::LogStream& out, const std::string prefix*/)
    {
      // suppress output, since we are not codim 0
      Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
      typedef typename EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
      for (int i = 0; i < entity.template count< c >(); ++i)
      {
        const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
        const GeometryType& geometryType = codimCentityPtr->type();
        const IndexType globalIndex = msGrid.globalGridPart_->indexSet().index(*codimCentityPtr);
        msGrid.addGeometryAndIndex(geometryMap, geometryType, globalIndex, devnull, "");
      }
      Add< c + 1, d >::subEntities(msGrid, entity, geometryMap/*, out, prefix*/);
    }
  };

public:
  Default(GridType& grid)
    : grid_(grid),
      globalGridPart_(Dune::shared_ptr< GlobalGridPartType >(new GlobalGridPartType(grid_))),
      globalGridView_(Dune::shared_ptr< GlobalGridViewType >(new GlobalGridViewType(globalGridPart_->gridView()))),
      finalized_(false),
      size_(0)
  {}

  Dune::shared_ptr< const GlobalGridPartType > globalGridPart() const
  {
    return globalGridPart_;
  }

  unsigned int size() const
  {
    return size_;
  }

  void prepare()
  {
    return;
  }

  void add(const EntityType& entity, const unsigned int subdomain, Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    // preparations
    assert(!finalized_);
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    IndexType globalIndex = globalGridPart_->indexSet().index(entity);
    debug << prefix << id << ".add entity " << globalIndex << " to subdomain " << subdomain << ":" << std::endl;
    // add subdomain to this entity index
    typename IndexToSubdomainMapType::iterator indexIt = indexToSubdomainMap_.find(globalIndex);
    if (indexIt == indexToSubdomainMap_.end()) {
      indexToSubdomainMap_.insert(std::pair< IndexType, std::set< unsigned int > >(globalIndex, std::set< unsigned int >()));
      indexIt = indexToSubdomainMap_.find(globalIndex);
    }
    indexIt->second.insert(subdomain);
    // create geometry map for this subdomain if needed (doing this explicitly instead of just using insert only to increment size)
    if (subdomainMap_.find(subdomain) == subdomainMap_.end()) {
      subdomainMap_.insert(std::pair< unsigned int, Dune::shared_ptr< GeometryMapType > >(
        subdomain,
        Dune::shared_ptr< GeometryMapType >(new GeometryMapType())));
      ++size_;
    }
    // get geometry map for this subdomain
    GeometryMapType& geometryMap = *(subdomainMap_.find(subdomain)->second);
    // add geometry and global index of this codim 0 entity
    GeometryType geometryType = entity.type();
    addGeometryAndIndex(geometryMap, geometryType, globalIndex, debug, prefix);
    // add all remaining codims
    Add< 1, dim >::subEntities(*this, entity, geometryMap/*, debug, prefix*/);
  } // void add(const EntityType& entity, const unsigned int subdomain, Dune::ParameterTree paramTree = Dune::ParameterTree())

  void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    // debug output
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    const unsigned int boundaryId = paramTree.get("boundaryId", 5);
    debug << prefix << "finalizing " << std::flush;
    // test for consecutive numbering of subdomains and create map for the neighbor informatio
    for (unsigned int subdomain = 0; subdomain < size(); ++subdomain) {
      if (subdomainMap_.find(subdomain) == subdomainMap_.end()) {
        std::stringstream msg;
        msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
        DUNE_THROW(Dune::InvalidStateException, msg.str());
      } else {
        subdomainNeighborInfoVector_.push_back(Dune::shared_ptr< IndexToNeighborInfoMapType >(new IndexToNeighborInfoMapType()));
      }
    }
    debug << size() << " subdomains:" << std::endl;
    // walk the global grid part to collect information about the subdomain boundaries
    debug << prefix << "  generating neighbor information:" << std::endl;
    for (typename GlobalGridPartType::template Codim< 0 >::IteratorType entityIt = globalGridPart_->template begin< 0 >();
         entityIt != globalGridPart_->template end< 0 >();
         ++entityIt) {
      // find the subdomains this entity lives in
      const EntityType& entity = *entityIt;
      const IndexType& globalEntityIndex = globalGridPart_->indexSet().index(entity);
      const std::set< unsigned int >& entitySubdomains = getSubdomainsOf(globalEntityIndex);
      // walk the neighbors
      for (typename GlobalGridPartType::IntersectionIteratorType intersectionIt = globalGridPart_->ibegin(entity);
           intersectionIt != globalGridPart_->iend(entity);
           ++intersectionIt) {
        typedef typename GlobalGridPartType::IntersectionIteratorType::Intersection IntersectionType;
        const IntersectionType& intersection = *intersectionIt;
        if (intersection.neighbor()) {
          typedef typename IntersectionType::EntityPointer EntityPtrType;
          const EntityPtrType neighborPtr = intersection.outside();
          const EntityType& neighbor = *neighborPtr;
          const IndexType& globalNeighborIndex = globalGridPart_->indexSet().index(neighbor);
          const std::set< unsigned int >& neighborSubdomains = getSubdomainsOf(globalNeighborIndex);
          // treat each subdomain of the entity seperately
          for (std::set< unsigned int >::const_iterator subdomainIt = entitySubdomains.begin();
               subdomainIt != entitySubdomains.end();
               ++subdomainIt) {
            const unsigned int entitySubdomain = *subdomainIt;
            // if neighbor is not contained in this subdomain
            if (neighborSubdomains.find(entitySubdomain) == neighborSubdomains.end()) {
              // guess neighbor subdomain as the first this one belongs to (not good for big overlapping)
              unsigned int neighborSubdomain = *(neighborSubdomains.begin());
              // get local index of intersection
              const int localIntersectionIndex = intersection.indexInInside();
              // then this entity is at the boundary
              debug << prefix << "    entity " << globalEntityIndex
                    << " lies at the boundary of subdomain " << entitySubdomain << std::endl;
              debug << prefix << "      at intersection " << localIntersectionIndex
                    << " with neighbor " << globalNeighborIndex
                    << ", (probably) of subdomain " << neighborSubdomain << std::endl;
              // get the correct map
              IndexToNeighborInfoMapType& indexToNeighborMap = *(subdomainNeighborInfoVector_[entitySubdomain]);
              // get the entry for this entity
              NeighborInfoContainerType& neighborInfoContainer = indexToNeighborMap[globalEntityIndex];
              // add the local intersection id and the boundary id this intersection shall have
              neighborInfoContainer.insert(std::pair< int, int >(localIntersectionIndex, boundaryId));
            } // if neighbor is not contained in this subdomain
          } // treat each subdomain of the entity seperately
        } // if (intersection.neighbor())
      } // walk the neighbors
    } // walk the global grid part to collect information about the subdomain boundaries

    // create local grid parts and report
    debug << prefix << "  creating local grid parts:" << std::endl;
    for (typename SubdomainMapType::iterator subdomainIterator = subdomainMap_.begin();
         subdomainIterator != subdomainMap_.end();
         ++subdomainIterator) {
      const unsigned int subdomain = subdomainIterator->first;
      // compute number of codim 0 entities
      unsigned int subdomainSize = 0;
      for (typename GeometryMapType::iterator geometryIterator = subdomainIterator->second->begin();
           geometryIterator != subdomainIterator->second->end();
           ++geometryIterator) {
        if (geometryIterator->first.dim() == dim)
          subdomainSize += geometryIterator->second.size();
      } // compute number of codim 0 entities
      debug << prefix << "    subdomain " << subdomain << " of size " << subdomainSize << "... " << std::flush;
      localGridParts_.push_back(Dune::shared_ptr< const LocalGridPartType >(
        new LocalGridPartType(*globalGridPart_, subdomainIterator->second, subdomainNeighborInfoVector_[subdomain])));
      localGridViews_.push_back(Dune::shared_ptr< const LocalGridViewType >(
        new LocalGridViewType(localGridParts_[subdomain]->gridView())));
      debug << "done" << std::endl;
    } // create local grid parts and report

    // finalize
    finalized_ = true;
    return;
  } // void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())

  const Dune::shared_ptr< const LocalGridPartType > localGridPart(const unsigned int subdomain) const
  {
    assert(finalized_);
    assert(subdomain < size());
    return localGridParts_[subdomain];
  }

  void visualize(Dune::ParameterTree paramTree = Dune::ParameterTree()) const
  {
    // preparations
    assert(finalized_);
    const std::string filename = paramTree.get("filename", "msGrid_visualization");

    // vtk writer
    Dune::VTKWriter< GlobalGridViewType > vtkwriter(*globalGridView_);

    // subdomain id
    std::vector< double > subdomainId = generateSubdomainVisualization();
    vtkwriter.addCellData(subdomainId, "subdomainId");

    // boundary id
    std::vector< double > boundaryId = generateBoundaryIdVisualization();
    vtkwriter.addCellData(boundaryId, "boundaryId");

    // codim 0 entity id
    std::vector< double > entityId = generateEntityVisualization();
    vtkwriter.addCellData(entityId, "entityId");

    // write
    vtkwriter.write(filename, Dune::VTK::ascii);
  } // void visualize(Dune::ParameterTree paramTree = Dune::ParameterTree()) const

private:
  template< int, int >
  friend class Add;

  void addGeometryAndIndex(GeometryMapType& geometryMap,
                           const GeometryType& geometryType,
                           const IndexType& globalIndex,
                           Dune::Stuff::Common::LogStream& out,
                           const std::string prefix)
  {
    // make sure the map to this geometry type exists
    typename GeometryMapType::iterator iterator = geometryMap.find(geometryType);
    if (iterator == geometryMap.end())
      geometryMap.insert(std::pair< GeometryType, IndexMapType >(geometryType, IndexMapType()));
    // add global (and local) index [we make use of the fact that a map does not change existing pairs on insert]
    IndexMapType& indexMap = geometryMap.find(geometryType)->second;
    if (indexMap.find(globalIndex) == indexMap.end()) {
      const IndexType localIndex = indexMap.size();
      indexMap.insert(std::pair< IndexType, IndexType >(globalIndex, localIndex));
      out << prefix << "- added " << geometryType << " with global index " << globalIndex << " and local index " << localIndex << std::endl;
    }
    return;
  } // void addGeometryAndIndex

  const std::set< unsigned int >& getSubdomainsOf(const IndexType& globalIndex) const
  {
    const typename IndexToSubdomainMapType::const_iterator result = indexToSubdomainMap_.find(globalIndex);
    if (result == indexToSubdomainMap_.end()) {
      std::stringstream msg;
      msg << "Error in " << id << ": missing information for entity " << globalIndex << "in indexToSubdomainMap_!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // unsigned int getSubdomainOf(const IndexType& globalIndex) const

  std::vector< double > generateSubdomainVisualization() const
  {
    std::vector< double > data(globalGridView_->indexSet().size(0));
    // walk the grid
    for (typename GlobalGridViewType::template Codim< 0 >::Iterator it = globalGridView_->template begin< 0 >();
         it != globalGridView_->template end< 0 >();
         ++it)
    {
      const EntityType& entity = *it;
      const IndexType& index = globalGridView_->indexSet().index(entity);
      const std::set< unsigned int >& subdomains = getSubdomainsOf(index);
      data[index] = *(subdomains.begin());
    } // walk the grid
    return data;
  } // std::vector< double > generateSubdomainVisualization() const

  std::vector< double > generateEntityVisualization() const
  {
    std::vector< double > data(globalGridView_->indexSet().size(0));
    // walk the grid
    for (typename GlobalGridViewType::template Codim< 0 >::Iterator it = globalGridView_->template begin< 0 >();
         it != globalGridView_->template end< 0 >();
         ++it)
    {
      const EntityType& entity = *it;
      const IndexType& index = globalGridView_->indexSet().index(entity);
      data[index] = double(index);
    } // walk the grid
    return data;
  } // std::vector< double > generateEntityVisualization() const

  std::vector< double > generateBoundaryIdVisualization() const
  {
    std::vector< double > data(globalGridView_->indexSet().size(0));
    // walk the grid
    for (typename GlobalGridViewType::template Codim< 0 >::Iterator it = globalGridView_->template begin< 0 >();
         it != globalGridView_->template end< 0 >();
         ++it)
    {
      const EntityType& entity = *it;
      const IndexType& index = globalGridView_->indexSet().index(entity);
      data[index] = 0.0;
      int numberOfBoundarySegments = 0;
      bool isOnBoundary = false;
      for (typename GlobalGridViewType::IntersectionIterator intersectionIt = globalGridView_->ibegin(entity);
           intersectionIt != globalGridView_->iend(entity);
           ++intersectionIt) {
        if (!intersectionIt->neighbor() && intersectionIt->boundary()){
          isOnBoundary = true;
          numberOfBoundarySegments += 1;
          data[index] += double(intersectionIt->boundaryId());
        }
      }
      if (isOnBoundary) {
        data[index] /= double(numberOfBoundarySegments);
      }
    } // walk the grid
    return data;
  } // std::vector< double > generateBoundaryIdVisualization() const

  GridType& grid_;
  Dune::shared_ptr< GlobalGridPartType > globalGridPart_;
  Dune::shared_ptr< GlobalGridViewType > globalGridView_;
  bool finalized_;
  unsigned int size_;
  SubdomainMapType subdomainMap_;
  IndexToSubdomainMapType indexToSubdomainMap_;
  SubdomainNeighborInfoVectorType subdomainNeighborInfoVector_;
  LocalGridPartVectorType localGridParts_;
  LocalGridViewVectorType localGridViews_;
}; // class Default

template< class GridType >
const std::string Default< GridType >::id = "grid.multiscale.default";

//! specialization to stop the recursion
template< class GridType >
template< int c >
struct Default< GridType >::Add< c, c >
{
  static void subEntities(Default< GridType >& msGrid,
                          const typename Default< GridType >::EntityType& entity,
                          typename Default< GridType >::GeometryMapType& geometryMap/*,
                          Dune::Stuff::Common::LogStream& out,
                          const std::string prefix*/)
  {
    // suppress output, since we are not codim 0
    Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
    typedef typename Default< GridType >::EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
    for (int i = 0; i < entity.template count< c >(); ++i)
    {
      const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
      const Default< GridType >::GeometryType& geometryType = codimCentityPtr->type();
      const typename Default< GridType >::IndexType globalIndex = msGrid.globalGridPart_->indexSet().index(*codimCentityPtr);
      msGrid.addGeometryAndIndex(geometryMap, geometryType, globalIndex, devnull, "");
    }
  }
};

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_DEFAULT_HH
