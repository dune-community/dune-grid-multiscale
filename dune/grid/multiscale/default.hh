
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
 *  \brief      Default implementation of a Dune::grid::Multiscale::Interface
 *
 *              The basic assumption is, that each entity of the global grid part is added to at most one local grid
 *              part. Overlapping or the like is then handled afterwards by addind additional overlapping regions to
 *              an existing local grid part.
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

  typedef Dune::grid::Part::Local::IndexBased::ConstCoupling< GlobalGridPartType > CouplingGridPartType;

  typedef typename CouplingGridPartType::GridViewType CouplingGridViewType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

private:
  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef Dune::GeometryType GeometryType;

  // for the indices per subdomain
  typedef std::map< IndexType, IndexType > IndexMapType;

  typedef std::map< GeometryType, IndexMapType > GeometryMapType;

  typedef std::map< unsigned int, Dune::shared_ptr< GeometryMapType > > SubdomainMapType;

  // for the local grid parts and views
  typedef std::vector< Dune::shared_ptr< const LocalGridPartType > > LocalGridPartVectorType;

  typedef std::vector< Dune::shared_ptr< const LocalGridViewType > > LocalGridViewVectorType;

  // for the coupling grid parts and views
  typedef std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > CouplingGridPartMapType;

  typedef std::map< unsigned int, Dune::shared_ptr< const CouplingGridViewType > > CouplingGridViewMapType;

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
  //! map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map< IndexType, unsigned int > EntityToSubdomainMapType;

  //! set type which contains all neighbors of a aubdomain
  typedef std::set< unsigned int > NeighboringSubdomainsSetType;

private:
  // for the coupling neighbor information
  typedef std::vector< Dune::shared_ptr< NeighboringSubdomainsSetType > > NeighboringSubdomainsInfoVectorType;

public:
  Default(GridType& grid)
    : grid_(grid),
      globalGridPart_(Dune::shared_ptr< GlobalGridPartType >(new GlobalGridPartType(grid_))),
      globalGridView_(Dune::shared_ptr< GlobalGridViewType >(new GlobalGridViewType(globalGridPart_->gridView()))),
      finalized_(false),
      prepared_(false),
      size_(0)
  {}

  Dune::shared_ptr< const GlobalGridPartType > globalGridPart() const
  {
    return globalGridPart_;
  }

  Dune::shared_ptr< const GlobalGridViewType > globalGridView() const
  {
    return globalGridView_;
  }

  unsigned int size() const
  {
    assert(finalized_);
    return size_;
  } // unsigned int size() const

  void prepare()
  {
    assert(!finalized_);
    if (!prepared_) {
      entityToSubdomainMap_ = Dune::shared_ptr< EntityToSubdomainMapType >(new EntityToSubdomainMapType());
      prepared_ = true;
    } // if (!prepared_)
  } // void prepare()

  void add(const EntityType& entity, const unsigned int subdomain, Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    // preparations
    assert(prepared_);
    assert(!finalized_);
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    IndexType globalIndex = globalGridPart_->indexSet().index(entity);
    debug << prefix << id << ".add entity " << globalIndex << " to subdomain " << subdomain << ":" << std::endl;
    // add subdomain to this entity index
    typename EntityToSubdomainMapType::iterator indexIt = entityToSubdomainMap_->find(globalIndex);
    if (indexIt == entityToSubdomainMap_->end()) {
      entityToSubdomainMap_->insert(std::pair< IndexType, unsigned int >(globalIndex, subdomain));
    } else {
      if (indexIt->second != subdomain) {
        std::stringstream msg;
        msg << "Error in " << id << ": can not add entity to more than one subdomain!";
        DUNE_THROW(Dune::InvalidStateException, msg.str());
      }
    }
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
    const GeometryType& geometryType = entity.type();
    addGeometryAndIndex(geometryMap, geometryType, globalIndex, debug, prefix);
    // add all remaining codims
    Add< 1, dim >::subEntities(*this, entity, geometryMap/*, debug, prefix*/);
  } // void add(const EntityType& entity, const unsigned int subdomain, Dune::ParameterTree paramTree = Dune::ParameterTree())

  void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    assert(prepared_);
    if (!finalized_) {
      // debug output
      const std::string prefix = paramTree.get("prefix", "");
      Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
      const unsigned int boundaryId = paramTree.get("boundaryId", 7);
      debug << prefix << "finalizing " << std::flush;

      // prepare data structures
      // entity indices for the coupling grid parts
      typedef std::vector< SubdomainMapType > CouplingMapsVectorType;
      CouplingMapsVectorType couplingMapsVector;
      // entity indices and intersection indices for the coupling grid parts
      typedef std::set< int > IntersectionSetType;
      typedef std::map< IndexType, IntersectionSetType > IntersectionSetMapType;
      typedef std::map< unsigned int, Dune::shared_ptr< IntersectionSetMapType > > CouplingIntersectionMapType;
      std::vector< CouplingIntersectionMapType > couplingIntersectionMaps;
      // entity indices and intersection indices for the boundary information
      typedef std::map< int, int > NeighborInfoContainerType;
      typedef std::map< IndexType, NeighborInfoContainerType > IndexToNeighborInfoMapType;
      std::vector< Dune::shared_ptr< IndexToNeighborInfoMapType > > subdomainNeighborInfoVector;

      // test for consecutive numbering of subdomains and init data structures
      for (unsigned int subdomain = 0; subdomain < size_; ++subdomain) {
        if (subdomainMap_.find(subdomain) == subdomainMap_.end()) {
          std::stringstream msg;
          msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
          DUNE_THROW(Dune::InvalidStateException, msg.str());
        } else {
          subdomainNeighborInfoVector.push_back(Dune::shared_ptr< IndexToNeighborInfoMapType >(new IndexToNeighborInfoMapType()));
          couplingMapsVector.push_back(SubdomainMapType());
          couplingIntersectionMaps.push_back(CouplingIntersectionMapType());
          neighboringSubdomainsInfoVector_.push_back(Dune::shared_ptr< std::set< unsigned int > >(new std::set< unsigned int >()));
        }
      } // test for consecutive numbering of subdomains and create map for the neighbor information
      debug << size_ << " subdomains:" << std::endl;

      // compute number of codim 0 entities per subdomain
      std::map< unsigned int, unsigned int > localSizes;
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
        localSizes.insert(std::pair< unsigned int, unsigned int >(subdomain, subdomainSize));
      } // compute number of codim 0 entities per subdomain

      // walk the global grid part
      debug << prefix << "  generating neighbor information:" << std::endl;
      for (typename GlobalGridPartType::template Codim< 0 >::IteratorType entityIt = globalGridPart_->template begin< 0 >();
           entityIt != globalGridPart_->template end< 0 >();
           ++entityIt) {
        // find the subdomains this entity lives in
        const EntityType& entity = *entityIt;
        const IndexType& globalEntityIndex = globalGridPart_->indexSet().index(entity);
        const unsigned int entitySubdomain = getSubdomainOf(globalEntityIndex);
        // get the correct map for the boundary information
        IndexToNeighborInfoMapType& indexToNeighborMap = *(subdomainNeighborInfoVector[entitySubdomain]);
        // get the entry for this entity
        NeighborInfoContainerType& neighborInfoContainer = indexToNeighborMap[globalEntityIndex];
        // get the correct set for the neighbor information
        NeighboringSubdomainsSetType& neighborsOfSubdomain = *(neighboringSubdomainsInfoVector_[entitySubdomain]);
        bool isConnected = false;
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
            const unsigned int neighborSubdomain = getSubdomainOf(globalNeighborIndex);
            // if neighbor is not contained in this subdomain
            if (neighborSubdomain != entitySubdomain) {
              // get local index of intersection
              const int localIntersectionIndex = intersection.indexInInside();
              // then this entity is at the boundary
              debug << prefix << "    entity " << globalEntityIndex
                    << " lies at the boundary of subdomain " << entitySubdomain << std::endl;
              debug << prefix << "      at intersection " << localIntersectionIndex
                    << " with neighbor " << globalNeighborIndex
                    << " of subdomain " << neighborSubdomain << std::endl;
              // add the local intersection id and the boundary id this intersection shall have
              neighborInfoContainer.insert(std::pair< int, int >(localIntersectionIndex, boundaryId));
              // add the neighboring subdomain
              neighborsOfSubdomain.insert(neighborSubdomain);
              // == start coupling entity information
              // get coupling map for this subdomain
              SubdomainMapType& couplingMap = couplingMapsVector[entitySubdomain];
              typename SubdomainMapType::iterator result = couplingMap.find(neighborSubdomain);
              // add map for this subdomain / neighbor combination if necessary
              if (result == couplingMap.end()) {
                couplingMap.insert(std::pair< unsigned int, Dune::shared_ptr< GeometryMapType > >(
                  neighborSubdomain, Dune::shared_ptr< GeometryMapType >(new GeometryMapType())));
                result = couplingMap.find(neighborSubdomain);
              }
              // get the geometry map for this subdomain/ neighbor combination
              GeometryMapType& couplingGeometryMap = *(result->second);
              const GeometryType& geometryType = entity.type();
              // add the index of the entity
              addGeometryAndIndex(couplingGeometryMap, geometryType, globalEntityIndex);
              // add all remaining codims
              Add< 1, dim >::subEntities(*this, entity, couplingGeometryMap);
              // == stop coupling entity information
              // == start coupling intersection information
              // get map for this subdomain
              CouplingIntersectionMapType& couplingIntersectionMap = couplingIntersectionMaps[entitySubdomain];
              typename CouplingIntersectionMapType::iterator couplingIntersectionMapResult = couplingIntersectionMap.find(neighborSubdomain);
              // add map for this subdomain / neighbor combination if necessary
              if (couplingIntersectionMapResult == couplingIntersectionMap.end()) {
                couplingIntersectionMap.insert(std::pair< unsigned int, Dune::shared_ptr< IntersectionSetMapType > >(
                  neighborSubdomain, Dune::shared_ptr< IntersectionSetMapType >(new IntersectionSetMapType())));
                couplingIntersectionMapResult = couplingIntersectionMap.find(neighborSubdomain);
              }
              // get the map for this subdomain / neighbor combination
              IntersectionSetMapType& intersectionSetMap = *(couplingIntersectionMapResult->second);
              // add entry for this entity if needed
              typename IntersectionSetMapType::iterator intersectionSetMapResult = intersectionSetMap.find(globalEntityIndex);
              if (intersectionSetMapResult == intersectionSetMap.end()) {
                intersectionSetMap.insert(std::pair< IndexType, IntersectionSetType >(globalEntityIndex, IntersectionSetType()));
                intersectionSetMapResult = intersectionSetMap.find(globalEntityIndex);
              }
              // get intersection index map for this entity
              IntersectionSetType& intersectionSet = intersectionSetMapResult->second;
              // add local intersection index
              intersectionSet.insert(localIntersectionIndex);
              // == stop coupling intersection information
            } else { // if neighbor is contained in this subdomain
              isConnected = true;
            } // if neighbor is contained in this subdomain
          } // if (intersection.neighbor())
        } // walk the neighbors
        // ensure that this entity is connected to the other entities of this subdomain
        if (localSizes[entitySubdomain] != 1)
          assert(isConnected);
      } // walk the global grid part

      // create local grid parts and report
      debug << prefix << "  creating local grid parts:" << std::endl;
      for (typename SubdomainMapType::iterator subdomainIterator = subdomainMap_.begin();
           subdomainIterator != subdomainMap_.end();
           ++subdomainIterator) {
        const unsigned int subdomain = subdomainIterator->first;
        unsigned int subdomainSize  = localSizes[subdomain];
        debug << prefix << "    subdomain " << subdomain << " of size " << subdomainSize << "... " << std::flush;
        localGridParts_.push_back(Dune::shared_ptr< const LocalGridPartType >(
          new LocalGridPartType(*globalGridPart_, subdomainIterator->second, subdomainNeighborInfoVector[subdomain])));
        localGridViews_.push_back(Dune::shared_ptr< const LocalGridViewType >(
          new LocalGridViewType(localGridParts_[subdomain]->gridView())));
        debug << "done" << std::endl;
      } // create local grid parts and report

      // create coupling grid parts and report
      debug << prefix << "  creating coupling grid parts:" << std::endl;
      for (unsigned int subdomain = 0; subdomain < couplingMapsVector.size(); ++subdomain) {
        const SubdomainMapType& couplingMap = couplingMapsVector[subdomain];
        const CouplingIntersectionMapType& couplingIntersectionMap = couplingIntersectionMaps[subdomain];
        couplingGridPartMaps_.push_back(CouplingGridPartMapType());
        couplingGridViewMaps_.push_back(CouplingGridViewMapType());
        CouplingGridPartMapType& couplingGridPartMap = couplingGridPartMaps_[subdomain];
        CouplingGridViewMapType& couplingGridViewMap = couplingGridViewMaps_[subdomain];
        // loop over all neighbors
        for (typename SubdomainMapType::const_iterator neighborIt = couplingMap.begin();
             neighborIt != couplingMap.end();
             ++neighborIt) {
          const unsigned int neighbor = neighborIt->first;
          debug << prefix << "    subdomain " << subdomain << ", neighbor " << neighbor <<  "... " << std::flush;
          // get entity indices map for this neighbor
          const typename SubdomainMapType::const_iterator indicesMap = couplingMap.find(neighbor);
          assert(indicesMap != couplingMap.end());
          // get intersection map for this neighbor
          const typename CouplingIntersectionMapType::const_iterator intersectionIndicesMap = couplingIntersectionMap.find(neighbor);
          assert(intersectionIndicesMap != couplingIntersectionMap.end());
          const Dune::shared_ptr< const CouplingGridPartType > couplingGridPartPtr(new CouplingGridPartType(*globalGridPart_, indicesMap->second, intersectionIndicesMap->second));
          couplingGridPartMap.insert(std::pair< unsigned int, Dune::shared_ptr< const CouplingGridPartType > >(
            neighbor, couplingGridPartPtr));
          couplingGridViewMap.insert(std::pair< unsigned int, Dune::shared_ptr< const CouplingGridViewType > >(
            neighbor, Dune::shared_ptr< const CouplingGridViewType >(new CouplingGridViewType(couplingGridPartPtr->gridView()))));
          debug << "done" << std::endl;
        } // loop over all neighbors
      } // create coupling grid parts and report

      // finalize
      finalized_ = true;
    } // if (!finalized_)
  } // void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())

  Dune::shared_ptr< const LocalGridPartType > localGridPart(const unsigned int subdomain) const
  {
    assert(finalized_);
    assert(subdomain < size());
    return localGridParts_[subdomain];
  }

  Dune::shared_ptr< const LocalGridViewType > localGridView(const unsigned int subdomain) const
  {
    assert(finalized_);
    assert(subdomain < size());
    return localGridViews_[subdomain];
  }

  Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const
  {
    assert(finalized_);
    assert(subdomain < size());
    assert(neighbor < size());
    const CouplingGridPartMapType& couplingGridPartMap = couplingGridPartMaps_[subdomain];
    const typename CouplingGridPartMapType::const_iterator result = couplingGridPartMap.find(neighbor);
    if (result == couplingGridPartMap.end()) {
      std::stringstream msg;
      msg << "Error in " << id << ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const

  Dune::shared_ptr< const CouplingGridViewType > couplingGridView(const unsigned int subdomain, const unsigned int neighbor) const
  {
    assert(finalized_);
    assert(subdomain < size());
    assert(neighbor < size());
    const CouplingGridViewMapType& couplingGridViewMap = couplingGridViewMaps_[subdomain];
    const typename CouplingGridViewMapType::const_iterator result = couplingGridViewMap.find(neighbor);
    if (result == couplingGridViewMap.end()) {
      std::stringstream msg;
      msg << "Error in " << id << ": subdomain " << neighbor << " is not a neighbor of subdomain " << subdomain << "!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // Dune::shared_ptr< const CouplingGridPartType > couplingGridPart(const unsigned int subdomain, const unsigned int neighbor) const

  Dune::shared_ptr< const EntityToSubdomainMapType > entityToSubdomainMap() const
  {
    return entityToSubdomainMap_;
  }

  Dune::shared_ptr< const NeighboringSubdomainsSetType > neighborsOf(const unsigned int subdomain) const
  {
    assert(finalized_);
    assert(subdomain < size());
    return neighboringSubdomainsInfoVector_[subdomain];
  } //  Dune::shared_ptr< const NeighboringSubdomainsSetType > neighborsOf(const unsigned int subdomain) const

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
                           Dune::Stuff::Common::LogStream& out = Dune::Stuff::Common::Logger().devnull(),
                           const std::string prefix = "")
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

  unsigned int getSubdomainOf(const IndexType& globalIndex) const
  {
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "Error in " << id << ": missing information for entity " << globalIndex << "in entityToSubdomainMap_!";
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
      const unsigned int subdomain = getSubdomainOf(index);
      data[index] = subdomain;
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
  bool prepared_;
  unsigned int size_;
  SubdomainMapType subdomainMap_;
  Dune::shared_ptr< EntityToSubdomainMapType > entityToSubdomainMap_;
  NeighboringSubdomainsInfoVectorType neighboringSubdomainsInfoVector_;
  LocalGridPartVectorType localGridParts_;
  LocalGridViewVectorType localGridViews_;
  std::vector< CouplingGridPartMapType > couplingGridPartMaps_;
  std::vector< CouplingGridViewMapType > couplingGridViewMaps_;
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
