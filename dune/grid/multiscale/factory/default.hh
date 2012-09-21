
#ifndef DUNE_GRID_MULTISCALE_FACTORY_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_FACTORY_DEFAULT_HH

// system
#include <vector>
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune-geometry
#include <dune/geometry/type.hh>

// dune-grid-multiscale
#include <dune/grid/part/leaf.hh>
#include <dune/grid/part/local/indexbased.hh>
#include <dune/grid/multiscale/default.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace Factory {

template< class GridImp >
class Default
{
public:
  typedef GridImp GridType;

  typedef Default< GridType > ThisType;

  static const std::string id;

  static const unsigned int dim = GridType::dimension;

  typedef Dune::grid::Multiscale::Default< GridType > MsGridType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

private:
  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

  typedef typename MsGridType::LocalGridPartType LocalGridPartType;

  typedef typename MsGridType::BoundaryGridPartType BoundaryGridPartType;

  typedef typename MsGridType::CouplingGridPartType CouplingGridPartType;

  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef Dune::GeometryType GeometryType;

  // i.e. maps a local to a globl index
  typedef std::map< IndexType, IndexType > IndexMapType;

  // i.e. maps a GeometryType to a map of local and global indices
  typedef std::map< GeometryType, IndexMapType > GeometryMapType;

  // i.e. contains a GeometryMap for each subdomain
  typedef std::map< unsigned int, Dune::shared_ptr< GeometryMapType > > SubdomainMapType;

  // map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map< IndexType, unsigned int > EntityToSubdomainMapType;

  typedef Dune::FieldVector< unsigned int, dim > CodimSizesType;

  // for the neighbor information between the subdomains
  typedef std::set< unsigned int > NeighboringSubdomainsSetType;

  template< int c, int d >
  struct Add
  {
    static void subEntities(ThisType& factory,
                            const EntityType& entity,
                            GeometryMapType& geometryMap,
                            CodimSizesType& localCodimSizes)
    {
      // suppress output, since we are not codim 0
      Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
      // loop over all codim c subentities of the entity
      typedef typename EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
      for (int i = 0; i < entity.template count< c >(); ++i) {
        const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
        const GeometryType& geometryType = codimCentityPtr->type();
        const IndexType globalIndex = factory.globalGridPart_->indexSet().index(*codimCentityPtr);
        factory.addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex, "", devnull);
      } // loop over all codim c subentities of the entity
      // add all codim c + 1 subentities
      Add< c + 1, d >::subEntities(factory, entity, geometryMap, localCodimSizes);
    } // static void subEntities()
  }; // struct Add

public:
  Default(const GridType& grid, const int boundaryId = 7)
    : grid_(Dune::stackobject_to_shared_ptr(grid))
    , boundaryId_(boundaryId)
    , prepared_(false)
    , finalized_(false)
    , size_(0)
  {}

  Default(const Dune::shared_ptr< const GridType > grid, const int boundaryId = 7)
    : grid_(grid)
    , boundaryId_(boundaryId)
    , prepared_(false)
    , finalized_(false)
    , size_(0)
  {}

  void prepare()
  {
    if (!prepared_) {
      globalGridPart_ = Dune::shared_ptr< const GlobalGridPartType >(new GlobalGridPartType(*grid_));
      entityToSubdomainMap_ = Dune::shared_ptr< EntityToSubdomainMapType >(new EntityToSubdomainMapType());
      prepared_ = true;
    } // if (!prepared_)
  } // void prepare()

  const Dune::shared_ptr< const GlobalGridPartType > globalGridPart() const
  {
    assert(prepared_ && "Please call prepare() before calling globalGridPart()!");
    return globalGridPart_;
  }

  void add(const EntityType& entity,
           const unsigned int subdomain,
           const std::string prefix = "",
           std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    // prepare
    assert(prepared_ && "Please call prepare() before calling add()!");
    assert(!finalized_ && "Do not call add() after calling finalized()!");
    const IndexType globalIndex = globalGridPart_->indexSet().index(entity);
    out << prefix << "adding entity " << globalIndex << " to subdomain " << subdomain << ":" << std::endl;

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
    } // add subdomain to this entity index

    // create geometry map for this subdomain if needed (doing this explicitly (instead of just using insert()) only to increment size)
    if (subdomainToEntityMap_.find(subdomain) == subdomainToEntityMap_.end()) {
      subdomainToEntityMap_.insert(std::pair< unsigned int, Dune::shared_ptr< GeometryMapType > >(
          subdomain,
          Dune::shared_ptr< GeometryMapType >(new GeometryMapType())));
      ++size_;
    } // create geometry map for this subdomain if needed
    // create local codim sizes for this subdomain
    localCodimSizes_.insert(std::pair< unsigned int, CodimSizesType >(
        subdomain,
        CodimSizesType(0)));

    // add this entity and all subentities to the geometry map (geometry map has to exist, see above)
    assert(subdomainToEntityMap_.find(subdomain) != subdomainToEntityMap_.end() && "This should not happen!");
    GeometryMapType& geometryMap = *(subdomainToEntityMap_.find(subdomain)->second);
    // add geometry and global index of this codim 0 entity
    const GeometryType& geometryType = entity.type();
    CodimSizesType& localCodimSizes = localCodimSizes_.find(subdomain)->second;
    addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex, prefix, out);
    // add all remaining codims
    Add< 1, dim >::subEntities(*this, entity, geometryMap, localCodimSizes);
  } // void add()

  void finalize(const std::string prefix = "", std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    assert(prepared_ && "Please call prepare() and add() before calling finalize()!");
    if (!finalized_) {
      // prepare
      out << prefix << "finalizing " << std::flush;
      // init data structures
      // for the subdomains inner boundaries
      //   * to map the local intersection index to the desired fake boundary id
      typedef std::map< int, int > IntersectionToBoundaryIdMapType;
      //   * to map the global entity index to one of those maps
      typedef std::map< IndexType, IntersectionToBoundaryIdMapType > EntityToIntersectionInfoMapType;
      //   * to hold one of those maps for each subdomain
      std::vector< Dune::shared_ptr< EntityToIntersectionInfoMapType > > subdomainInnerBoundaryInfos(size_);
      // for the coupling grid parts
      //   * vector to hold a (neighboring subdomain -> entity) map for each subdomain
      typename std::vector< SubdomainMapType > couplingMaps(size_, SubdomainMapType());
      //   * vector to hold a map of coupling sizes
      std::vector< std::map< unsigned int, CodimSizesType > > couplingCodimSizeMaps(size_, std::map< unsigned int, CodimSizesType >());
      //   * set of local coupling intersections
      typedef std::set< int > IntersectionInfoSetType;
      //   * map to hold the above information for each coupling entity
      typedef std::map< IndexType, IntersectionInfoSetType > EntityToIntersectionSetMapType;
      //   * map to hold the above map for each neighboring subdomain
      typedef std::map< unsigned int, Dune::shared_ptr< EntityToIntersectionSetMapType > > CouplingIntersectionMapType;
      //   * vector to hold the above map for each subdomain
      std::vector< CouplingIntersectionMapType > couplingBoundaryInfos(size_, CouplingIntersectionMapType());
      // for the boundary grid parts
      //   * a vector to hold the global entity ids for each subdomain
      typename std::vector< Dune::shared_ptr< GeometryMapType > > boundaryGeometryMaps(size_);
      //   * a vector to hold the codim sizes
      typename std::vector< CodimSizesType > boundaryCodimSizesVector(size_, CodimSizesType(0));
      //   * a vector to hold the intersection informations
      std::vector< Dune::shared_ptr< EntityToIntersectionSetMapType > > boundaryInfos(size_);
      // for the neighboring information
      neighboringSubdomainSets_
          = Dune::shared_ptr< std::vector< NeighboringSubdomainsSetType > >(new std::vector< NeighboringSubdomainsSetType >(size_, NeighboringSubdomainsSetType()));
      std::vector< NeighboringSubdomainsSetType >& neighboringSubdomainSets = *neighboringSubdomainSets_;

      // loop over all subdomains
      //   * to test for consecutive numbering
      //   * to compute the number of codim 0 entities per subdomain
      //   * to initialize data structures
      std::vector< unsigned int > subdomainSizes(size_, 0);
      for (unsigned int subdomain = 0; subdomain < size_; ++subdomain) {
        // test if this subdomain exists
        typename SubdomainMapType::iterator subdomainToEntityMapIt = subdomainToEntityMap_.find(subdomain);
        if (subdomainToEntityMapIt == subdomainToEntityMap_.end()) {
          std::stringstream msg;
          msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize()!";
          DUNE_THROW(Dune::InvalidStateException, msg.str());
        } else {
          // compute number of codim 0 entities
          const GeometryMapType& subdomainGeometryMap = *(subdomainToEntityMapIt->second);
          // loop over all geometry types of this subdomain
          for (typename GeometryMapType::const_iterator subdomainGeometryMapIt = subdomainGeometryMap.begin();
               subdomainGeometryMapIt != subdomainGeometryMap.end();
               ++subdomainGeometryMapIt) {
            // get geometry
            const GeometryType& geometryType = subdomainGeometryMapIt->first;
            // if this is a codim 0 geometry
            if (geometryType.dim() == dim) {
              // get the index map of this geometry type
              const IndexMapType& indexMap = subdomainGeometryMapIt->second;
              subdomainSizes[subdomain] += indexMap.size();
            } // if this is a codim 0 geometry
          } // loop over all geometry types of this subdomain
          // init data structures
          subdomainInnerBoundaryInfos[subdomain] = Dune::shared_ptr< EntityToIntersectionInfoMapType >(new EntityToIntersectionInfoMapType());
          boundaryGeometryMaps[subdomain] = Dune::shared_ptr< GeometryMapType >(new GeometryMapType());
          boundaryInfos[subdomain] = Dune::shared_ptr< EntityToIntersectionSetMapType >(new EntityToIntersectionSetMapType());
        } // test if this subdomain exists
      } // loop over all subdomains
      out << size_ << " subdomains:" << std::endl;

      // walk the global grid part
      //   * to generate the information which sudomains neighbor each other
      out << prefix << "  generating boundary information:" << std::endl;
      for (typename GlobalGridPartType::template Codim< 0 >::IteratorType entityIt = globalGridPart_->template begin< 0 >();
           entityIt != globalGridPart_->template end< 0 >();
           ++entityIt) {
        // find the subdomains this entity lives in
        const EntityType& entity = *entityIt;
        const IndexType entityGlobalIndex = globalGridPart_->indexSet().index(entity);
        const unsigned int entitySubdomain = getSubdomainOf(entityGlobalIndex);
        // get the set of this subdomains neighbors
        NeighboringSubdomainsSetType& neighborsOfSubdomain = neighboringSubdomainSets[entitySubdomain];
        // get the boundary info map for this subdomain
        EntityToIntersectionInfoMapType& subdomainInnerBoundaryInfo
            = *(subdomainInnerBoundaryInfos[entitySubdomain]);
        // walk the neighbors
        bool subdomainsEntitiesAreConnected = false;
        for (typename GlobalGridPartType::IntersectionIteratorType intersectionIt = globalGridPart_->ibegin(entity);
             intersectionIt != globalGridPart_->iend(entity);
             ++intersectionIt) {
          typedef typename GlobalGridPartType::IntersectionIteratorType::Intersection IntersectionType;
          const IntersectionType& intersection = *intersectionIt;
          // check the type of this intersection
          if (intersection.boundary() && !intersection.neighbor()) {
            // get local index of the intersection
            const int intersectionLocalIndex = intersection.indexInInside();
            // report
            out << prefix << "    entity " << entityGlobalIndex
                  << " lies at the domain boundary of subdomain " << entitySubdomain << std::endl;
            // for the boundary grid part
            //   * get the map for this subdomain
            GeometryMapType& boundaryGeometryMap = *(boundaryGeometryMaps[entitySubdomain]);
            //   * and add geometry and global index of this codim 0 entity
            const GeometryType& entityGeometryType = entity.type();
            CodimSizesType& boundaryCodimSizes = boundaryCodimSizesVector[entitySubdomain];
            addGeometryAndIndex(boundaryGeometryMap, boundaryCodimSizes, entityGeometryType, entityGlobalIndex);
            //   * and of all remaining codims
            Add< 1, dim >::subEntities(*this, entity, boundaryGeometryMap, boundaryCodimSizes);
            // for the intersection information of the boundary grid part
            //   * get the map for this subdomain
            EntityToIntersectionSetMapType& boundaryInfo = *(boundaryInfos[entitySubdomain]);
            //   * and get the entry for this entity
            IntersectionInfoSetType& entityBoundaryInfo = boundaryInfo[entityGlobalIndex];
            //   * and add this local intersection
            entityBoundaryInfo.insert(intersectionLocalIndex);
          } else if (intersection.neighbor()) {
            // then this entity lies inside the domain
            // and has a neighbor
            typedef typename IntersectionType::EntityPointer EntityPtrType;
            const EntityPtrType neighborPtr = intersection.outside();
            const EntityType& neighbor = *neighborPtr;
            const IndexType& neighborGlobalIndex = globalGridPart_->indexSet().index(neighbor);
            const unsigned int neighborSubdomain = getSubdomainOf(neighborGlobalIndex);
            // check if neighbor is in another same subdomain
            if (neighborSubdomain != entitySubdomain) {
              // get local index of the intersection
              const int intersectionLocalIndex = intersection.indexInInside();
              // report
              out << prefix << "    entity " << entityGlobalIndex
                    << " lies at an inner boundary of subdomain " << entitySubdomain << std::endl;
              out << prefix << "      at intersection " << intersectionLocalIndex
                    << " with neighbor " << neighborGlobalIndex
                    << " of subdomain " << neighborSubdomain << std::endl;
              // for the neighbor information between the subdomains
              //   * the subdomain of the neighbor is a neighboring subdomain of the entities subdomain
              neighborsOfSubdomain.insert(neighborSubdomain);
              // for the subdomain grid part
              //   * get the boundary info map for this entity
              IntersectionToBoundaryIdMapType& entityInnerBoundaryInfo = subdomainInnerBoundaryInfo[entityGlobalIndex];
              //   * and add the local intersection id and its desired fake boundary id to this entities map
              entityInnerBoundaryInfo.insert(std::pair< int, int >(intersectionLocalIndex, boundaryId_));
              // for the coupling grid part
              //   * get the coupling map for this subdomain
              SubdomainMapType& couplingMap = couplingMaps[entitySubdomain];
              //   * create an entry for the neighboring subdomain (if needed)
              if (couplingMap.find(neighborSubdomain) == couplingMap.end())
                couplingMap.insert(std::pair< unsigned int, Dune::shared_ptr< GeometryMapType > >(
                                     neighborSubdomain,
                                     Dune::shared_ptr< GeometryMapType >(new GeometryMapType)));
              //   * and get it
              GeometryMapType& couplingGeometryMap = *(couplingMap[neighborSubdomain]);
              //   * and add geometry and global index of this codim 0 entity
              const GeometryType& entityGeometryType = entity.type();
              std::map< unsigned int, CodimSizesType >& couplingCodimSizeMap = couplingCodimSizeMaps[entitySubdomain];
              //   * create codim sizes vector for this neighbor (if needed)
              if (couplingCodimSizeMap.find(neighborSubdomain) == couplingCodimSizeMap.end())
                couplingCodimSizeMap.insert(std::pair< unsigned int, CodimSizesType >(neighborSubdomain, CodimSizesType(0)));
              //  * and get it
              CodimSizesType& couplingCodimSizes = couplingCodimSizeMap[neighborSubdomain];
              addGeometryAndIndex(couplingGeometryMap, couplingCodimSizes, entityGeometryType, entityGlobalIndex);
              //   * and of all remaining codims
              Add< 1, dim >::subEntities(*this, entity, couplingGeometryMap, couplingCodimSizes);
              // for the intersection information of the coupling
              //   * get the map for this subdomain
              CouplingIntersectionMapType& couplingBoundaryInfo = couplingBoundaryInfos[entitySubdomain];
              //   * and create it for this neighbor (if needed)
              if (couplingBoundaryInfo.find(neighborSubdomain) == couplingBoundaryInfo.end())
                couplingBoundaryInfo.insert(
                      std::pair< unsigned int, Dune::shared_ptr< EntityToIntersectionSetMapType > >(
                        neighborSubdomain,
                        Dune::shared_ptr< EntityToIntersectionSetMapType >(new EntityToIntersectionSetMapType())));
              //   * and get it
              EntityToIntersectionSetMapType& couplingBoundaryInfoMap = *(couplingBoundaryInfo[neighborSubdomain]);
              //   * get the entry for this entity
              IntersectionInfoSetType& entityCouplingBoundaryInfo = couplingBoundaryInfoMap[entityGlobalIndex];
              //   * and add this local intersection
              entityCouplingBoundaryInfo.insert(intersectionLocalIndex);
            } else { // if neighbor is contained in this subdomain
              subdomainsEntitiesAreConnected = true;
            } // check if neighbor is in another subdomain
          } // check the type of this intersection
        } // walk the neighbors
        // check if this entity is connected to the other entities of this subdomain
        if (subdomainSizes[entitySubdomain] != 1 && !subdomainsEntitiesAreConnected) {
          std::stringstream msg;
          msg << "Error in " << id << ": at least one entity of subdomain " << entitySubdomain
              << " is not connected to entity " << entityGlobalIndex << " (connected)!";
          DUNE_THROW(Dune::InvalidStateException, msg.str());
        } // check if this entity is connected to the other entities of this subdomain
      } // walk the global grid part

      // walk the subdomains
      //   * to create the local grid parts
      localGridParts_ = Dune::shared_ptr< std::vector< Dune::shared_ptr< const LocalGridPartType > > >(
            new std::vector< Dune::shared_ptr< const LocalGridPartType > >(size_));
      std::vector< Dune::shared_ptr< const LocalGridPartType > >& localGridParts = *localGridParts_;
      //   * to create the boundary grid parts
      boundaryGridParts_ = Dune::shared_ptr< std::vector< Dune::shared_ptr< const BoundaryGridPartType > > >(
            new std::vector< Dune::shared_ptr< const BoundaryGridPartType > >(size_));
      std::vector< Dune::shared_ptr< const BoundaryGridPartType > >& boundaryGridParts = *boundaryGridParts_;
      out << prefix << "creating local and boundary grid parts:" << std::endl;
      for (typename SubdomainMapType::const_iterator subdomainIterator = subdomainToEntityMap_.begin();
           subdomainIterator != subdomainToEntityMap_.end();
           ++subdomainIterator) {
        // report
        const unsigned int subdomain = subdomainIterator->first;
        const unsigned int subdomainSize  = subdomainSizes[subdomain];
        out << prefix << "  subdomain " << subdomain << " (of size " << subdomainSize << ")... " << std::flush;
        // for the local grid part
        //   * get the geometry map
        const Dune::shared_ptr< const GeometryMapType > localGeometryMap = subdomainIterator->second;
        //   * get the boundary info map
        const Dune::shared_ptr< const EntityToIntersectionInfoMapType > localBoundaryInfo = subdomainInnerBoundaryInfos[subdomain];
        //   * and create the local grid part
        localGridParts[subdomain] = Dune::shared_ptr< const LocalGridPartType >(
              new LocalGridPartType(globalGridPart_,
                                    localGeometryMap,
                                    localBoundaryInfo));
        // for the boundary grid part
        //   * get the geometry map
        const Dune::shared_ptr< const GeometryMapType > boundaryGeometryMap = boundaryGeometryMaps[subdomain];
        //   * get the boundary info map
        const Dune::shared_ptr< const EntityToIntersectionSetMapType > boundaryBoundaryInfo = boundaryInfos[subdomain];
        //   * and create the boundary grid part
        boundaryGridParts[subdomain] = Dune::shared_ptr< const BoundaryGridPartType >(
              new BoundaryGridPartType(globalGridPart_,
                                       boundaryGeometryMap,
                                       boundaryBoundaryInfo,
                                       localGridParts[subdomain]));
        out << "done" << std::endl;
      } // walk the subdomains

      // walk the subdomains
      //   * to create the coupling grid parts
      couplingGridPartsMaps_ = Dune::shared_ptr< std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > > >(
            new std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > >(size_));
      std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > >& couplingGridPartsMaps = *couplingGridPartsMaps_;
      out << prefix << "creating coupling grid parts:" << std::endl;
      for (unsigned int subdomain = 0; subdomain < couplingMaps.size(); ++subdomain) {
        // get the coupling map for this subdomain
        const SubdomainMapType& couplingMap = couplingMaps[subdomain];
        // get the intersection information map for this subdomain
        const CouplingIntersectionMapType& couplingBoundaryInfo = couplingBoundaryInfos[subdomain];
        // get the target map for this subdomain
        std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > >& couplingGridPartsMap = couplingGridPartsMaps[subdomain];
        // loop over all neighbors
        for (typename SubdomainMapType::const_iterator neighborIt = couplingMap.begin();
             neighborIt != couplingMap.end();
             ++neighborIt) {
          const unsigned int neighbor = neighborIt->first;
          out << prefix << "  subdomain " << subdomain << ", neighbor " << neighbor <<  "... " << std::flush;
          // get the geometry map
          const Dune::shared_ptr< const GeometryMapType > couplingGeometryMap = neighborIt->second;
          // get the boundary info map
          typename CouplingIntersectionMapType::const_iterator result = couplingBoundaryInfo.find(neighbor);
          assert(result != couplingBoundaryInfo.end() && "This should not happen (see above)!");
          const Dune::shared_ptr< const EntityToIntersectionSetMapType > couplingBoundaryInfo = result->second;
          // and create the coupling grid part
          couplingGridPartsMap.insert(std::pair< unsigned int, Dune::shared_ptr< const CouplingGridPartType > >(
                                        neighbor,
                                        Dune::shared_ptr< const CouplingGridPartType >(
                                          new CouplingGridPartType(globalGridPart_,
                                                                   couplingGeometryMap,
                                                                   couplingBoundaryInfo,
                                                                   localGridParts[subdomain],
                                                                   localGridParts[neighbor]))));
          out << "done" << std::endl;
        } // loop over all neighbors
      } // walk the subdomains

      // done
      finalized_ = true;
    } // if (!finalized_)
  } // void finalize()

  const Dune::shared_ptr< const MsGridType > createMsGrid() const
  {
    assert(finalized_ && "Please call finalize() before calling createMsGrid()!");
    const Dune::shared_ptr< const MsGridType > msGrid(new MsGridType(grid_,
                                                                     globalGridPart_,
                                                                     size_,
                                                                     neighboringSubdomainSets_,
                                                                     entityToSubdomainMap_,
                                                                     localGridParts_,
                                                                     boundaryGridParts_,
                                                                     couplingGridPartsMaps_));
    return msGrid;
  } // const Dune::shared_ptr< const MsGridType > createMsGrid() const

private:
  void addGeometryAndIndex(GeometryMapType& geometryMap,
                           CodimSizesType& localCodimSizes,
                           const GeometryType& geometryType,
                           const IndexType& globalIndex,
                           const std::string& prefix = "",
                           std::ostream& out = Dune::Stuff::Common::Logger().devnull())
  {
    // get the map to this geometry type
    typename GeometryMapType::mapped_type& indexMap = geometryMap[geometryType];
    // add if needed
    if (indexMap.find(globalIndex) == indexMap.end()) {
      const unsigned int codim = dim - geometryType.dim();
      const IndexType localIndex = localCodimSizes[codim];
      indexMap.insert(std::pair< IndexType, IndexType >(globalIndex, localIndex));
      // increase count for this codim
      ++(localCodimSizes[codim]);
      // report
      out << prefix << "- added " << geometryType << " with global index " << globalIndex << " and local index " << localIndex << std::endl;
    } // add if needed
  } // void addGeometryAndIndex()

  unsigned int getSubdomainOf(const IndexType& globalIndex) const
  {
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "Error in " << id << ": entity " << globalIndex << " not added to any subdomain!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // unsigned int getSubdomainOf(const IndexType& globalIndex) const

  // friends
  template< int, int >
  friend class Add;

  // members
  const Dune::shared_ptr< const GridType > grid_;
  const int boundaryId_;
  bool prepared_;
  bool finalized_;
  unsigned int size_;
  Dune::shared_ptr< const GlobalGridPartType > globalGridPart_;
  // for the entity <-> subdomain relations
  Dune::shared_ptr< EntityToSubdomainMapType > entityToSubdomainMap_;
  SubdomainMapType subdomainToEntityMap_;
  // for the neighboring information
  Dune::shared_ptr< std::vector< NeighboringSubdomainsSetType > > neighboringSubdomainSets_;
  // for the local grid parts
  std::map< unsigned int, CodimSizesType > localCodimSizes_;
  Dune::shared_ptr< std::vector< Dune::shared_ptr< const LocalGridPartType > > > localGridParts_;
  // for the boundary grid parts
  Dune::shared_ptr< std::vector< Dune::shared_ptr< const BoundaryGridPartType > > > boundaryGridParts_;
  // for the coupling grid parts
  Dune::shared_ptr< std::vector< std::map< unsigned int, Dune::shared_ptr< const CouplingGridPartType > > > > couplingGridPartsMaps_;
}; // class Default

template< class GridType >
const std::string Default< GridType >::id = "grid.multiscale.factory.default";

//! specialization to stop the recursion
template< class GridType >
template< int c >
struct Default< GridType >::Add< c, c >
{
  static void subEntities(Default< GridType >& factory,
                          const typename Default< GridType >::EntityType& entity,
                          typename Default< GridType >::GeometryMapType& geometryMap,
                          typename Default< GridType >::CodimSizesType& localCodimSizes)
  {
    // suppress output, since we are not codim 0
    Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
    // loop over all codim c subentities of this entity
    typedef typename Default< GridType >::EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
    for (int i = 0; i < entity.template count< c >(); ++i) {
      const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
      const Default< GridType >::GeometryType& geometryType = codimCentityPtr->type();
      const typename Default< GridType >::IndexType globalIndex = factory.globalGridPart_->indexSet().index(*codimCentityPtr);
      factory.addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex, "", devnull);
    } // loop over all codim c subentities of this entity
  } // static void subEntities()
}; // struct Default< GridType >::Add< c, c >

} // namespace Factory

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_FACTORY_DEFAULT_HH
