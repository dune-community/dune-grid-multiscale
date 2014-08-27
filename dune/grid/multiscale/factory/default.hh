// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_FACTORY_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_FACTORY_DEFAULT_HH

#include <memory>
#include <vector>
#include <map>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/sgrid.hh>
#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#include <dune/grid/part/leaf.hh>
#include <dune/grid/part/local/indexbased.hh>
#include <dune/grid/multiscale/default.hh>

#include <dune/stuff/common/logging.hh>

namespace Dune {
namespace grid {
namespace Multiscale {
namespace Factory {


template< class GridImp >
class NeighborRecursionLevel
{
  static_assert(AlwaysFalse< GridImp >::value, "Please add an appropriate specialization for this GridImp!");
public:
  static size_t compute() = delete;
};

template<>
class NeighborRecursionLevel< SGrid< 2, 2 > >
{
public:
  static size_t compute()
  {
    return 1;
  }
};

#if HAVE_ALUGRID
template<>
class NeighborRecursionLevel< ALUConformGrid< 2, 2 > >
{
public:
  static size_t compute()
  {
    return 3;
  }
};

template<>
class NeighborRecursionLevel< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming, Dune::No_Comm > >
{
public:
  static size_t compute()
  {
    return 3;
  }
};

template<>
class NeighborRecursionLevel< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming, Dune::No_Comm > >
{
public:
  static size_t compute()
  {
    return 1;
  }
};
#endif


template< class GridImp >
class Default
{
public:
  typedef GridImp GridType;

  typedef Default< GridType > ThisType;

  static const unsigned int dim = GridType::dimension;

  typedef Dune::grid::Multiscale::Default< GridType > MsGridType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

  static const std::string id()
  {
    return "grid.multiscale.factory.default";
  }

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
  typedef std::map< size_t, std::shared_ptr< GeometryMapType > > SubdomainMapType;

  // map type which maps from an entity index (of the global grid parts index set) to a subdomain
  typedef std::map< IndexType, size_t > EntityToSubdomainMapType;

  typedef Dune::FieldVector< size_t, dim + 1 > CodimSizesType;

  // for the neighbor information between the subdomains
  typedef std::set< size_t > NeighboringSubdomainsSetType;

  template< int c, int d >
  struct Add
  {
    static void subEntities(ThisType& factory,
                            const EntityType& entity,
                            GeometryMapType& geometryMap,
                            CodimSizesType& localCodimSizes,
                            const std::string prefix = "",
                            std::ostream& out = Dune::Stuff::Common::Logger().devnull())
    {
//      // suppress output, since we are not codim 0
//      Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
      // loop over all codim c subentities of the entity
      typedef typename EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
      for (int i = 0; i < entity.template count< c >(); ++i) {
        const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
        const GeometryType& geometryType = codimCentityPtr->type();
        const IndexType globalIndex = factory.globalGridPart_->indexSet().index(*codimCentityPtr);
        factory.addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex, prefix, out);
      } // loop over all codim c subentities of the entity
      // add all codim c + 1 subentities
      Add< c + 1, d >::subEntities(factory, entity, geometryMap, localCodimSizes, prefix, out);
    } // static void subEntities()
  }; // struct Add

public:
  Default(const GridType& grid, const int boundaryId = 7)
    : grid_(Dune::stackobject_to_shared_ptr(grid))
    , boundaryId_(boundaryId)
    , prepared_(false)
    , finalized_(false)
    , size_(0)
    , oversampled_(false)
  {}

  Default(const std::shared_ptr< const GridType > grid, const int boundaryId = 7)
    : grid_(grid)
    , boundaryId_(boundaryId)
    , prepared_(false)
    , finalized_(false)
    , size_(0)
    , oversampled_(false)
  {}

  void prepare()
  {
    if (!prepared_) {
      globalGridPart_ = std::make_shared< const GlobalGridPartType >(const_cast< GridType& >(*grid_));
      entityToSubdomainMap_ = std::shared_ptr< EntityToSubdomainMapType >(new EntityToSubdomainMapType());
      prepared_ = true;
    } // if (!prepared_)
  } // void prepare()

  const std::shared_ptr< const GlobalGridPartType > globalGridPart() const
  {
    assert(prepared_ && "Please call prepare() before calling globalGridPart()!");
    return globalGridPart_;
  }

  void add(const EntityType& entity,
           const size_t subdomain,
           const std::string prefix = "",
           std::ostream& out = Dune::Stuff::Common::Logger().debug())
  {
    // prepare
    assert(prepared_ && "Please call prepare() before calling add()!");
    assert(!finalized_ && "Do not call add() after calling finalized()!");
    const IndexType globalIndex = globalGridPart_->indexSet().index(entity);
#ifndef NDEBUG
    out << prefix << "adding entity " << globalIndex << " to subdomain " << subdomain << ":" << std::endl;
#endif
    // add subdomain to this entity index
    typename EntityToSubdomainMapType::iterator indexIt = entityToSubdomainMap_->find(globalIndex);
    if (indexIt == entityToSubdomainMap_->end()) {
      entityToSubdomainMap_->insert(std::pair< IndexType, size_t >(globalIndex, subdomain));
    } else {
      if (indexIt->second != subdomain) {
        std::stringstream msg;
        msg << "Error in " << id()<< ": can not add entity to more than one subdomain!";
        DUNE_THROW(Dune::InvalidStateException, msg.str());
      }
    } // add subdomain to this entity index
    // create geometry map for this subdomain if needed (doing this explicitly (instead of just using insert()) only to increment size)
    if (subdomainToEntityMap_.find(subdomain) == subdomainToEntityMap_.end()) {
      subdomainToEntityMap_.insert(std::pair< size_t, std::shared_ptr< GeometryMapType > >(
          subdomain,
          std::shared_ptr< GeometryMapType >(new GeometryMapType())));
      ++size_;
    } // create geometry map for this subdomain if needed
    // create local codim sizes for this subdomain
    localCodimSizes_.insert(std::pair< size_t, CodimSizesType >(
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
    Add< 1, dim >::subEntities(*this, entity, geometryMap, localCodimSizes, prefix, out);
  } // void add()

  void finalize(const size_t oversamplingLayers = 0,
                const size_t neighbor_recursion_level = NeighborRecursionLevel< GridType >::compute(),
                const std::string prefix = "",
                std::ostream& out = Dune::Stuff::Common::Logger().debug(),
                bool assert_connected = true)
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
      std::vector< std::shared_ptr< EntityToIntersectionInfoMapType > > subdomainInnerBoundaryInfos(size_);
      // for the coupling grid parts
      //   * vector to hold a (neighboring subdomain -> entity) map for each subdomain
      typename std::vector< SubdomainMapType > couplingMaps(size_, SubdomainMapType());
      //   * vector to hold a map of coupling sizes
      std::vector< std::map< size_t, CodimSizesType > > couplingCodimSizeMaps(size_, std::map< size_t, CodimSizesType >());
      //   * set of local coupling intersections
      typedef std::vector< int > IntersectionInfoSetType;
      //   * map to hold the above information for each coupling entity
      typedef std::map< IndexType, IntersectionInfoSetType > EntityToIntersectionSetMapType;
      //   * map to hold the above map for each neighboring subdomain
      typedef std::map< size_t, std::shared_ptr< EntityToIntersectionSetMapType > > CouplingIntersectionMapType;
      //   * vector to hold the above map for each subdomain
      std::vector< CouplingIntersectionMapType > couplingBoundaryInfos(size_, CouplingIntersectionMapType());
      // for the boundary grid parts
      //   * a map to hold the global entity ids for each subdomain
      typename std::map< size_t, std::shared_ptr< GeometryMapType > > boundaryGeometryMapMap;
      //   * a map to hold the codim sizes
      typename std::map< size_t, CodimSizesType > boundaryCodimSizesMap;
      //   * a map to hold the intersection informations
      std::map< size_t, std::shared_ptr< EntityToIntersectionSetMapType > > boundaryInfoMap;
      // for the neighboring information
      neighboringSubdomainSets_
          = std::shared_ptr< std::vector< NeighboringSubdomainsSetType > >(new std::vector< NeighboringSubdomainsSetType >(size_, NeighboringSubdomainsSetType()));
      std::vector< NeighboringSubdomainsSetType >& neighboringSubdomainSets = *neighboringSubdomainSets_;
      // loop over all subdomains
      //   * to test for consecutive numbering
      //   * to compute the number of codim 0 entities per subdomain
      //   * to initialize data structures
      std::vector< size_t > subdomainSizes(size_, 0);
      for (size_t subdomain = 0; subdomain < size_; ++subdomain) {
        // test if this subdomain exists
        typename SubdomainMapType::iterator subdomainToEntityMapIt = subdomainToEntityMap_.find(subdomain);
        if (subdomainToEntityMapIt == subdomainToEntityMap_.end()) {
          std::stringstream msg;
          msg << "Error in " << id()<< ": numbering of subdomains has to be consecutive upon calling finalize()!";
          DUNE_THROW(InvalidStateException, msg.str());
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
          subdomainInnerBoundaryInfos[subdomain] = std::shared_ptr< EntityToIntersectionInfoMapType >(new EntityToIntersectionInfoMapType());
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
        const size_t entitySubdomain = getSubdomainOf(entityGlobalIndex);
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
#ifndef NDEBUG
            out << prefix << "    entity " << entityGlobalIndex
                  << " lies at the domain boundary of subdomain " << entitySubdomain << std::endl;
#endif
            // for the boundary grid part
            //   * get the maps for this subdomain (and create them, if necessary)
            if (boundaryGeometryMapMap.find(entitySubdomain) == boundaryGeometryMapMap.end())
              boundaryGeometryMapMap.insert(std::pair< size_t, std::shared_ptr< GeometryMapType > >(entitySubdomain, std::shared_ptr< GeometryMapType >(new GeometryMapType())));
            GeometryMapType& boundaryGeometryMap = *(boundaryGeometryMapMap[entitySubdomain]);
            if (boundaryCodimSizesMap.find(entitySubdomain) == boundaryCodimSizesMap.end())
              boundaryCodimSizesMap.insert(std::pair< size_t, CodimSizesType >(entitySubdomain, CodimSizesType(0)));
            CodimSizesType& boundaryCodimSizes = boundaryCodimSizesMap[entitySubdomain];
            //   * and add geometry and global index of this codim 0 entity
            const GeometryType& entityGeometryType = entity.type();
            addGeometryAndIndex(boundaryGeometryMap, boundaryCodimSizes, entityGeometryType, entityGlobalIndex);
            //   * and of all remaining codims
            Add< 1, dim >::subEntities(*this, entity, boundaryGeometryMap, boundaryCodimSizes);
            // for the intersection information of the boundary grid part
            //   * get the map for this subdomain (and create it, if necessary)
            if (boundaryInfoMap.find(entitySubdomain) == boundaryInfoMap.end())
              boundaryInfoMap.insert(std::pair< size_t, std::shared_ptr< EntityToIntersectionSetMapType > >(entitySubdomain, std::shared_ptr< EntityToIntersectionSetMapType >(new EntityToIntersectionSetMapType())));
            EntityToIntersectionSetMapType& boundaryInfo = *(boundaryInfoMap[entitySubdomain]);
            //   * and get the entry for this entity
            IntersectionInfoSetType& entityBoundaryInfo = boundaryInfo[entityGlobalIndex];
            //   * and add this local intersection
            entityBoundaryInfo.push_back(intersectionLocalIndex);
          } else if (intersection.neighbor()) {
            // then this entity lies inside the domain
            // and has a neighbor
            typedef typename IntersectionType::EntityPointer EntityPtrType;
            const EntityPtrType neighborPtr = intersection.outside();
            const EntityType& neighbor = *neighborPtr;
            const IndexType& neighborGlobalIndex = globalGridPart_->indexSet().index(neighbor);
            const size_t neighborSubdomain = getSubdomainOf(neighborGlobalIndex);
            // check if neighbor is in another or in the same subdomain
            if (neighborSubdomain != entitySubdomain) {
              // get local index of the intersection
              const int intersectionLocalIndex = intersection.indexInInside();
#ifndef NDEBUG
              // report
              out << prefix << "    entity " << entityGlobalIndex
                    << " lies at an inner boundary of subdomain " << entitySubdomain << std::endl;
              out << prefix << "      at intersection " << intersectionLocalIndex
                    << " with neighbor " << neighborGlobalIndex
                    << " of subdomain " << neighborSubdomain << std::endl;
#endif
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
                couplingMap.insert(std::pair< size_t, std::shared_ptr< GeometryMapType > >(
                                     neighborSubdomain,
                                     std::shared_ptr< GeometryMapType >(new GeometryMapType)));
              //   * and get it
              GeometryMapType& couplingGeometryMap = *(couplingMap[neighborSubdomain]);
              //   * and add geometry and global index of this codim 0 entity
              const GeometryType& entityGeometryType = entity.type();
              std::map< size_t, CodimSizesType >& couplingCodimSizeMap = couplingCodimSizeMaps[entitySubdomain];
              //   * create codim sizes vector for this neighbor (if needed)
              if (couplingCodimSizeMap.find(neighborSubdomain) == couplingCodimSizeMap.end())
                couplingCodimSizeMap.insert(std::pair< size_t, CodimSizesType >(neighborSubdomain, CodimSizesType(0)));
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
                      std::pair< size_t, std::shared_ptr< EntityToIntersectionSetMapType > >(
                        neighborSubdomain,
                        std::shared_ptr< EntityToIntersectionSetMapType >(new EntityToIntersectionSetMapType())));
              //   * and get it
              EntityToIntersectionSetMapType& couplingBoundaryInfoMap = *(couplingBoundaryInfo[neighborSubdomain]);
              //   * get the entry for this entity
              IntersectionInfoSetType& entityCouplingBoundaryInfo = couplingBoundaryInfoMap[entityGlobalIndex];
              //   * and add this local intersection
              entityCouplingBoundaryInfo.push_back(intersectionLocalIndex);
            } else { // if neighbor is contained in this subdomain
              subdomainsEntitiesAreConnected = true;
            } // check if neighbor is in another subdomain
          } // check the type of this intersection
        } // walk the neighbors
        // check if this entity is connected to the other entities of this subdomain
        if (assert_connected && subdomainSizes[entitySubdomain] != 1 && !subdomainsEntitiesAreConnected) {
          std::stringstream msg;
          msg << "Error in " << id()<< ": at least one entity of subdomain " << entitySubdomain
              << " is not connected to entity " << entityGlobalIndex << " (connected)!";
          DUNE_THROW(Dune::InvalidStateException, msg.str());
        } // check if this entity is connected to the other entities of this subdomain
      } // walk the global grid part
      // walk the subdomains
      //   * to create the local grid parts
      localGridParts_ = std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > >(
            new std::vector< std::shared_ptr< const LocalGridPartType > >(size_));
      std::vector< std::shared_ptr< const LocalGridPartType > >& localGridParts = *localGridParts_;
      out << prefix << "creating local grid parts:" << std::endl;
      for (typename SubdomainMapType::const_iterator subdomainIterator = subdomainToEntityMap_.begin();
           subdomainIterator != subdomainToEntityMap_.end();
           ++subdomainIterator) {
        // report
        const size_t subdomain = subdomainIterator->first;
        const size_t subdomainSize  = subdomainSizes[subdomain];
#ifndef NDEBUG
        out << prefix << "  subdomain " << subdomain << " (of size " << subdomainSize << ")... " << std::flush;
#endif
        // for the local grid part
        //   * get the geometry map
        const std::shared_ptr< const GeometryMapType > localGeometryMap = subdomainIterator->second;
        //   * get the boundary info map
        const std::shared_ptr< const EntityToIntersectionInfoMapType > localBoundaryInfo = subdomainInnerBoundaryInfos[subdomain];
        //   * and create the local grid part
        localGridParts[subdomain] = std::shared_ptr< const LocalGridPartType >(
              new LocalGridPartType(globalGridPart_,
                                    localGeometryMap,
                                    localBoundaryInfo));
#ifndef NDEBUG
        out << "done" << std::endl;
#endif
      } // walk the subdomains

      // walk those subdomains which have a boundary grid part
      out << prefix << "creating boundary grid parts:" << std::endl;
      //   * to create the boundary grid parts
      boundaryGridParts_ = std::shared_ptr< std::map< size_t, std::shared_ptr< const BoundaryGridPartType > > >(
            new std::map< size_t, std::shared_ptr< const BoundaryGridPartType > >());
      std::map< size_t, std::shared_ptr< const BoundaryGridPartType > >& boundaryGridParts = *boundaryGridParts_;
      typename std::map< size_t, CodimSizesType >::const_iterator boundaryCodimSizesMapIt = boundaryCodimSizesMap.begin();
      typename std::map< size_t, std::shared_ptr< EntityToIntersectionSetMapType > >::const_iterator boundaryInfoMapIt = boundaryInfoMap.begin();
      for (typename std::map< size_t, std::shared_ptr< GeometryMapType > >::const_iterator boundaryGeometryMapMapIt = boundaryGeometryMapMap.begin();
           boundaryGeometryMapMapIt != boundaryGeometryMapMap.end();
           ++boundaryGeometryMapMapIt,
           ++boundaryCodimSizesMapIt,
           ++boundaryInfoMapIt) {
        const size_t boundarySubdomain = boundaryGeometryMapMapIt->first;
        assert(boundarySubdomain == boundaryCodimSizesMapIt->first && "We should not get here: we are in big trouble, if these maps do not correspond to each other!");
        assert(boundarySubdomain == boundaryInfoMapIt->first && "We should not get here: we are in big trouble, if these maps do not correspond to each other!");
#ifndef NDEBUG
        out << prefix << "  subdomain " << boundarySubdomain << " (of size " << boundaryCodimSizesMapIt->second.operator[](0) << ")... " << std::flush;
#endif
        // for the boundary grid part
        //   * get the geometry map
        const std::shared_ptr< const GeometryMapType > boundaryGeometryMap = boundaryGeometryMapMapIt->second;
        //   * get the boundary info map
        const std::shared_ptr< const EntityToIntersectionSetMapType > boundaryBoundaryInfo = boundaryInfoMapIt->second;
        //   * and create the boundary grid part
        boundaryGridParts.insert(std::pair< size_t, std::shared_ptr< const BoundaryGridPartType > >(
                                   boundarySubdomain,
                                   std::shared_ptr< const BoundaryGridPartType >(
                                     new BoundaryGridPartType(globalGridPart_,
                                                              boundaryGeometryMap,
                                                              boundaryBoundaryInfo,
                                                              localGridParts[boundarySubdomain]))));
#ifndef NDEBUG
        out << "done" << std::endl;
#endif
      } // walk those subdomains which have a boundary grid part
      // walk the subdomains
      //   * to create the coupling grid parts
      couplingGridPartsMaps_ = std::shared_ptr< std::vector< std::map< size_t, std::shared_ptr< const CouplingGridPartType > > > >(
            new std::vector< std::map< size_t, std::shared_ptr< const CouplingGridPartType > > >(size_));
      std::vector< std::map< size_t, std::shared_ptr< const CouplingGridPartType > > >& couplingGridPartsMaps = *couplingGridPartsMaps_;
      out << prefix << "creating coupling grid parts:" << std::endl;
      for (size_t subdomain = 0; subdomain < couplingMaps.size(); ++subdomain) {
        // get the coupling map for this subdomain
        const SubdomainMapType& couplingMap = couplingMaps[subdomain];
        // get the intersection information map for this subdomain
        const CouplingIntersectionMapType& couplingBoundaryInfo = couplingBoundaryInfos[subdomain];
        // get the target map for this subdomain
        std::map< size_t, std::shared_ptr< const CouplingGridPartType > >& couplingGridPartsMap = couplingGridPartsMaps[subdomain];
        // loop over all neighbors
        for (typename SubdomainMapType::const_iterator neighborIt = couplingMap.begin();
             neighborIt != couplingMap.end();
             ++neighborIt) {
          const size_t neighbor = neighborIt->first;
#ifndef NDEBUG
          out << prefix << "  subdomain " << subdomain << ", neighbor " << neighbor <<  "... " << std::flush;
#endif
          // get the geometry map
          const std::shared_ptr< const GeometryMapType > couplingGeometryMap = neighborIt->second;
          // get the boundary info map
          typename CouplingIntersectionMapType::const_iterator result = couplingBoundaryInfo.find(neighbor);
          assert(result != couplingBoundaryInfo.end() && "This should not happen (see above)!");
          const std::shared_ptr< const EntityToIntersectionSetMapType > coupling_boundary_info = result->second;
          // and create the coupling grid part
          couplingGridPartsMap.insert(std::pair< size_t, std::shared_ptr< const CouplingGridPartType > >(
                                        neighbor,
                                        std::shared_ptr< const CouplingGridPartType >(
                                          new CouplingGridPartType(globalGridPart_,
                                                                   couplingGeometryMap,
                                                                   coupling_boundary_info,
                                                                   localGridParts[subdomain],
                                                                   localGridParts[neighbor]))));
#ifndef NDEBUG
          out << "done" << std::endl;
#endif
        } // loop over all neighbors
      } // walk the subdomains

      // create the first layer of oversampling
      if (oversamplingLayers > 0) {
        oversampledLocalGridParts_ = addOneLayerOfOverSampling(subdomainToEntityMap_,
                                                               *localGridParts_,
                                                               neighbor_recursion_level,
                                                               subdomainToOversamplingEntitiesMap_);
        oversampled_ = true;
      }
      // and the rest
      for (size_t ii = 1; ii < oversamplingLayers; ++ii) {
        std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > >
            tmpOversapledGridParts = oversampledLocalGridParts_;
        oversampledLocalGridParts_ = std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > >(
                                       new std::vector< std::shared_ptr< const LocalGridPartType > >(size_));
        SubdomainMapType tmpSubdomainToOversamplingEntitiesMap = subdomainToOversamplingEntitiesMap_;
        subdomainToOversamplingEntitiesMap_ = SubdomainMapType();
        oversampledLocalGridParts_ = addOneLayerOfOverSampling(tmpSubdomainToOversamplingEntitiesMap,
                                                               *tmpOversapledGridParts,
                                                               neighbor_recursion_level,
                                                               subdomainToOversamplingEntitiesMap_);
      }

      // done
      finalized_ = true;
    } // if (!finalized_)
  } // void finalize()

  const std::shared_ptr< const MsGridType > createMsGrid() const
  {
    assert(finalized_ && "Please call finalize() before calling createMsGrid()!");
    if (oversampled_)
      return Dune::make_shared< MsGridType >(grid_,
                                             globalGridPart_,
                                             size_,
                                             neighboringSubdomainSets_,
                                             entityToSubdomainMap_,
                                             localGridParts_,
                                             boundaryGridParts_,
                                             couplingGridPartsMaps_,
                                             oversampledLocalGridParts_);
    else
      return Dune::make_shared< MsGridType >(grid_,
                                             globalGridPart_,
                                             size_,
                                             neighboringSubdomainSets_,
                                             entityToSubdomainMap_,
                                             localGridParts_,
                                             boundaryGridParts_,
                                             couplingGridPartsMaps_);
  } // const std::shared_ptr< const MsGridType > createMsGrid() const

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
      const size_t codim = dim - geometryType.dim();
      const IndexType localIndex = boost::numeric_cast< IndexType >(localCodimSizes[codim]);
      indexMap.insert(std::pair< IndexType, IndexType >(globalIndex, localIndex));
      // increase count for this codim
      ++(localCodimSizes[codim]);
      // report
#ifndef NDEBUG
      out << prefix << "- " << geometryType << ", global index " << globalIndex << ", local index " << localIndex << std::endl;
#endif
    } // add if needed
  } // void addGeometryAndIndex()

  size_t getSubdomainOf(const IndexType& globalIndex) const
  {
    const typename EntityToSubdomainMapType::const_iterator result = entityToSubdomainMap_->find(globalIndex);
    if (result == entityToSubdomainMap_->end()) {
      std::stringstream msg;
      msg << "Error in " << id()<< ": entity " << globalIndex << " not added to any subdomain!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }
    return result->second;
  } // size_t getSubdomainOf(const IndexType& globalIndex) const

      std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > >
  addOneLayerOfOverSampling(const SubdomainMapType& subdomainToEntityMap,
                            const std::vector< std::shared_ptr< const LocalGridPartType > >& localGridParts,
                            const size_t neighbor_recursion_level,
                            SubdomainMapType& subdomainToOversamplingEntitiesMap)
  {
    // init data structures
    // for the subdomains inner boundaries
    //   * to map the local intersection index to the desired fake boundary id
    typedef std::map< int, int > IntersectionToBoundaryIdMapType;
    //   * to map the global entity index to one of those maps
    typedef std::map< IndexType, IntersectionToBoundaryIdMapType > EntityToIntersectionInfoMapType;
    //   * to hold one of those maps for each subdomain
    std::vector< std::shared_ptr< EntityToIntersectionInfoMapType > > oversamplingSubdomainInnerBoundaryInfos(size_);
    // walk the subdomains to create the oversampling
    for (auto subdomainIt = subdomainToEntityMap.begin();
         subdomainIt != subdomainToEntityMap.end();
         ++subdomainIt) {
      const size_t subdomain = subdomainIt->first;
      const GeometryMapType& geometryMap = *(subdomainIt->second);
      // * therefore, hardcopy the existing map we want to extend,
      std::shared_ptr< GeometryMapType > geometryMapCopy = Dune::make_shared< GeometryMapType >(geometryMap);
      // * and create an empty local boundary info map for later use
      oversamplingSubdomainInnerBoundaryInfos[subdomain] = Dune::make_shared< EntityToIntersectionInfoMapType >();
      // * then walk the local grid part to find the local boundary entities
      const LocalGridPartType& localGridPart = *(localGridParts[subdomain]);
      for (auto entityIt = localGridPart.template begin< 0 >();
           entityIt != localGridPart.template end< 0 >();
           ++entityIt) {
        // get the entity index
        const EntityType& entity = *entityIt;
//        const IndexType entityGlobalIndex = globalGridPart_->indexSet().index(entity);
        // lets see if this is a boundary entity of the local grid part
        bool isOnLocalBoundary = false;
        for (auto intersectionIt = localGridPart.ibegin(entity);
             intersectionIt != localGridPart.iend(entity);
             ++intersectionIt) {
          const auto& intersection = *intersectionIt;
          if (intersection.boundary())
            isOnLocalBoundary = true;
        }
        if (isOnLocalBoundary) {
          // add all the "neighbors"
          // * therefore, iterate over the intersections in the global grid part
          for (auto intersectionIt = globalGridPart_->ibegin(entity);
               intersectionIt != globalGridPart_->iend(entity);
               ++intersectionIt) {
            const auto& intersection = *intersectionIt;
            // if this intersection is not on the domain boundary
            if (intersection.neighbor()) {
              // get the neighbor
              const auto neighborPtr = intersection.outside();
              const auto& neighbor = *neighborPtr;
              const IndexType& neighborGlobalIndex = globalGridPart_->indexSet().index(neighbor);
//              const size_t neighborSubdomain = getSubdomainOf(neighborGlobalIndex);
              // if the neighbor is not in the subdomain
              bool neighborIsNotInThisSubdomain = true;
              if (geometryMap.find(neighbor.type()) != geometryMap.end()) {
                if (geometryMap.find(neighbor.type())->second.find(neighborGlobalIndex)
                    != geometryMap.find(neighbor.type())->second.end()) {
                  neighborIsNotInThisSubdomain = false;
                }
              }
              if (neighborIsNotInThisSubdomain) {
                // add him to the oversampling
                // * therefore we can use the old localCodimSizes,
                CodimSizesType& localCodimSizes = localCodimSizes_.find(subdomain)->second;
                // * add the neighbor
                addGeometryAndIndex(*geometryMapCopy, localCodimSizes, neighbor.type(), neighborGlobalIndex);
                // * and all remaining codims entities
                Add< 1, dim >::subEntities(*this, neighbor, *geometryMapCopy, localCodimSizes);
                // and also check all its neighbors
                if (neighbor_recursion_level > 0)
                  add_neighbors_neighbors_recursively(entity,
                                                      neighbor,
                                                      geometryMap,
                                                      neighbor_recursion_level,
                                                      localCodimSizes,
                                                      *geometryMapCopy);
              } // if the neighbor is not in the subdomain
            } // if this intersection is not on the domain boundary
          } // iterate over the intersections in the global grid part
        } // lets see if this is a boundary entity of the local grid part
      } // then walk the local grid part to find the local boundary entities
      subdomainToOversamplingEntitiesMap.insert(std::make_pair(subdomain, geometryMapCopy));
    } // walk the subdomains to create the oversampling

    // now we need to create the local boundary info for the oversampling, so walk the global grid part
    for (auto entityIt = globalGridPart_->template begin< 0 >();
         entityIt != globalGridPart_->template end< 0 >();
         ++entityIt) {
      const auto& entity = *entityIt;
      const IndexType entityIndex = globalGridPart_->indexSet().index(entity);
      // now we find all the oversampled subdomains this entity is a part of
      for (auto subdomainToOversamplingEntitiesMapIt : subdomainToOversamplingEntitiesMap) {
        const auto geometryMapIt = subdomainToOversamplingEntitiesMapIt.second->find(entity.type());
        if (geometryMapIt != subdomainToOversamplingEntitiesMapIt.second->end()) {
          const auto& geometryMap = geometryMapIt->second;
          if (geometryMap.find(entityIndex) != geometryMap.end()) {
            // this entity is a part of this subdomain!
            const size_t entitySubdomain = subdomainToOversamplingEntitiesMapIt.first;
            // then walk the neighbors
            for (auto intersectionIt = globalGridPart_->ibegin(entity);
                 intersectionIt != globalGridPart_->iend(entity);
                 ++intersectionIt) {
              const auto& intersection = *intersectionIt;
              if (intersection.neighbor()) {
                const auto neighborPtr = intersection.outside();
                const auto& neighbor = *neighborPtr;
                const IndexType neighborIndex = globalGridPart_->indexSet().index(neighbor);
                // and check, if the neighbor is in the same subdomain
                bool isInSame = false;
                const auto neighborGeometryMapIt = subdomainToOversamplingEntitiesMapIt.second->find(neighbor.type());
                if (neighborGeometryMapIt != subdomainToOversamplingEntitiesMapIt.second->end()) {
                  const auto& neighborGeometryMap = neighborGeometryMapIt->second;
                  if (neighborGeometryMap.find(neighborIndex) != neighborGeometryMap.end()) {
                    isInSame = true;
                  }
                } // and check, if the neighbor is in the same subdomain
                if (!isInSame) {
                  // the neighbor is not part of this oversampled subdomain,
                  // so the entity in question is on the boundary!
                  const int intersectionLocalIndex = intersection.indexInInside();
                  auto& localBoundaryInfo = *(oversamplingSubdomainInnerBoundaryInfos[entitySubdomain]);
                  // get the boundary info map for this entity
                  IntersectionToBoundaryIdMapType& entityBoundaryInfo = localBoundaryInfo[entityIndex];
                  // and add the local intersection id and its desired fake boundary id to this entities map
                  entityBoundaryInfo.insert(std::pair< int, int >(intersectionLocalIndex, boundaryId_));
                } // if (!isInSame)
              }
            } // then walk the neighbors
          }
        }
      } // now we find all the oversampled subdomains this entity is a part of
    } // walk the global grid part

    // and create the oversampled local grid parts
    std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > >
      oversampledLocalGridPartsRet = std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > >(
          new std::vector< std::shared_ptr< const LocalGridPartType > >(size_));
    std::vector< std::shared_ptr< const LocalGridPartType > >& oversampledLocalGridParts = *oversampledLocalGridPartsRet;
    for (typename SubdomainMapType::const_iterator subdomainIterator = subdomainToOversamplingEntitiesMap.begin();
         subdomainIterator != subdomainToOversamplingEntitiesMap.end();
         ++subdomainIterator) {
      // report
      const size_t subdomain = subdomainIterator->first;
      // for the local grid part
      //   * get the geometry map
      const std::shared_ptr< const GeometryMapType > localGeometryMap = subdomainIterator->second;
      //   * get the boundary info map
      const std::shared_ptr< const EntityToIntersectionInfoMapType > localBoundaryInfo = oversamplingSubdomainInnerBoundaryInfos[subdomain];
      //   * and create the local grid part
      oversampledLocalGridParts[subdomain] = std::shared_ptr< const LocalGridPartType >(
            new LocalGridPartType(globalGridPart_,
                                  localGeometryMap,
                                  localBoundaryInfo));
    } // and crete the oversampled local grid parts
    return oversampledLocalGridPartsRet;
  } // ... addOneLayerOfOverSampling(...)

  template< class EntityType, class NeighborType >
  void add_neighbors_neighbors_recursively(const EntityType& entity,
                                           const NeighborType& neighbor,
                                           const GeometryMapType& geometryMap,
                                           size_t recursion_level,
                                           CodimSizesType& localCodimSizes,
                                           GeometryMapType& geometryMapCopy)
  {
    // loop over all the neighbors of the neighbor
    for (auto neighborIntersectionIt = globalGridPart_->ibegin(neighbor);
         neighborIntersectionIt != globalGridPart_->iend(neighbor);
         ++neighborIntersectionIt) {
      const auto& neighborIntersection = *neighborIntersectionIt;
      if (neighborIntersection.neighbor()) {
        // get the neighbors neighbor
        const auto neighborsNeighborPtr = neighborIntersection.outside();
        const auto& neighborsNeighbor = *neighborsNeighborPtr;
        const IndexType neighborsNeighborGlobalIndex = globalGridPart_->indexSet().index(neighborsNeighbor);
        bool neighborsNeighborIsNotInThisSubdomain = true;
        if (geometryMap.find(neighborsNeighbor.type()) != geometryMap.end()) {
          if (geometryMap.find(neighborsNeighbor.type())->second.find(neighborsNeighborGlobalIndex)
              != geometryMap.find(neighborsNeighbor.type())->second.end()) {
            neighborsNeighborIsNotInThisSubdomain = false;
          }
        }
        // if the neighbor is not in the subdomain
        if (neighborsNeighborIsNotInThisSubdomain) {
          // check, if he intersects the entity
          const auto& entityGeometry = entity.geometry();
          const auto& neighborsNeighborGeometry = neighborsNeighbor.geometry();
          // * therefore loop over all corners of the entity
          for (int ii = 0; ii < entityGeometry.corners(); ++ii) {
            const auto entityCorner = entityGeometry.corner(ii);
            // then loop over all the corners of the neighbors neighbor
            for (int jj = 0; jj < neighborsNeighborGeometry.corners(); ++jj) {
              const auto neighborsNeighborCorner = neighborsNeighborGeometry.corner(jj);
              // and check for equality
              if (entityCorner == neighborsNeighborCorner) {
                // then add the neighbors neighbor
                addGeometryAndIndex(geometryMapCopy,
                                    localCodimSizes,
                                    neighborsNeighbor.type(),
                                    neighborsNeighborGlobalIndex);
                Add< 1, dim >::subEntities(*this, neighborsNeighbor, geometryMapCopy, localCodimSizes);
              }
            }
          }
        } // if the neighbor is not in the subdomain
        // call this function on the neighbours neighbor
        if (recursion_level > 0)
          add_neighbors_neighbors_recursively(entity,
                                              neighborsNeighbor,
                                              geometryMap,
                                              --recursion_level,
                                              localCodimSizes,
                                              geometryMapCopy);
      }
    } // loop over all the neighbors of the neighbor
  } // ... add_neighbors_neighbors_recursively(...)

  // friends
  template< int, int >
  friend struct Add;

  // members
  const std::shared_ptr< const GridType > grid_;
  const int boundaryId_;
  bool prepared_;
  bool finalized_;
  size_t size_;
  std::shared_ptr< const GlobalGridPartType > globalGridPart_;
  // for the entity <-> subdomain relations
  std::shared_ptr< EntityToSubdomainMapType > entityToSubdomainMap_;
  SubdomainMapType subdomainToEntityMap_;
  SubdomainMapType subdomainToOversamplingEntitiesMap_;
  // for the neighboring information
  std::shared_ptr< std::vector< NeighboringSubdomainsSetType > > neighboringSubdomainSets_;
  // for the local grid parts
  std::map< size_t, CodimSizesType > localCodimSizes_;
  std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > > localGridParts_;
  std::shared_ptr< std::vector< std::shared_ptr< const LocalGridPartType > > > oversampledLocalGridParts_;
  // for the boundary grid parts
  std::shared_ptr< std::map< size_t, std::shared_ptr< const BoundaryGridPartType > > > boundaryGridParts_;
  // for the coupling grid parts
  std::shared_ptr< std::vector< std::map< size_t, std::shared_ptr< const CouplingGridPartType > > > > couplingGridPartsMaps_;
  bool oversampled_;
}; // class Default

//! specialization to stop the recursion
template< class GridType >
template< int c >
struct Default< GridType >::Add< c, c >
{
  static void subEntities(Default< GridType >& factory,
                          const typename Default< GridType >::EntityType& entity,
                          typename Default< GridType >::GeometryMapType& geometryMap,
                          typename Default< GridType >::CodimSizesType& localCodimSizes,
                          const std::string prefix = "",
                          std::ostream& out = Dune::Stuff::Common::Logger().devnull())
  {
//    // suppress output, since we are not codim 0
//    Dune::Stuff::Common::LogStream& devnull = Dune::Stuff::Common::Logger().devnull();
    // loop over all codim c subentities of this entity
    typedef typename Default< GridType >::EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
    for (int i = 0; i < entity.template count< c >(); ++i) {
      const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
      const Default< GridType >::GeometryType& geometryType = codimCentityPtr->type();
      const typename Default< GridType >::IndexType globalIndex = factory.globalGridPart_->indexSet().index(*codimCentityPtr);
      factory.addGeometryAndIndex(geometryMap, localCodimSizes, geometryType, globalIndex, prefix, out);
    } // loop over all codim c subentities of this entity
  } // static void subEntities()
}; // struct Default< GridType >::Add< c, c >

} // namespace Factory
} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_FACTORY_DEFAULT_HH
