
#ifndef DUNE_GRID_MULTISCALE_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_DEFAULT_HH

// system
#include <vector>
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

// dune-geometry
#include <dune/geometry/type.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/gridpart/leaf.hh>
//#include <dune/grid/multiscale/gridpart/local.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

template< class GridImp >
class Default
{
public:
  typedef GridImp GridType;

  typedef Default< GridType > ThisType;

  static const std::string id;

  static const unsigned int dim = GridType::dimension;

  typedef Dune::grid::Multiscale::GridPart::Leaf< GridType > GlobalGridPartType;

//  typedef Dune::grid::Multiscale::GridPart::Local< GridType > LocalGridPartType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

private:
  typedef typename GlobalGridPartType::IndexSetType::IndexType IndexType;

  typedef Dune::GeometryType GeometryType;

  typedef std::map< IndexType, IndexType > IndexMapType;

  typedef std::map< GeometryType, IndexMapType > GeometryMapType;

  typedef std::map< unsigned int, Dune::shared_ptr< GeometryMapType > > SubdomainMapType;

//  typedef std::vector< Dune::shared_ptr< LocalGridPartType > > LocalGridPartVectorType;

  template< int c, int d >
  struct Add
  {
    static void subEntities(ThisType& msGrid, const EntityType& entity, GeometryMapType& geometryMap, std::ostream& out, const std::string prefix)
    {
      typedef typename EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
      for (int i = 0; i < entity.template count< c >(); ++i)
      {
        const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
        const GeometryType& geometryType = codimCentityPtr->type();
        const IndexType globalIndex = msGrid.globalGridPart_->indexSet().index(*codimCentityPtr);
        msGrid.addGeometryAndIndex(geometryMap, geometryType, globalIndex, out, prefix);
      }
      Add< c + 1, d >::subEntities(msGrid, entity, geometryMap, out, prefix);
    }
  };

public:
  Default(GridType& grid)
    : grid_(grid),
      globalGridPart_(Dune::shared_ptr< GlobalGridPartType >(new GlobalGridPartType(grid_))),
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
    debug << prefix << "adding entity " << globalIndex << " to subdomain " << subdomain << ":" << std::endl;

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
    Add< 1, dim >::subEntities(*this, entity, geometryMap, debug, prefix);
  } // void add(const EntityType& entity, const unsigned int subdomain, Dune::ParameterTree paramTree = Dune::ParameterTree())

#if 0
  void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    //// debug output
    //const std::string prefix = paramTree.get("prefix", "");
    //Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    //debug << prefix << id << ".finalize: " << std::flush;

    // test for consecutive numbering of subdomains and finalize subgrids
    bool consecutive = true;
    for (unsigned int subdomain = 0; subdomain < size(); ++subdomain) {
      if (subdomainToIndexPairMap_.find(subdomain) == subdomainToIndexPairMap_.end())
        consecutive = false;
    }
    if (!consecutive) {
      std::stringstream msg;
      msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }

    // create local grid parts and report
    //debug << "found " << size() << " subdomains" << std::endl;
    for (typename SubdomainToIndexPairMapType::iterator pairOfSubdomainAndIndexPair = subdomainToIndexPairMap_.begin();
         pairOfSubdomainAndIndexPair != subdomainToIndexPairMap_.end();
         ++pairOfSubdomainAndIndexPair) {
//      const unsigned int subdomain = pairOfSubdomainAndIndexPair->first;
//      const unsigned int subdomainSize = pairOfSubdomainAndIndexPair->second->size();
      //debug << prefix << "- subdomain " << subdomain << ", " << subdomainSize << " entities:" << std::endl;
      //debug << prefix << "  " << std::flush;
      unsigned int counter = 0;
      for (typename GlobalToLocalIndexMapType::iterator indexPair = pairOfSubdomainAndIndexPair->second->begin();
           indexPair != pairOfSubdomainAndIndexPair->second->end();
           ++indexPair) {
        //debug << indexPair->first;
        //if (counter < subdomainSize - 1)
          //debug << ", ";
        ++counter;
      }

      //debug << std::endl;
      localGridParts_.push_back(Dune::shared_ptr< LocalGridPartType >(new LocalGridPartType(*globalGridPart_, pairOfSubdomainAndIndexPair->second)));
    } // create local grid parts and report

    // finalize
    finalized_ = true;
    return;
  } // void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())
#endif

//  const Dune::shared_ptr< const LocalGridPartType > localGridPart(const unsigned int subdomain) const
//  {
//    assert(finalized_);
//    assert(subdomain < size());
//    return localGridParts_[subdomain];
//  }

private:
  template< int, int >
  friend class Add;

  void addGeometryAndIndex(GeometryMapType& geometryMap,
                           const GeometryType& geometryType,
                           const IndexType& globalIndex,
                           std::ostream& out,
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
  }

  GridType& grid_;
  Dune::shared_ptr< GlobalGridPartType > globalGridPart_;
  bool finalized_;
  unsigned int size_;
  SubdomainMapType subdomainMap_;
//  LocalGridPartVectorType localGridParts_;
}; // class Default

template< class GridType >
const std::string Default< GridType >::id = "grid.multiscale.default";

template< class GridType >
template< int c >
struct Default< GridType >::Add< c, c >
{
  static void subEntities(Default< GridType >& msGrid,
                          const typename Default< GridType >::EntityType& entity,
                          typename Default< GridType >::GeometryMapType& geometryMap,
                          std::ostream& out,
                          const std::string prefix)
  {
    typedef typename Default< GridType >::EntityType::template Codim< c >::EntityPointer CodimCentityPtrType;
    for (int i = 0; i < entity.template count< c >(); ++i)
    {
      const CodimCentityPtrType codimCentityPtr = entity.template subEntity< c >(i);
      const Default< GridType >::GeometryType& geometryType = codimCentityPtr->type();
      const typename Default< GridType >::IndexType globalIndex = msGrid.globalGridPart_->indexSet().index(*codimCentityPtr);
      msGrid.addGeometryAndIndex(geometryMap, geometryType, globalIndex, out, prefix);
    }
  }
};

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_DEFAULT_HH
