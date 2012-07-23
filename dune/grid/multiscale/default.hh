
#ifndef DUNE_GRID_MULTISCALE_DEFAULT_HH
#define DUNE_GRID_MULTISCALE_DEFAULT_HH

// system
#include <vector>
#include <set>
#include <map>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

// dune-stuff
#include <dune/stuff/common/logging.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/gridpart/leaf.hh>
#include <dune/grid/multiscale/gridpart/local.hh>

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

  typedef Dune::grid::Multiscale::GridPart::Leaf< GridType > GlobalGridPartType;

  typedef Dune::grid::Multiscale::GridPart::Local< GridType > LocalGridPartType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

private:
  typedef typename GlobalGridPartType::IndexSetType GlobalIndexSetType;

  typedef typename GlobalIndexSetType::IndexType GlobalIndexType;

  typedef std::set< GlobalIndexType > LocalIndicesSetType;

  typedef std::map< unsigned int, LocalIndicesSetType > SubdomainMapType;

public:
  Default(GridType& grid)
    : grid_(grid),
      globalGridPart_(Dune::shared_ptr< const GlobalGridPartType >(new GlobalGridPartType(grid_))),
      finalized_(false),
      size_(0)
  {}

  const GlobalGridPartType& globalGridPart() const
  {
    return *(globalGridPart_);
  }

  unsigned int size() const
  {
    return size_;
  }

  void prepare()
  {
    return;
  }

  void add(const EntityType& entity, unsigned int subdomain)
  {
    assert(!finalized_);
    // create subdomain map if needed
    if (subdomainMap_.find(subdomain) == subdomainMap_.end()) {
      subdomainMap_[subdomain] = LocalIndicesSetType();
      ++size_;
    }
    // add index of entity
    subdomainMap_.find(subdomain)->second.insert(globalGridPart_->indexSet().index(entity));
  }

  void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())
  {
    // debug output
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    debug << prefix << id << ".finalize: " << std::flush;

    // test for consecutive numbering of subdomains and finalize subgrids
    bool consecutive = true;
    for (unsigned int subdomain = 0; subdomain < size(); ++subdomain) {
      if (subdomainMap_.find(subdomain) == subdomainMap_.end())
        consecutive = false;
    }
    if (!consecutive) {
      std::stringstream msg;
      msg << "Error in " << id << ": numbering of subdomains has to be consecutive upon calling finalize!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }

    // report
    debug << "found " << size() << " subdomains" << std::endl;
    for (typename SubdomainMapType::iterator subdomain = subdomainMap_.begin(); subdomain != subdomainMap_.end(); ++subdomain) {
      const unsigned int subdomainSize = subdomain->second.size();
      debug << prefix << "- subdomain " << subdomain->first << ", " << subdomainSize << " entities:" << std::endl;
      debug << prefix << "  " << std::flush;
      unsigned int counter = 0;
      for (typename LocalIndicesSetType::iterator localIndex = subdomain->second.begin(); localIndex != subdomain->second.end(); ++localIndex) {
        debug << *localIndex;
        if (counter < subdomainSize - 1)
          debug << ", ";
        ++counter;
      }
      debug << std::endl;
    } // report
    return;
  } // void finalize(Dune::ParameterTree paramTree = Dune::ParameterTree())

private:
  GridType& grid_;
  const Dune::shared_ptr< const GlobalGridPartType > globalGridPart_;
  bool finalized_;
  unsigned int size_;
  SubdomainMapType subdomainMap_;
}; // class Default

template< class GridType >
const std::string Default< GridType >::id = "grid.multiscale.default";

} // namespace Multiscale

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_DEFAULT_HH
