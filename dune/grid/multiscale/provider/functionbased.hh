// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_FUNCTIONBASED_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_FUNCTIONBASED_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/sgrid.hh>

#include <dune/stuff/grid/provider.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/functions.hh>

#include <dune/grid/multiscale/factory/default.hh>

#include "interface.hh"

namespace Dune {
namespace grid {
namespace Multiscale {


#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridImp = Dune::GridSelector::GridType >
#else
template< class GridImp = Dune::SGrid< 2, 2 > >
#endif
class ProviderFunctionbased
  : public ProviderInterface< GridImp >
{
public:
  typedef ProviderFunctionbased< GridImp >              ThisType;
  typedef Dune::Stuff::GridProviderInterface< GridImp > NonMultiscaleType;
  typedef ProviderInterface< GridImp >                  InterfaceType;
  typedef GridImp                                       GridType;

  static const unsigned int dim = NonMultiscaleType::dim;

  typedef typename NonMultiscaleType::ctype ctype;
  typedef typename NonMultiscaleType::CoordinateType CoordinateType;
  typedef double RangeFieldType;
  // until now only implemented for scalar functions
  typedef Dune::Stuff::FunctionInterface< ctype, dim, RangeFieldType, 1, 1 > FunctionType;
private:
  typedef Dune::grid::Multiscale::Factory::Default< GridType > MsGridFactoryType;

public:
  typedef typename MsGridFactoryType::MsGridType MsGridType;

private:
  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

public:
  static const std::string id()
  {
    return InterfaceType::id() + ".functionbased";
  }

  ProviderFunctionbased(std::shared_ptr<const GridType > grid,
                        std::shared_ptr<const FunctionType > function,
                        std::vector<double> partitions)
    : grid_(grid)
  {
    // prepare
    MsGridFactoryType factory(grid);
    factory.prepare();
    // debug output
    const std::string prefix = "";
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().devnull();
    // global grid part
    typedef typename MsGridType::GlobalGridPartType GridPartType;
    const Dune::shared_ptr< const GridPartType> gridPart = factory.globalGridPart();
    // walk the grid
    for (typename GridPartType::template Codim< 0 >::IteratorType it = gridPart->template begin< 0 >();
        it != gridPart->template end< 0 >();
        ++it) {
      // get center of entity
      const auto& entity = *it;
      const CoordinateType center = entity.geometry().center();
      debug << prefix << "  entity (" << center << "):" << std::endl;
      // decide on the subdomain this entity shall belong to
      unsigned int subdomain = 0;
      int size = partitions.size();
      for (int ii = 1; ii < size; ii++){
        if ( ((function->evaluate(center)) > partitions[ii-1]) && ((function->evaluate(center)) < partitions[ii] )){
          subdomain = ii;
        }
      }
      if ((function->evaluate(center)) > partitions[size-1]){
        subdomain = size;
      }
      // add entity to subdomain
      factory.add(entity, subdomain, prefix + "  ", debug);

    } // walk the grid
    const unsigned int numberOfSubdomains = partitions.size() + 1;
    debug << prefix << "  the grid has " << numberOfSubdomains << " subdomains." << std::endl;
    // finalize
    factory.finalize(0, prefix + "  ", debug, false);
    //debug << std::flush;

    // be done with it
    msGrid_ = factory.createMsGrid();
  }

  ProviderFunctionbased(const ThisType& other)
    : NonMultiscaleType(other)
    , msGrid_(other.msGrid_)
  {}

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    typedef Stuff::GridProviders< GridType > Providers;
    typedef Stuff::Functions< ctype, dim, RangeFieldType, 1, 1 > Functions;
    Dune::ParameterTree description;
    description["provider"] = Providers::available()[0];
    description["function"] = Functions::available()[0];
    description["partitions"] = "[0.5; 1.0]";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::Configuration extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... createSampleDescription(...)


  static ThisType* create(const Dune::ParameterTree& description, const std::string subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::Configuration extendedParamTree;
    if (description.hasSub(subName))
      extendedParamTree = description.sub(subName);
    else
      extendedParamTree = description;

    // get grid provider
    NonMultiscaleType gridProvider;
    if (extendedParamTree.hasKey("provider")) {
        typedef Stuff::GridProviders< GridType > Providers;
        gridProvider = extendedParamTree.get("provider", Providers::available()[0]);
    } else {
        DUNE_THROW(Dune::IOError, "Error: no key 'provider' given in " << id() << std::endl);
    }
    std::shared_ptr< GridType > grid(gridProvider->grid());
    // get function
    FunctionType function;
    if (extendedParamTree.hasKey("function")) {
        typedef Stuff::Functions< ctype, dim, RangeFieldType, 1, 1 > Functions;
        gridProvider = extendedParamTree.get("function", Functions::available()[0]);
    } else {
        DUNE_THROW(Dune::IOError, "Error: no key 'function' given in " << id() << std::endl);
    }
    // get partitions
    std::vector< RangeFieldType > partitions;
    if (extendedParamTree.hasVector("partitions")) {
      partitions = extendedParamTree.getVector< RangeFieldType >("partitions", 1);
    } else if (extendedParamTree.hasKey("partitions")) {
      partitions = std::vector< RangeFieldType >(1, extendedParamTree.get("partitions", 1u));
    } else {
      std::cout << "WARNING in " << id() << ": neither vector nor key 'partitions' given, defaulting to 1!" << std::endl;
      partitions = std::vector< RangeFieldType >(1, 1u);
    }
    return new ThisType(grid, function, partitions);
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

  const ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      NonMultiscaleType::operator=(other);
      msGrid_ = other.msGrid();
    }
    return this;
  } // ThisType& operator=(ThisType& other)

  virtual const Dune::shared_ptr< const GridType > grid() const
  {
    return grid_;
  }

  virtual const Dune::shared_ptr< const MsGridType > msGrid() const
  {
    return msGrid_;
  }

  using InterfaceType::visualize;

private:
  Dune::shared_ptr< const GridType > grid_;
  Dune::shared_ptr< const MsGridType > msGrid_;
}; // class ProviderFunctionbased


} // namespace Multiscale
} // namespace Grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_FUNCTIONBASED_HH
