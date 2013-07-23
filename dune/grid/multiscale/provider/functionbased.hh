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

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/functions/fromfile.hh>

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
  typedef int RangeFieldType;
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

//  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
//  {
//    Dune::ParameterTree description;
//    description["lowerLeft"] = "[0.0; 0.0; 0.0]";
//    description["upperRight"] = "[1.0; 1.0; 1.0]";
//    description["numElements"] = "[4; 4; 4]";
//    description["partitions"] = "[2; 2; 2]";
//    description["oversamplingLayers"] = "2";
//    if (subName.empty())
//      return description;
//    else {
//      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
//      extendedDescription.add(description, subName);
//      return extendedDescription;
//    }
//  } // ... createSampleDescription(...)


  //! \todo TODO Use NonMultiscaleType::createFromParamTree() internally to avoid code duplication!
//  static ThisType* create(const Dune::ParameterTree& description, const std::string subName = id())
//  {
//    // get correct paramTree
//    Dune::Stuff::Common::ExtendedParameterTree extendedParamTree;
//    if (description.hasSub(subName))
//      extendedParamTree = description.sub(subName);
//    else
//      extendedParamTree = description;
//    // get lower left
//    std::vector< ctype > lowerLefts;
//    if (extendedParamTree.hasVector("lowerLeft")) {
//      lowerLefts = extendedParamTree.getVector("lowerLeft", ctype(0), dim);
//      assert(lowerLefts.size() >= dim && "Given vector too short!");
//    } else if (extendedParamTree.hasKey("lowerLeft")) {
//        const ctype lowerLeft = extendedParamTree.get("lowerLeft", ctype(0));
//        lowerLefts = std::vector< ctype >(dim, lowerLeft);
//    } else {
//      std::cout << "WARNING in " << id() << ": neither vector nor key 'lowerLeft' given, defaulting to 0.0!" << std::endl;
//      lowerLefts = std::vector< ctype >(dim, ctype(0));
//    }
//    // get upper right
//    std::vector< ctype > upperRights;
//    if (extendedParamTree.hasVector("upperRight")) {
//      upperRights = extendedParamTree.getVector("upperRight", ctype(1), dim);
//      assert(upperRights.size() >= dim && "Given vector too short!");
//    } else if (extendedParamTree.hasKey("upperRight")) {
//        const ctype upperRight = extendedParamTree.get("upperRight", ctype(1));
//        upperRights = std::vector< ctype >(dim, upperRight);
//    } else {
//      std::cout << "WARNING in " << id() << ": neither vector nor key 'upperRight' given, defaulting to 1.0!" << std::endl;
//      upperRights = std::vector< ctype >(dim, ctype(1));
//    }
//    // get number of elements
//    std::vector< unsigned int > tmpNumElements;
//    if (extendedParamTree.hasVector("numElements")) {
//      tmpNumElements = extendedParamTree.getVector("numElements", 4u, dim);
//      assert(tmpNumElements.size() >= dim && "Given vector too short!");
//    } else if (extendedParamTree.hasKey("numElements")) {
//        const unsigned int numElement = extendedParamTree.get("numElements", 4u);
//        tmpNumElements = std::vector< unsigned int >(dim, numElement);
//    } else {
//      std::cout << "WARNING in " << id() << ": neither vector nor key 'numElements' given, defaulting to 4!" << std::endl;
//      tmpNumElements = std::vector< unsigned int >(dim, 4u);
//    }
//    // get partitions
//    std::vector< unsigned int > tmpPartitions;
//    if (extendedParamTree.hasVector("partitions")) {
//      tmpPartitions = extendedParamTree.getVector("partitions", 2u, dim);
//      assert(tmpPartitions.size() >= dim && "Given vector too short!");
//    } else if (extendedParamTree.hasKey("partitions")) {
//        const unsigned int partitions = extendedParamTree.get("partitions", 2u);
//        tmpPartitions = std::vector< unsigned int >(dim, partitions);
//    } else {
//      std::cout << "WARNING in " << id() << ": neither vector nor key 'partitions' given, defaulting to 2!" << std::endl;
//      tmpPartitions = std::vector< unsigned int >(dim, 2u);
//    }
//    // get oversampling size
//    const size_t oversamplingLayers = extendedParamTree.get< size_t >("oversamplingLayers", 0);
//    // check and save
//    CoordinateType lowerLeft;
//    CoordinateType upperRight;
//    std::vector< unsigned int> numElements;
//    for (unsigned int d = 0; d < dim; ++d) {
//      assert(lowerLefts[d] < upperRights[d]
//             && "Given 'upperRight' hast to be elementwise larger than given 'lowerLeft'!");
//      lowerLeft[d] = lowerLefts[d];
//      upperRight[d] = upperRights[d];
//      assert(tmpNumElements[d] > 0 && "Given 'numElements' has to be elementwise positive!");
//      numElements.push_back(tmpNumElements[d]);
//      assert(tmpPartitions[d] > 0 && "Given 'partitions' has to be elementwise positive!");
//    }
//    // get filename for the function from file
//    const std::string filename = extendedParamTree.get< std::string >("filename", 0);
//    FunctionType function(filename, lowerLeft, upperRight, numElements);
//    return new ThisType(function, /*lowerLeft, upperRight, numElements,*/ tmpPartitions, oversamplingLayers);
//  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

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
