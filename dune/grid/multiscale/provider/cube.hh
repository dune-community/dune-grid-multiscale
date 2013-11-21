// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH

#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/sgrid.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/print.hh>

#include <dune/grid/multiscale/factory/default.hh>

#include "interface.hh"

namespace Dune {
namespace grid {
namespace Multiscale {


template< class GridImp >
class ProviderCube
  : public ProviderInterface< GridImp >
{
  typedef Dune::Stuff::GridProviderCube< GridImp > CubeGridProvider;
  typedef ProviderInterface< GridImp > InterfaceType;
public:
  static const unsigned int dimension = InterfaceType::dimension;
  static const unsigned int dim = InterfaceType::dimension;
  typedef typename InterfaceType::ctype           ctype;
  typedef typename InterfaceType::CoordinateType  CoordinateType;
  typedef GridImp GridType;
private:
  typedef Dune::grid::Multiscale::Factory::Default< GridType > MsGridFactoryType;
public:
  typedef typename MsGridFactoryType::MsGridType MsGridType;
private:
  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

public:
  static std::string static_id()
  {
    return InterfaceType::static_id() + ".cube";
  }

  virtual std::string id() const
  {
    return InterfaceType::static_id() + ".cube";
  }

  ProviderCube(const std::vector< ctype > lower_left       = std::vector< ctype >({0.0, 0.0, 0.0}),
               const std::vector< ctype > upper_right      = std::vector< ctype >({1.0, 1.0, 1.0}),
               const std::vector< size_t > num_elements    = std::vector< size_t >({8, 8, 8}),
               const std::vector< size_t > num_partittions = std::vector< size_t >({2, 2, 2}),
               const size_t num_oversampling_layers        = 0,
               std::ostream& out = DSC_LOG.devnull(), const std::string prefix = "")
  {
    if (lower_left.size() < dim)
      DUNE_THROW(Dune::RangeError,
                 "lower_left has to be at least of size " << dim << " (is " << lower_left.size() << ")!");
    if (upper_right.size() < dim)
      DUNE_THROW(Dune::RangeError,
                 "upper_right has to be at least of size " << dim << " (is " << upper_right.size() << ")!");
    if (num_elements.size() < dim)
      DUNE_THROW(Dune::RangeError,
                 "num_elements has to be at least of size " << dim << " (is " << num_elements.size() << ")!");
    if (num_partittions.size() < dim)
      DUNE_THROW(Dune::RangeError,
                 "num_partittions has to be at least of size " << dim << " (is " << num_partittions.size() << ")!");
    for (size_t ii = 0; ii < dim; ++ii) {
      if (num_partittions[ii] < num_elements[ii])
        DUNE_THROW(Dune::RangeError,
                   num_partittions[ii] << " = num_partittions[" << ii << "] has to be smaller than num_elements[" << ii
                   << "] = " << num_elements[ii] << "!)");
    }
    grid_ = CubeGridProvider(lower_left, upper_right, num_elements).grid();
    setup(lower_left, upper_right, num_partittions, num_oversampling_layers, out, prefix);
  }

  ProviderCube(const std::shared_ptr< const GridType > grd,
               const std::vector< ctype > lower_left       = std::vector< ctype >({0.0, 0.0, 0.0}),
               const std::vector< ctype > upper_right      = std::vector< ctype >({1.0, 1.0, 1.0}),
               const std::vector< size_t > num_partittions = std::vector< size_t >({2, 2, 2}),
               const size_t num_oversampling_layers        = 0,
               std::ostream& out = DSC_LOG.devnull(), const std::string prefix = "")
    : grid_(grd)
  {
    if (lower_left.size() < dim)
      DUNE_THROW(Dune::RangeError,
                 "lower_left has to be at least of size " << dim << " (is " << lower_left.size() << ")!");
    if (upper_right.size() < dim)
      DUNE_THROW(Dune::RangeError,
                 "upper_right has to be at least of size " << dim << " (is " << upper_right.size() << ")!");
    if (num_partittions.size() < dim)
      DUNE_THROW(Dune::RangeError,
                 "num_partittions has to be at least of size " << dim << " (is " << num_partittions.size() << ")!");
    for (size_t ii = 0; ii < dim; ++ii) {
      if (lower_left[ii] >= upper_right[ii])
        DUNE_THROW(Dune::RangeError,
                   lower_left[ii] << " = lower_left[" << ii << "] has to be smaller than upper_right[" << ii
                   << "] = " << upper_right[ii] << "!)");
    }
    setup(lower_left, upper_right, num_partittions, num_oversampling_layers, out, prefix);
  }

//  ProviderCube(const ThisType& other)
//    : BaseType(other)
//    , msGrid_(other.msGrid_)
//  {}

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

//  //! \todo TODO Use BaseType::createFromParamTree() internally to avoid code duplication!
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
//    Dune::array< unsigned int, dim > numElements;
//    for (unsigned int d = 0; d < dim; ++d) {
//      assert(lowerLefts[d] < upperRights[d]
//             && "Given 'upperRight' hast to be elementwise larger than given 'lowerLeft'!");
//      lowerLeft[d] = lowerLefts[d];
//      upperRight[d] = upperRights[d];
//      assert(tmpNumElements[d] > 0 && "Given 'numElements' has to be elementwise positive!");
//      numElements[d] = tmpNumElements[d];
//      assert(tmpPartitions[d] > 0 && "Given 'partitions' has to be elementwise positive!");
//    }
//    return new ThisType(lowerLeft, upperRight, numElements, tmpPartitions, oversamplingLayers);
//  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

//  ThisType& operator=(const ThisType& other)
//  {
//    if (this != &other) {
//      BaseType::operator=(other);
//      msGrid_ = other.msGrid();
//    }
//    return this;
//  } // ThisType& operator=(ThisType& other)

  std::shared_ptr< const GridType > grid() const
  {
    return grid_;
  }

  std::shared_ptr< const MsGridType > msGrid() const
  {
    return ms_grid_;
  }

  using InterfaceType::visualize;

private:
  void setup(const std::vector< ctype >& lower_left,
             const std::vector< ctype >& upper_right,
             const std::vector< size_t >& num_partitions,
             const size_t num_oversampling_layers,
             std::ostream& out = DSC_LOG.devnull(), const std::string prefix = "")
  {
    const size_t neighbor_recursion_level = Factory::NeighborRecursionLevel< GridType >::compute();
    // prepare
    MsGridFactoryType factory(grid_);
    factory.prepare();
#ifndef NDEBUG
    // debug output
    out << prefix << static_id()<< ":" << std::endl;
    Stuff::Common::print(lower_left, "lower_left", out, prefix);
    Stuff::Common::print(upper_right, "upper_right", out, prefix);
    Stuff::Common::print(num_partitions, "num_partitions", out, prefix);
#endif // NDEBUG
    // global grid part
    typedef typename MsGridType::GlobalGridPartType GridPartType;
    const auto global_grid_part = factory.globalGridPart();
    // walk the grid
    const auto entity_it_end = global_grid_part->template end< 0 >();
    for (auto entity_it = global_grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      // get center of entity
      const auto& entity = *entity_it;
      const auto center = entity.geometry().center();
#ifndef NDEBUG
      const size_t entity_index = global_grid_part->indexSet().index(entity);
      Stuff::Common::print(center, "entity (" + Stuff::Common::toString(entity_index) + ")", out, prefix);
#endif // NDEBUG
      // decide on the subdomain this entity shall belong to
      std::vector< size_t > whichPartition(dim, 0);
      for (size_t dd = 0; dd < dim; ++dd)
        whichPartition[dd] = (std::min((size_t)(std::floor(num_partitions[dd]*((center[dd] - lower_left[dd])/(upper_right[dd] - lower_left[dd])))),
                                       num_partitions[dd] - 1));
      size_t subdomain = 0;
      if (dim == 1)
        subdomain = whichPartition[0];
      else if (dim == 2)
        subdomain = whichPartition[0] + whichPartition[1]*num_partitions[0];
      else if (dim == 3)
        subdomain = whichPartition[0] + whichPartition[1]*num_partitions[0] + whichPartition[2]*num_partitions[1]*num_partitions[0];
      else
        DUNE_THROW(Dune::NotImplemented,
                   "ERROR in " << static_id()<< ": not implemented for grid dimensions other than 1, 2 or 3!");
      // add entity to subdomain
      factory.add(entity, subdomain, prefix + "  ", out);
    } // walk the grid
    // finalize
    factory.finalize(num_oversampling_layers, neighbor_recursion_level, prefix + "  ", out);
//    debug << std::flush;
    // be done with it
    ms_grid_ = factory.createMsGrid();
  } // void setup(const Dune::ParameterTree& paramTree)

  std::shared_ptr< const GridType > grid_;
  std::shared_ptr< const MsGridType > ms_grid_;
}; // class ProviderCube


} // namespace Multiscale
} // namespace Grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
