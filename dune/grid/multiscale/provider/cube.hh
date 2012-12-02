#ifndef DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/sgrid.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/parameter/tree.hh>

#include <dune/grid/multiscale/factory/default.hh>

namespace Dune {
namespace grid {
namespace Multiscale {
namespace Provider {

#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridImp = Dune::GridSelector::GridType >
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template< class GridImp = Dune::SGrid< 2, 2 > >
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class Cube
  : public Dune::Stuff::Grid::Provider::Cube< GridImp >
{
public:
  typedef GridImp GridType;

  typedef Dune::Stuff::Grid::Provider::Cube< GridType > BaseType;

  typedef Cube< GridType > ThisType;

  static const unsigned int dim = BaseType::dim;

  typedef typename BaseType::ctype ctype;
  typedef typename BaseType::CoordinateType CoordinateType;

private:
  typedef Dune::grid::Multiscale::Factory::Default< GridType > MsGridFactoryType;

public:
  typedef typename MsGridFactoryType::MsGridType MsGridType;

private:
  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

public:
  static const std::string id()
  {
    return "grid.multiscale.provider.cube";
  }

  Cube(const double lowerLeft = 0.0,
       const double upperRight = 1.0,
       const unsigned int numElements = 4u,
       const unsigned int partitions = 2u)
    : BaseType(lowerLeft, upperRight, numElements)
  {
    setup(std::vector< unsigned int >(dim, partitions));
  }

  Cube(const CoordinateType& lowerLeft,
       const CoordinateType& upperRight,
       const unsigned int numElements = 4u,
       const unsigned int partitions = 2u)
    : BaseType(lowerLeft, upperRight, numElements)
  {
    setup(std::vector< unsigned int >(dim, partitions));
  }

  template< class ElementsContainerType, class PartitionsContainerType >
  Cube(const CoordinateType& lowerLeft,
       const CoordinateType& upperRight,
       const ElementsContainerType numElements
        = boost::assign::list_of< typename ElementsContainerType::value_type>().repeat(dim,
                                                                                       typename ElementsContainerType::value_type(4u)),
       const PartitionsContainerType partitions
        = boost::assign::list_of< typename PartitionsContainerType::value_type>().repeat(dim,
                                                                                         typename PartitionsContainerType::value_type(2u)))
    : BaseType(lowerLeft, upperRight, numElements)
  {
    std::vector< unsigned int > tmpPartitions(dim);
    static_assert(std::is_unsigned< typename PartitionsContainerType::value_type >::value
                  && std::is_integral< typename PartitionsContainerType::value_type >::value,
                  "Only unsigned integral number of partitions per dimension allowed!");
    std::copy(partitions.begin(), partitions.end(), tmpPartitions.begin());
    setup(tmpPartitions);
  }

  Cube(const ThisType& other)
    : BaseType(other)
    , msGrid_(other.msGrid_)
  {}

  //! \todo TODO Use BaseType::createFromParamTree() internally to avoid duplication!
  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree extendedParamTree;
    if (paramTree.hasSub(subName))
      extendedParamTree = paramTree.sub(subName);
    else
      extendedParamTree = paramTree;
    // get lower left
    std::vector< ctype > lowerLefts;
    if (extendedParamTree.hasVector("lowerLeft")) {
      lowerLefts = extendedParamTree.getVector("lowerLeft", ctype(0));
      assert(lowerLefts.size() >= dim && "Given vector too short!");
    } else if (extendedParamTree.hasKey("lowerLeft")) {
        const ctype lowerLeft = extendedParamTree.get("lowerLeft", ctype(0));
        lowerLefts = std::vector< ctype >(dim, lowerLeft);
    } else {
      std::cout << "Warning in " << id() << ": neither vector nor key 'lowerLeft' given, defaulting to 0.0!" << std::endl;
      lowerLefts = std::vector< ctype >(dim, ctype(0));
    }
    // get upper right
    std::vector< ctype > upperRights;
    if (extendedParamTree.hasVector("upperRight")) {
      upperRights = extendedParamTree.getVector("upperRight", ctype(1));
      assert(upperRights.size() >= dim && "Given vector too short!");
    } else if (extendedParamTree.hasKey("upperRight")) {
        const ctype upperRight = extendedParamTree.get("upperRight", ctype(1));
        upperRights = std::vector< ctype >(dim, upperRight);
    } else {
      std::cout << "Warning in " << id() << ": neither vector nor key 'upperRight' given, defaulting to 1.0!" << std::endl;
      upperRights = std::vector< ctype >(dim, ctype(1));
    }
    // get number of elements
    std::vector< unsigned int > tmpNumElements;
    if (extendedParamTree.hasVector("numElements")) {
      tmpNumElements = extendedParamTree.getVector("numElements", 4u);
      assert(tmpNumElements.size() >= dim && "Given vector too short!");
    } else if (extendedParamTree.hasKey("numElements")) {
        const unsigned int numElement = extendedParamTree.get("numElements", 4u);
        tmpNumElements = std::vector< unsigned int >(dim, numElement);
    } else {
      std::cout << "Warning in " << id() << ": neither vector nor key 'numElements' given, defaulting to 4!" << std::endl;
      tmpNumElements = std::vector< unsigned int >(dim, 4u);
    }
    // get partitions
    std::vector< unsigned int > tmpPartitions;
    if (extendedParamTree.hasVector("partitions")) {
      tmpPartitions = extendedParamTree.getVector("partitions", 2u);
      assert(tmpPartitions.size() >= dim && "Given vector too short!");
    } else if (extendedParamTree.hasKey("partitions")) {
        const unsigned int partitions = extendedParamTree.get("partitions", 2u);
        tmpPartitions = std::vector< unsigned int >(dim, partitions);
    } else {
      std::cout << "Warning in " << id() << ": neither vector nor key 'partitions' given, defaulting to 2!" << std::endl;
      tmpPartitions = std::vector< unsigned int >(dim, 2u);
    }
    // check and save
    CoordinateType lowerLeft;
    CoordinateType upperRight;
    Dune::array< unsigned int, dim > numElements;
    for (unsigned int d = 0; d < dim; ++d) {
      assert(lowerLefts[d] < upperRights[d]
             && "Given 'upperRight' hast to be elementwise larger than given 'lowerLeft'!");
      lowerLeft[d] = lowerLefts[d];
      upperRight[d] = upperRights[d];
      assert(tmpNumElements[d] > 0 && "Given 'numElements' has to be elementwise positive!");
      numElements[d] = tmpNumElements[d];
      assert(tmpPartitions[d] > 0 && "Given 'partitions' has to be elementwise positive!");
    }
    return Cube(lowerLeft, upperRight, numElements, tmpPartitions);
  } // static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())

  const ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      BaseType::operator=(other);
      msGrid_ = other.msGrid();
    }
    return this;
  } // ThisType& operator=(ThisType& other)

  const Dune::shared_ptr< const MsGridType > msGrid() const
  {
    return msGrid_;
  }

private:
  void setup(std::vector< unsigned int >& partitions)
  {
    // prepare
    const CoordinateType& lowerLeft = BaseType::lowerLeft();
    const CoordinateType& upperRight = BaseType::upperRight();
    MsGridFactoryType factory(BaseType::gridPtr());
    factory.prepare();

    // debug output
    const std::string prefix = "";
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    debug << prefix << id << ":" << std::endl;
    debug << prefix << "  lowerLeft: " << lowerLeft << std::endl;
    debug << prefix << "  upperRight: " << upperRight << std::endl;
    debug << prefix << "  partitions per dim: ";
    for (unsigned int d = 0; d < dim; ++d)
      debug << partitions[d] << " ";
    debug << std::endl;

    // global grid part
    typedef typename MsGridType::GlobalGridPartType GridPartType;
    const Dune::shared_ptr< const GridPartType> gridPart = factory.globalGridPart();

    // walk the grid
    for (typename GridPartType::template Codim< 0 >::IteratorType it = gridPart->template begin< 0 >();
        it != gridPart->template end< 0 >();
        ++it) {
      // get center of entity
      typename GridType::LeafGridView::template Codim< 0 >::Iterator::Entity& entity = *it;
      const CoordinateType center = entity.geometry().global(entity.geometry().center());
//      debug << prefix << "  entity (" << center << "):" << std::endl;

      // decide on the subdomain this entity shall belong to
      std::vector< unsigned int > whichPartition;
      for (unsigned int d = 0; d < dim; ++d) {
        whichPartition.push_back(std::floor(partitions[d]*((center[d] - lowerLeft[d])/(upperRight[d] - lowerLeft[d]))));
      }
      unsigned int subdomain = 0;
      if (dim == 1)
        subdomain = whichPartition[0];
      else if (dim == 2)
        subdomain = whichPartition[0] + whichPartition[1]*partitions[0];
      else if (dim == 3)
        subdomain = whichPartition[0] + whichPartition[1]*partitions[0] + whichPartition[2]*partitions[1]*partitions[0];
      else {
        std::stringstream msg;
        msg << "Error in " << id << ": not implemented for grid dimensions other than 1, 2 or 3!";
        DUNE_THROW(Dune::NotImplemented, msg.str());
      } // decide on the subdomain this entity shall belong to

      // add entity to subdomain
      factory.add(entity, subdomain, prefix + "  ", debug);
    } // walk the grid

    // finalize
    factory.finalize(prefix + "  ", debug);
    debug << std::flush;

    msGrid_ = factory.createMsGrid();
  } // void setup(const Dune::ParameterTree& paramTree)

  Dune::shared_ptr< const MsGridType > msGrid_;
}; // class Cube

} // namespace Provider
} // namespace Multiscale
} // namespace Grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
