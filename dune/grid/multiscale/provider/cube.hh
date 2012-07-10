#ifndef DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH

//#ifdef HAVE_DUNE_STUFF

// system
#include <vector>
#include <sstream>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>

// dune-stuff
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/logging.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/factory/subgrid.hh>

namespace Dune {

namespace grid {

namespace Multiscale {

namespace Provider {

template< class GridImp >
class Cube
  : public Dune::Stuff::Grid::Provider::Cube< GridImp >
{
public:
  typedef Dune::Stuff::Grid::Provider::Cube< GridImp > BaseType;

  typedef typename BaseType::GridType GridType;

  typedef Cube< GridType > ThisType;

  static const std::string id;

  static const unsigned int dim = BaseType::dim;

  typedef typename BaseType::CoordinateType CoordinateType;

private:
  typedef Dune::grid::Multiscale::Factory::Subgrid::FromFineGrid< GridType > MsGridFactoryType;

public:
  typedef typename MsGridFactoryType::MsGridType MsGridType;

  Cube(Dune::ParameterTree paramTree = Dune::ParameterTree())
    : BaseType(paramTree)
  {
    buildMsGrid(paramTree);
  }

  MsGridType& msGrid()
  {
    return *msGrid_;
  }

  const MsGridType& msGrid() const
  {
    return *msGrid_;
  }

private:
  void buildMsGrid(Dune::ParameterTree paramTree)
  {
    // preparations
    const CoordinateType& lowerLeft = BaseType::lowerLeft();
    const CoordinateType& upperRight = BaseType::upperRight();

    // debug output
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::Logging::LogStream& debug = Dune::Stuff::Common::Logger().Dbg();
    debug << prefix << id << ".buildMsGrid():" << std::endl;
    debug << prefix << "  lowerLeft: " << lowerLeft << std::endl;
    debug << prefix << "  upperRight: " << upperRight << std::endl;
    debug << prefix << "  partitions per dim: ";

    // get number of subdomains per direction
    std::vector< unsigned int > partitions;
    for (unsigned int d = 0; d < dim; ++d) {
      std::stringstream key;
      key << "partitions." << d;
      partitions.push_back(paramTree.get(key.str(), 1));
      debug << partitions[d] << " ";
    }
    debug << std::endl;


    // create factory
    MsGridFactoryType msGridFactory(BaseType::grid());
    msGridFactory.prepare();

    // walk the host grid
    typename GridType::LeafGridView::template Codim< 0 >::Iterator itEnd = BaseType::grid().leafView().template end< 0 >();
    for (typename GridType::LeafGridView::template Codim< 0 >::Iterator it = BaseType::grid().leafView().template begin< 0 >();
         it != itEnd;
         ++it) {
      // get center of entity
      typename GridType::LeafGridView::template Codim< 0 >::Iterator::Entity& entity = *it;
      const CoordinateType center = entity.geometry().global(entity.geometry().center());
      debug << prefix << "  - entity: " << center << std::endl;

      // decide on the subdomain this entity shall belong to
      std::vector< unsigned int > whichPartition;
      for (unsigned int d = 0; d < dim; ++d) {
          whichPartition.push_back(std::floor(partitions[d]*((center[d] - lowerLeft[d])/(upperRight[d] - lowerLeft[d]))));
          debug << prefix << "    partition[" << d  << "]: " << whichPartition[d] << std::endl;
      }
      unsigned int subdomain = 0;
      if (dim == 1)
        subdomain = whichPartition[0];
      else if (dim == 2)
        subdomain = whichPartition[0] + whichPartition[1]*partitions[0];
      else if (dim == 3)
        subdomain = whichPartition[0] + whichPartition[1]*partitions[0] + whichPartition[2]*partitions[1];
      else {
        std::stringstream msg;
        msg << "Error in " << id << ": not implemented for grid dimensions other than 1, 2 or 3!";
        DUNE_THROW(Dune::InvalidStateException, msg.str());
      } // decide on the subdomain this entity shall belong to
      debug << prefix << "    subdomain: " << subdomain << std::endl;
    } // walk the host grid
  } // void buildMsGrid(paramTree)

  Dune::shared_ptr< MsGridType > msGrid_;
}; // class Cube

template< class GridType >
const std::string Cube< GridType >::id = "grid.multiscale.provider.cube";


} // namespace Provider

} // namespace Multiscale

} // namespace Grid

} // namespace Dune

//#endif // HAVE_DUNE_STUFF

#endif // DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
