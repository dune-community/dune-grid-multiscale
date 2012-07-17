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

// dune-grid
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// dune-stuff
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/entity.hh>

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
  typedef Dune::grid::Multiscale::Factory::Subgrid::FromCoarseGrid< GridType > MsGridFactoryType;

  template< int dim >
  struct P0Layout
  {
    template< class GeometryType >
    bool contains(GeometryType& geometry)
    {
      if (geometry.dim() == dim)
        return true;
      return false;
    }
  }; // layout class for codim 0 mapper

public:
  typedef typename MsGridFactoryType::MsGridType MsGridType;

  Cube(Dune::ParameterTree paramTree = Dune::ParameterTree())
    : BaseType(paramTree)
  {
    buildMsGridFromCoarseGrid(paramTree);
  }

  MsGridType& msGrid()
  {
    return *msGrid_;
  }

  const MsGridType& msGrid() const
  {
    return *msGrid_;
  }

  void visualize(const Dune::ParameterTree& paramTree) const
  {
    // debug output
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::Logging::LogStream& debug = Dune::Stuff::Common::Logger().Dbg();

    const std::string filename = paramTree.get("visualize", id + ".msgrid");

    // mapper
    Dune::LeafMultipleCodimMultipleGeomTypeMapper< GridType, P0Layout > mapper(BaseType::grid());
    std::vector< double > data(mapper.size());

    // walk the subdomains
    for (unsigned int subdomain = 0; subdomain < msGrid_->numSubdomains(); ++subdomain) {
      // walk the local grid part
      typename MsGridType::LocalGridPartType::template Codim< 0 >::IteratorType itEnd = msGrid_->localGridPart(subdomain).template end< 0 >();
      for (typename MsGridType::LocalGridPartType::template Codim< 0 >::IteratorType it = msGrid_->localGridPart(subdomain).template begin< 0 >();
          it != itEnd;
          ++it) {
        // local entity
        typename MsGridType::LocalEntityType& localEntity = *it;
        Dune::Stuff::Grid::Entity::print(localEntity, debug, "  ");


//        data[mapper.map(hostEntity)] = double(subdomain);
      } // walk the local grid
    } // walk the subdomains
    // write to vtk
    Dune::VTKWriter< typename GridType::LeafGridView > vtkwriter(BaseType::grid().leafView());
    vtkwriter.addCellData(data, "subdomain");
    vtkwriter.write(filename, Dune::VTK::ascii);
    return;
  } // void visualize(const Dune::ParameterTree& paramTree) const

private:
  void buildMsGridFromCoarseGrid(const Dune::ParameterTree& paramTree)
  {
    // preparations
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::Logging::LogStream& debug = Dune::Stuff::Common::Logger().Dbg();
    debug << prefix << id << ":" << std::endl;
    debug << prefix << "  creating " << BaseType::grid().levelView(0).size(0) << " local grids " << std::flush;
    const unsigned int refineLevel = paramTree.get("refineLevel", 1);
    BaseType::grid().globalRefine(refineLevel);
    debug << "with approx. " << BaseType::grid().leafView().size(0)/BaseType::grid().levelView(0).size(0)
          << " elements each... " << std::flush;

    // create factory
    MsGridFactoryType msGridFactory(BaseType::grid());
    msGridFactory.prepare();

    // walk the coarse grid
    unsigned int subdomain = 0;
    typename GridType::LevelGridView::template Codim< 0 >::Iterator itEnd = BaseType::grid().levelView(0).template end< 0 >();
    for (typename GridType::LevelGridView::template Codim< 0 >::Iterator it = BaseType::grid().levelView(0).template begin< 0 >();
         it != itEnd;
         ++it) {
      // create subgrid from coarse entity
      typename GridType::LevelGridView::template Codim< 0 >::Iterator::Entity& entity = *it;
      msGridFactory.add(entity, subdomain);
      ++subdomain;
    } // walk the coarse grid

    // finalize factory
    msGrid_ = msGridFactory.finalize();

    debug << "done" << std::flush;
    return;
  } // void buildMsGridFromCoarseGrid(const Dune::ParameterTree& paramTree)

  void buildMsGridFromFineGrid(const Dune::ParameterTree& paramTree)
  {
    // preparations
    const CoordinateType& lowerLeft = BaseType::lowerLeft();
    const CoordinateType& upperRight = BaseType::upperRight();

    // debug output
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::Logging::LogStream& debug = Dune::Stuff::Common::Logger().Dbg();
    debug << prefix << id << ".buildMsGridFromFineGrid:" << std::endl;
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

      // add entity to corresponding subdomain
      msGridFactory.add(entity, subdomain);
    } // walk the host grid

    // finalize
    msGridFactory.finalize();
    return;
  } // void buildMsGridFromFineGrid(const Dune::ParameterTree& paramTree)

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
