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
#include <dune/grid/multiscale/filtered/cube.hh>
#include <dune/grid/multiscale/filtered.hh>

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
  typedef Dune::grid::Multiscale::Filter::Cube< GridType > FilterType;

//  template< int dim >
//  struct P0Layout
//  {
//    template< class GeometryType >
//    bool contains(GeometryType& geometry)
//    {
//      if (geometry.dim() == dim)
//        return true;
//      return false;
//    }
//  }; // layout class for codim 0 mapper

public:
  typedef typename Dune::grid::Multiscale::Filtered< GridType, FilterType > MsGridType;

  Cube(Dune::ParameterTree paramTree = Dune::ParameterTree())
    : BaseType(paramTree)
  {
    buildMsGridFromFineGrid(paramTree);
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
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().dbg();

//    const std::string filename = paramTree.get("visualize", id + ".msgrid");

//    // mapper
//    Dune::LeafMultipleCodimMultipleGeomTypeMapper< GridType, P0Layout > mapper(BaseType::grid());
//    std::vector< double > data(mapper.size());

    // walk the subdomains
    typedef typename MsGridType::LocalGridPartType LocalGridPartType;
    for (unsigned int subdomain = 0; subdomain < msGrid_->size(); ++subdomain) {
        debug << prefix << "subdomain " << subdomain << std::endl;
      // walk the local grid part
      typename LocalGridPartType::template Codim< 0 >::IteratorType itEnd = msGrid_->localGridPart(subdomain)->template end< 0 >();
      for (typename LocalGridPartType::template Codim< 0 >::IteratorType it = msGrid_->localGridPart(subdomain)->template begin< 0 >();
          it != itEnd;
          ++it) {
        // local entity
        typename LocalGridPartType::template Codim< 0 >::IteratorType::Entity& entity = *it;
        Dune::Stuff::Grid::Entity::print(entity, debug, "  ");
//        data[mapper.map(hostEntity)] = double(subdomain);
      } // walk the local grid
    } // walk the subdomains
//    // write to vtk
//    Dune::VTKWriter< typename GridType::LeafGridView > vtkwriter(BaseType::grid().leafView());
//    vtkwriter.addCellData(data, "subdomain");
//    vtkwriter.write(filename, Dune::VTK::ascii);
    return;
  } // void visualize(const Dune::ParameterTree& paramTree) const

private:
//  void buildMsGridFromCoarseGrid(const Dune::ParameterTree& paramTree)
//  {
//    // preparations
//    const std::string prefix = paramTree.get("prefix", "");
//    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().dbg();
//    debug << prefix << id << ":" << std::endl;
//    debug << prefix << "  creating " << BaseType::grid().levelView(0).size(0) << " local grids " << std::flush;
//    const unsigned int refineLevel = paramTree.get("refineLevel", 1);
//    BaseType::grid().globalRefine(refineLevel);
//    debug << "with approx. " << BaseType::grid().leafView().size(0)/BaseType::grid().levelView(0).size(0)
//          << " elements each... " << std::flush;

//    // create factory
//    MsGridFactoryType msGridFactory(BaseType::grid());
//    msGridFactory.prepare();

//    // walk the coarse grid
//    unsigned int subdomain = 0;
//    typename GridType::LevelGridView::template Codim< 0 >::Iterator itEnd = BaseType::grid().levelView(0).template end< 0 >();
//    for (typename GridType::LevelGridView::template Codim< 0 >::Iterator it = BaseType::grid().levelView(0).template begin< 0 >();
//         it != itEnd;
//         ++it) {
//      // create subgrid from coarse entity
//      typename GridType::LevelGridView::template Codim< 0 >::Iterator::Entity& entity = *it;
//      msGridFactory.add(entity, subdomain);
//      ++subdomain;
//    } // walk the coarse grid

//    // finalize factory
//    msGrid_ = msGridFactory.finalize();

//    debug << "done" << std::flush;
//    return;
//  } // void buildMsGridFromCoarseGrid(const Dune::ParameterTree& paramTree)

  void buildMsGridFromFineGrid(const Dune::ParameterTree& paramTree)
  {
    // preparations
    const CoordinateType& lowerLeft = BaseType::lowerLeft();
    const CoordinateType& upperRight = BaseType::upperRight();
    msGrid_ = Dune::shared_ptr< MsGridType >(new MsGridType(BaseType::grid()));

    // debug output
    const std::string prefix = paramTree.get("prefix", "");
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    debug << prefix << id << ":" << std::endl;
    debug << prefix << "  lowerLeft: " << lowerLeft << std::endl;
    debug << prefix << "  upperRight: " << upperRight << std::endl;
    debug << prefix << "  partitions per dim: ";

    // get number of subdomains per direction
    std::vector< unsigned int > partitions(dim, 0);
    unsigned int totalSubdomains = 1;
    for (unsigned int d = 0; d < dim; ++d) {
      std::stringstream key;
      key << "partitions." << d;
      partitions[d] = paramTree.get(key.str(), 1);
      totalSubdomains *= partitions[d];
      debug << partitions[d] << " ";
    }
    debug << std::endl;

    // create subdomain cubes
    std::vector< CoordinateType > lowerLefts(totalSubdomains, CoordinateType(0.0));
    std::vector< CoordinateType > upperRights(totalSubdomains, CoordinateType(0.0));
    const CoordinateType length = upperRight - lowerLeft;
    unsigned int subdomain = 0;
    if (dim == 1) {
      for (unsigned int i = 0; i < partitions[0]; ++i) {
        lowerLefts[subdomain][0] = i*(length[0]/partitions[0]);
        upperRights[subdomain][0] = (i + 1)*(length[0]/partitions[0]);
        ++ subdomain;
      }
    } else if (dim == 2) {
      for (unsigned int i = 0; i < partitions[0]; ++i) {
        for (unsigned int j = 0; j < partitions[1]; ++j) {
          lowerLefts[subdomain][0] = i*(length[0]/partitions[0]);
          lowerLefts[subdomain][1] = j*(length[1]/partitions[1]);
          upperRights[subdomain][0] = (i + 1)*(length[0]/partitions[0]);
          upperRights[subdomain][1] = (j + 1)*(length[1]/partitions[1]);
          ++ subdomain;
        }
      }
    } else if (dim == 3) {
      for (unsigned int i = 0; i < partitions[0]; ++i) {
        for (unsigned int j = 0; j < partitions[1]; ++j) {
          for (unsigned int k = 0; k < partitions[2]; ++ k) {
            lowerLefts[subdomain][0] = i*(length[0]/partitions[0]);
            lowerLefts[subdomain][1] = j*(length[1]/partitions[1]);
            lowerLefts[subdomain][2] = k*(length[2]/partitions[2]);
            upperRights[subdomain][0] = (i + 1)*(length[0]/partitions[0]);
            upperRights[subdomain][1] = (j + 1)*(length[1]/partitions[1]);
            upperRights[subdomain][2] = (k + 1)*(length[2]/partitions[2]);
            ++ subdomain;
          }
        }
      }
    } else {
      std::stringstream msg;
      msg << "Error in " << id << ": only implemented for dimensions 1, 2 and 3!";
      DUNE_THROW(Dune::InvalidStateException, msg.str());
    }

    for (unsigned int s = 0; s < totalSubdomains; ++s) {
      debug << prefix << "  subdomain:  " << s << std::endl;
      debug << prefix << "  lowerLeft:  " << lowerLefts[s] << std::endl;
      debug << prefix << "  upperRight: " << upperRights[s] << std::endl;
    }

    // create filter for each subdomain and register at the factory
    for (unsigned int subdomain = 0; subdomain < totalSubdomains; ++subdomain) {
      Dune::shared_ptr< FilterType > filter = Dune::shared_ptr< FilterType >(new FilterType(lowerLefts[subdomain], upperRights[subdomain]));
      msGrid_->create(filter, subdomain);
    }
    msGrid_->finalize();

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
