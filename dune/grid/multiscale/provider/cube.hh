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
//#include <dune/stuff/common/logging.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/default.hh>

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

  typedef Dune::grid::Multiscale::Default< GridType > MsGridType;

  typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;

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

  Cube(Dune::ParameterTree paramTree = Dune::ParameterTree())
    : BaseType(paramTree)
  {
    buildMsGrid(paramTree);
  }

  const MsGridType& msGrid() const
  {
    return *msGrid_;
  }

private:
  void buildMsGrid(const Dune::ParameterTree& paramTree)
  {
    // preparations
    const CoordinateType& lowerLeft = BaseType::lowerLeft();
    const CoordinateType& upperRight = BaseType::upperRight();
    msGrid_ = Dune::shared_ptr< MsGridType >(new MsGridType(BaseType::grid()));
    msGrid_->prepare();

    //// debug output
    //const std::string prefix = paramTree.get("prefix", "");
    //Dune::ParameterTree indentTree;
    //indentTree["prefix"] = prefix + "  ";
    //Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();
    //debug << prefix << id << ":" << std::endl;
    //debug << prefix << "  lowerLeft: " << lowerLeft << std::endl;
    //debug << prefix << "  upperRight: " << upperRight << std::endl;
    //debug << prefix << "  partitions per dim: ";

    // get number of subdomains per direction
    std::vector< unsigned int > partitions(dim, 0);
    unsigned int totalSubdomains = 1;
    for (unsigned int d = 0; d < dim; ++d) {
      std::stringstream key;
      key << "partitions." << d;
      partitions[d] = paramTree.get(key.str(), 1);
      totalSubdomains *= partitions[d];
      //debug << partitions[d] << " ";
    }
    //debug << std::endl;

    // grid part
    typedef typename MsGridType::GlobalGridPartType GridPartType;
    const GridPartType& gridPart = msGrid_->globalGridPart();

    // walk the grid
    for (typename GridPartType::template Codim< 0 >::IteratorType it = gridPart.template begin< 0 >();
        it != gridPart.template end< 0 >();
        ++it) {
      // get center of entity
      typename GridType::LeafGridView::template Codim< 0 >::Iterator::Entity& entity = *it;
      const CoordinateType center = entity.geometry().global(entity.geometry().center());
      //debug << prefix << "  entity " << gridPart.indexSet().index(entity) << " (" << center << ")";

      // decide on the subdomain this entity shall belong to
      std::vector< unsigned int > whichPartition;
      for (unsigned int d = 0; d < dim; ++d) {
        whichPartition.push_back(std::floor(partitions[d]*((center[d] - lowerLeft[d])/(upperRight[d] - lowerLeft[d]))));
        //debug << prefix << "    partition[" << d  << "]: " << whichPartition[d] << std::endl;
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
      //debug << ", subdomain " << subdomain << std::endl;

      // add entity to subdomain
      msGrid_->add(entity, subdomain);

    } // walk the grid

    msGrid_->finalize(/*indentTree*/);
    return;
  } // void buildMsGrid(const Dune::ParameterTree& paramTree)

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
