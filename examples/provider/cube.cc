
#include "config.h"

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-stuff
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/parameter/tree.hh>

/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(std::string filename)
{
  // only write param file if there is none
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[helper-tools.grid.provider.cube]" << std::endl;
    file << "level = 4" << std::endl;
    file << "[detailed-solvers.stationary.linear.elliptic.model]" << std::endl;
    file << "diffusion.variable = x" << std::endl;
    file << "diffusion.expression.0 = 1.0"  << std::endl;
    file << "diffusion.expression.1 = 1.0"  << std::endl;
    file << "diffusion.expression.2 = 1.0"  << std::endl;
    file << "diffusion.order = 0"  << std::endl;
    file << "force.variable = x" << std::endl;
    file << "force.expression.0 = 1.0"  << std::endl;
    file << "force.expression.1 = 1.0"  << std::endl;
    file << "force.expression.2 = 1.0"  << std::endl;
    file << "force.order = 0"  << std::endl;
    file << "[detailed-solvers.stationary.linear.elliptic.continuousgalerkin]" << std::endl;
    file << "init.verbose = true" << std::endl;
    file << "solve.verbose = true" << std::endl;
    file << "solve.type = eigen.bicgstab.incompletelut" << std::endl;
    file << "solve.maxIter = 5000"  << std::endl;
    file << "solve.precision = 1e-12"  << std::endl;
    file << "visualize.verbose = true"  << std::endl;
    file << "visualize.name = solution"  << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);

    // parameter
    const std::string id = "cube";
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::ParameterTree paramTree = Dune::Stuff::Common::Parameter::Tree::init(argc, argv, filename);

    // timer
    Dune::Timer timer;

    // grid
    std::cout << "setting up grid:" << std::endl;
    typedef Dune::Stuff::Grid::Provider::UnitCube< Dune::GridSelector::GridType > GridProviderType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, GridProviderType::id, filename);
    GridProviderType gridProvider(paramTree.sub(GridProviderType::id));
    typedef GridProviderType::GridType GridType;
    GridType& grid = gridProvider.grid();
    typedef Dune::LeafGridPart< GridType > GridPartType;
    GridPartType gridPart(grid);
    typedef Dune::GridPartView< GridPartType > GridViewType;
    GridViewType gridView(gridPart);
    std::cout << "took " << timer.elapsed() << " sec, has " << gridView.size(0) << " entities" << std::endl;

    // if we came that far we can as well be happy about it
    return 0;
  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch( ... ) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
