
#include "config.h"

// system
#include <iostream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-stuff
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/filesystem.hh>

// dune-grid-multiscale
#include <dune/grid/multiscale/provider/cube.hh>

/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(std::string filename)
{
  Dune::Stuff::Common::Filesystem::testCreateDirectory(filename);
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[grid.multiscale.provider.cube]" << std::endl;
    file << "partitions.0 = 2" << std::endl;
    file << "partitions.1 = 2" << std::endl;
    file << "partitions.2 = 2" << std::endl;
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

    // logger
    Dune::Stuff::Common::Logger().Create(Dune::Stuff::Common::Logging::LOG_INFO |
                                         Dune::Stuff::Common::Logging::LOG_CONSOLE |
                                         Dune::Stuff::Common::Logging::LOG_DEBUG);
    Dune::Stuff::Common::Logging::LogStream& info = Dune::Stuff::Common::Logger().Info();
    Dune::Stuff::Common::Logging::LogStream& debug = Dune::Stuff::Common::Logger().Dbg();

    // timer
    Dune::Timer timer;

    // grid
    info << "setting up grid:" << std::endl;
    debug.Suspend();
    typedef Dune::grid::Multiscale::Provider::Cube< Dune::GridSelector::GridType > GridProviderType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, GridProviderType::id, filename);
    GridProviderType gridProvider(paramTree.sub(GridProviderType::id));
    typedef GridProviderType::GridType GridType;
    GridType& grid = gridProvider.grid();
    debug.Resume();
    info << "took " << timer.elapsed() << " sec, has " << grid.size(0) << " entities" << std::endl;

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
