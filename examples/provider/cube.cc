
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
#include <dune/stuff/common/string.hh>

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
    file << "numElements.0 = 6" << std::endl;
    file << "numElements.1 = 6" << std::endl;
    file << "numElements.2 = 6" << std::endl;
    file << "boundaryId = 7" << std::endl;                  // a cube from the factory gets the boundary ids 1 to 4 ind 2d and 1 to 6 ? in 3d
    file << "filename = msGrid_visualization" << std::endl;
    file << "partitions.0 = 2" << std::endl;
    file << "partitions.1 = 2" << std::endl;
    file << "partitions.2 = 2" << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

template< class GridPartType >
void measureTiming(const GridPartType& gridPart, Dune::Stuff::Common::LogStream& out, const std::string name = "")
{
  out << "  walking " << name << " grid part... " << std::flush;
  Dune::Timer timer;
  unsigned int elements = 0;
  for (typename GridPartType::template Codim< 0 >::IteratorType it = gridPart.template begin< 0 >();
       it != gridPart.template end< 0 >();
       ++it)
    ++elements;
  out << " done (has " << elements << " elements, took " << timer.elapsed() << " sek)" << std::endl;
  return;
}

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
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO |
                                         Dune::Stuff::Common::LOG_CONSOLE |
                                         Dune::Stuff::Common::LOG_DEBUG);
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().info();
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    // timer
    Dune::Timer timer;

    // grid
    info << "setting up grid: " << std::endl;
    debug.suspend();
    typedef Dune::grid::Multiscale::Provider::Cube< Dune::GridSelector::GridType > GridProviderType;
    Dune::Stuff::Common::Parameter::Tree::assertSub(paramTree, GridProviderType::id, id);
    GridProviderType gridProvider(paramTree.sub(GridProviderType::id));
    typedef GridProviderType::MsGridType MsGridType;
    const MsGridType& msGrid = gridProvider.msGrid();
    debug.resume();
    info << "  took " << timer.elapsed()
         << " sec (has " << msGrid.globalGridPart()->grid().size(0) << " elements, "
         << msGrid.size() << " subdomains)" << std::endl;

    info << "visualizing... " << std::flush;
    debug.suspend();
    timer.reset();
    msGrid.visualize();
    debug.resume();
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;

    // time grid parts
    info << "timing grid parts:" << std::endl;
    typedef MsGridType::GlobalGridPartType GlobalGridPartType;
    const Dune::shared_ptr< const GlobalGridPartType > globalGridPart = msGrid.globalGridPart();
    measureTiming(*globalGridPart, info, "global");
    typedef MsGridType::LocalGridPartType LocalGridPartType;
    const Dune::shared_ptr< const LocalGridPartType > firstLocalGridPart = msGrid.localGridPart(0);
    measureTiming(*firstLocalGridPart, info, "local (subdomain 0)");
    for (unsigned int subdomain = 1; subdomain < msGrid.size(); ++subdomain) {
      const Dune::shared_ptr< const LocalGridPartType > localGridPart = msGrid.localGridPart(subdomain);
      measureTiming(*localGridPart, debug, "local (subdomain " + Dune::Stuff::Common::String::convertTo(subdomain) + ")");
    }

    // if we came that far we can as well be happy about it
    return 0;
  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch(...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
