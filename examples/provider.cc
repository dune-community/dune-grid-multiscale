// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// system
#include <iostream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>

// dune-stuff
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/string.hh>

// dune-grid-multiscale
#include <dune/stuff/grid/provider.hh>
#include <dune/grid/multiscale/provider.hh>

const std::string id = "provider";

/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(std::string filename)
{
  Dune::Stuff::Common::testCreateDirectory(filename);
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[" << id << "]" << std::endl;
    file << "filename = " << id << ".grid" << std::endl;
    file << "grid = " << "grid.multiscale.provider.cube" << std::endl;
    file << "[grid.multiscale.provider.cube]" << std::endl;
    file << "lowerLeft = [0.0; 0.0; 0.0]" << std::endl;
    file << "upperRight = [1.0; 1.0; 1.0]" << std::endl;
    file << "numElements = [16; 16; 16]" << std::endl;
    file << "partitions = [2; 2; 2]" << std::endl;
    file << "oversamplingLayers = 1" << std::endl;
    file << "boundaryId = 7" << std::endl; // a cube from the factory gets the boundary ids 1 to 4 ind 2d and 1 to 6 in 3d?
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
} // void measureTiming()

template< class GlobalGridPartType, class CouplingGridPartType >
void testCoupling(const GlobalGridPartType& globalGridPart, const CouplingGridPartType& couplingGridPart, Dune::Stuff::Common::LogStream& out)
{
  // walk the grid part
  for (typename CouplingGridPartType::template Codim< 0 >::IteratorType entityIterator = couplingGridPart.template begin< 0 >();
       entityIterator != couplingGridPart.template end< 0 >();
       ++entityIterator) {
    const typename CouplingGridPartType::template Codim< 0 >::EntityType& entity = *entityIterator;
    const unsigned int entityIndex = globalGridPart.indexSet().index(entity);
    out << "entity " << entityIndex << ", neighbors " << std::flush;
    // walk the intersections
    for (typename CouplingGridPartType::IntersectionIteratorType intersectionIterator = couplingGridPart.ibegin(entity);
         intersectionIterator != couplingGridPart.iend(entity);
         ++intersectionIterator) {
      const typename CouplingGridPartType::IntersectionIteratorType::Intersection& intersection = *intersectionIterator;
      const typename CouplingGridPartType::IntersectionIteratorType::Intersection::EntityPointer neighborPtr = intersection.outside();
      const typename CouplingGridPartType::template Codim< 0 >::EntityType& neighbor = *neighborPtr;
      const unsigned int neighborIndex = globalGridPart.indexSet().index(neighbor);
      out << neighborIndex << " ";
    } // walk the intersections
    out << std::endl;
  } // walk the grid part
} // void testCoupling()

template< class GlobalGridPartType, class LocalGridPartType, class OutStreamType >
struct Inspect
{
  template< int dim, int codim >
  struct Codim
  {
    static void entities(const GlobalGridPartType& globalGridPart, const LocalGridPartType& localGridPart, const std::string prefix, OutStreamType& out)
    {
      // walk all codim entitites
      for (auto entityIt = localGridPart.template begin< codim >(); entityIt != localGridPart.template end< codim >(); ++entityIt) {
        const auto& entity = *entityIt;
        const auto& geometryType = entity.type();
        const auto globalIndex = globalGridPart.indexSet().index(entity);
        const auto localIndex = localGridPart.indexSet().index(entity);
        out << prefix << geometryType << ", global index " << globalIndex << ", local index " << localIndex << std::endl;
        out << prefix << "  corners: ";
        auto geometry = entity.geometry();
        for (int i = 0; i < (geometry.corners() - 1); ++i)
          out << "(" << geometry.corner(i) << "), ";
        out << "(" << geometry.corner(geometry.corners() - 1) << ")" << std::endl;
      } // walk all codim entitites
      // increase Codim
      Inspect< GlobalGridPartType, LocalGridPartType, OutStreamType >::Codim< dim, codim + 1 >::entities(globalGridPart, localGridPart, prefix, out);
    }
  }; // struct Codim
}; // struct Inspect

template< class GlobalGridPartType, class LocalGridPartType, class OutStreamType >
template< int codim >
struct Inspect< GlobalGridPartType, LocalGridPartType, OutStreamType >::Codim< codim, codim >
{
  static void entities(const GlobalGridPartType& globalGridPart, const LocalGridPartType& localGridPart, const std::string prefix, OutStreamType& out)
  {
    // walk all codim entitites
    for (auto entityIt = localGridPart.template begin< codim >(); entityIt != localGridPart.template end< codim >(); ++entityIt) {
      const auto& entity = *entityIt;
      const auto& geometryType = entity.type();
      const auto globalIndex = globalGridPart.indexSet().index(entity);
      const auto localIndex = localGridPart.indexSet().index(entity);
      out << prefix << geometryType << ", global index " << globalIndex << ", local index " << localIndex << std::endl;
      out << prefix << "  corners: (" << entity.geometry().center() << ")" << std::endl;
    } // walk all codim entitites
  }
};

int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);

    // parameter
    const std::string filename = id + ".param";
    ensureParamFile(filename);
    Dune::Stuff::Common::Configuration paramTree(argc, argv, filename);

    // logger
    Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO |
                                         Dune::Stuff::Common::LOG_CONSOLE |
                                         Dune::Stuff::Common::LOG_DEBUG);
    Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().debug();
    Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

    // timer
    Dune::Timer timer;

    // grid
    info << "setting up grid:" << std::endl;
    debug.suspend();

    typedef Dune::grid::Multiscale::ProviderFunctionbased<> MultiscaleGridProviderType;
    typedef typename MultiscaleGridProviderType::NonMultiscaleType NonMultiscaleType;
    typedef typename MultiscaleGridProviderType::GridType GridType;
    typedef double      DomainFieldType;
    const int           dimDomain = GridType::dimension;
    typedef int         RangeFieldType;
    typedef double      DiffusionRangeFieldType;
    typedef Dune::FieldMatrix< DiffusionRangeFieldType, 3, 3> RangeType;

    typedef typename Dune::Stuff::FunctionFromFile< DomainFieldType, dimDomain, RangeFieldType ,1 ,1 >
        FunctionFromFileScalarType;
    typedef typename Dune::Stuff::FunctionFromFile< DomainFieldType, dimDomain, DiffusionRangeFieldType ,dimDomain ,dimDomain >
        FunctionFromFileMatrixType;

    typedef Dune::Stuff::GridProviderCube< GridType > CubeProviderType;

    std::string fileForFunction = "brainmask.txt";
    std::string fileForDiffusiontensor = "tensors.txt";

    typedef typename Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
    DomainType lowerLeft(0.0);
    DomainType upperRight(117);
    upperRight[0]=117;
    upperRight[1]=142;
    upperRight[2]=124;
    std::vector< size_t > numElements = {117, 142, 124};

    // function representing the diffusion tensors
    std::shared_ptr< FunctionFromFileMatrixType > diffusion = std::make_shared< FunctionFromFileMatrixType >(fileForDiffusiontensor,
                                                                                                            lowerLeft,
                                                                                                            upperRight,
                                                                                                            numElements);

    // anlegen des Gitters, auf dem gerechnet werden soll mithilfe der Daten aus brain_mask
    std::shared_ptr< CubeProviderType > gridProvider = std::make_shared< CubeProviderType >(lowerLeft,
                                                                                            upperRight,
                                                                                            numElements);

    std::shared_ptr< GridType > grid(gridProvider->grid());
    gridProvider->visualize();

    // function representing the brainmask
    std::shared_ptr< FunctionFromFileScalarType > function = std::make_shared< FunctionFromFileScalarType >(fileForFunction,
                                                                                                            lowerLeft,
                                                                                                            upperRight,
                                                                                                            numElements);

    std::vector<double> partitions = {0.01};
    MultiscaleGridProviderType multiscaleProvider(grid, function, partitions);
    typedef MultiscaleGridProviderType::MsGridType MsGridType;
    const Dune::shared_ptr< const MsGridType > msGrid = multiscaleProvider.msGrid();
    //debug.resume();
    info << "  took " << timer.elapsed()
         << " sec (has " << msGrid->globalGridPart()->grid().size(0) << " elements, "
         << msGrid->size() << " subdomains)" << std::endl;

    info << "visualizing... " << std::flush;
    //debug.suspend();
    timer.reset();
    multiscaleProvider.visualize();
    info << " done (took " << timer.elapsed() << " sek)" << std::endl;
    //debug.resume();

    // visualization of the diffusion tensor
    // make a discrete function (make 9 functions) out of it and visualize it on the grid (which one?) with the vtkwriter
    // (auf dem grid vom gridProvider oder vom multiscaleGridProvider ? multi macht wohl mehr sinn)
    // evtl die diskreten Funktionen so anlegen, wie den multiscaleGridProvider?! d.h. über alle Entities laufen, im
    // center auswerten und die Werte der Matrix, die man erhält, in 9 Vektoren schreiben. d.h. als return value einen
    // Vektor von 9 Vektoren übergeben (leere Vektoren). (das ist aber wahrscheinlich sehr sehr langsam, wenn da mit so
    // vielen Schleifen und so vielen Zugriffen auf den Vektor gearbeitet wird.

    // write function diffusion into 9 vectors
    std::vector< double > diffusionVector1;
    std::vector< double > diffusionVector2;
    std::vector< double > diffusionVector3;
    std::vector< double > diffusionVector4;
    std::vector< double > diffusionVector5;
    std::vector< double > diffusionVector6;
    std::vector< double > diffusionVector7;
    std::vector< double > diffusionVector8;
    std::vector< double > diffusionVector9;

    typedef typename MsGridType::GlobalGridPartType GridPartType;
    const Dune::shared_ptr< const GridPartType> gridPart = msGrid->globalGridPart();
    RangeType ret;

    // walk the grid
    for (typename GridPartType::Codim< 0 >::IteratorType it = gridPart->begin< 0 >(); it != gridPart->end< 0 >(); ++it) {
      // get center of entity
      const auto& entity = *it;
      const DomainType center = entity.geometry().center();

      diffusion->evaluate(center, ret);
      // write the values of ret into the nine vectors --> this isn't very intelligent at the moment
      diffusionVector1.push_back(ret[0][0]);
      diffusionVector2.push_back(ret[0][1]);
      diffusionVector3.push_back(ret[0][2]);
      diffusionVector4.push_back(ret[1][0]);
      diffusionVector5.push_back(ret[1][1]);
      diffusionVector6.push_back(ret[1][2]);
      diffusionVector7.push_back(ret[2][0]);
      diffusionVector8.push_back(ret[2][1]);
      diffusionVector9.push_back(ret[2][2]);
    } // walk the grid

    typedef typename GridPartType::GridViewType GridViewType;
    std::cout << "visualizing diffusion tensors... " << std::flush;
    typedef Dune::VTKWriter< GridViewType > VTKWriterType;
    VTKWriterType vtkwriter((gridPart->gridView()));
    // write
    vtkwriter.addCellData(diffusionVector1, "diffusion_tensor[0][0]");
    vtkwriter.addCellData(diffusionVector2, "diffusion_tensor[0][1]");
    vtkwriter.addCellData(diffusionVector3, "diffusion_tensor[0][2]");
    vtkwriter.addCellData(diffusionVector4, "diffusion_tensor[1][0]");
    vtkwriter.addCellData(diffusionVector5, "diffusion_tensor[1][1]");
    vtkwriter.addCellData(diffusionVector6, "diffusion_tensor[1][2]");
    vtkwriter.addCellData(diffusionVector7, "diffusion_tensor[2][0]");
    vtkwriter.addCellData(diffusionVector8, "diffusion_tensor[2][1]");
    vtkwriter.addCellData(diffusionVector9, "diffusion_tensor[2][2]");
    vtkwriter.write("diffusiontensors", Dune::VTK::ascii);

//    info << "inspecting global grid part:" << std::endl;
//    const auto& globalGridPart = *(msGrid->globalGridPart());
//    Inspect< MsGridType::GlobalGridPartType, MsGridType::GlobalGridPartType, Dune::Stuff::Common::LogStream >
//        ::Codim< MultiscaleGridProviderType::dim, 0 >
//        ::entities(globalGridPart, globalGridPart, "  ", debug);

//    info << "inspecting local grid parts:" << std::endl;
//    for (unsigned int subdomain = 0; subdomain < msGrid->size(); ++subdomain) {
//      info << "subdomain " << subdomain << std::endl;
//      const auto& localGridPart = *(msGrid->localGridPart(subdomain));
//      Inspect< MsGridType::GlobalGridPartType, MsGridType::LocalGridPartType, Dune::Stuff::Common::LogStream >
//          ::Codim< MultiscaleGridProviderType::dim, 0 >
//          ::entities(globalGridPart, localGridPart, "  ", debug);
//    }

//    info << "inspecting boundary grid parts:" << std::endl;
//    for (unsigned int subdomain = 0; subdomain < msGrid->size(); ++subdomain) {
//      if (msGrid->boundary(subdomain)) {
//        info << "subdomain " << subdomain << std::endl;
//        const auto& boundaryGridPart = *(msGrid->boundaryGridPart(subdomain));
//        Inspect< MsGridType::GlobalGridPartType, MsGridType::BoundaryGridPartType, Dune::Stuff::Common::LogStream >
//            ::Codim< MultiscaleGridProviderType::dim, 0 >
//            ::entities(globalGridPart, boundaryGridPart, "  ", debug);
//      }
//    }

//    // time grid parts
//    info << "timing grid parts:" << std::endl;
//    typedef MsGridType::GlobalGridPartType GlobalGridPartType;
//    const Dune::shared_ptr< const GlobalGridPartType > globalGridPart = msGrid->globalGridPart();
//    measureTiming(*globalGridPart, info, "global");
//    typedef MsGridType::CouplingGridPartType CouplingGridPartType;
//    const unsigned int neighbor = *(msGrid->neighborsOf(0)->begin());
//    const Dune::shared_ptr< const CouplingGridPartType > couplingGridPart = msGrid->couplingGridPart(0, neighbor);
//    measureTiming(*couplingGridPart, info, "coupling (subdomain 0 with " + Dune::Stuff::Common::String::convertTo(neighbor) + ")");
//    typedef MsGridType::LocalGridPartType LocalGridPartType;
//    const Dune::shared_ptr< const LocalGridPartType > firstLocalGridPart = msGrid->localGridPart(0);
//    measureTiming(*firstLocalGridPart, info, "local (subdomain 0)");
//    for (unsigned int subdomain = 1; subdomain < msGrid->size(); ++subdomain) {
//      const Dune::shared_ptr< const LocalGridPartType > localGridPart = msGrid->localGridPart(subdomain);
//      measureTiming(*localGridPart, debug, "local (subdomain " + Dune::Stuff::Common::String::convertTo(subdomain) + ")");
//    }

//    // test coupling grid part
//    info << "testing coupling grid part:" << std::endl;
//    testCoupling(*globalGridPart, *couplingGridPart, info);

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
