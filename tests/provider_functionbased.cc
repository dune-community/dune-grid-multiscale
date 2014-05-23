#include <dune/stuff/test/test_common.hh>

#include <iostream>
#include <fstream>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/fromfile.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>

// dune-grid-multiscale
#include <dune/stuff/grid/provider.hh>
#include <dune/grid/multiscale/provider.hh>

using namespace Dune::Stuff;


TEST(PROVIDER_FUNCTIONBASED, All) {
  // logger
  Dune::Stuff::Common::Logger().create(Dune::Stuff::Common::LOG_INFO |
                                       Dune::Stuff::Common::LOG_CONSOLE |
                                       Dune::Stuff::Common::LOG_DEBUG);
  Dune::Stuff::Common::LogStream& info = Dune::Stuff::Common::Logger().debug();
  Dune::Stuff::Common::LogStream& debug = Dune::Stuff::Common::Logger().debug();

  // timer
  Dune::Timer timer;

  // create the textfiles from which the functions are gonna be created
  std::string fileForFunction = "brainmask.txt";
  std::string fileForDiffusiontensor = "tensors.txt";

  std::ofstream myfile;
  myfile.open (fileForFunction);
  myfile << "1 1 1 1\n";
  myfile << "1 1 2 0\n";
  myfile << "1 1 3 0\n";
  myfile << "1 2 1 0\n";
  myfile << "1 2 2 1\n";
  myfile << "1 2 3 0\n";
  myfile << "1 3 1 1\n";
  myfile << "1 3 2 0\n";
  myfile << "1 3 3 0\n";
  myfile << "2 1 1 1\n";
  myfile << "2 1 2 1\n";
  myfile << "2 1 3 1\n";
  myfile << "2 2 1 1\n";
  myfile << "2 2 2 1\n";
  myfile << "2 2 3 0\n";
  myfile << "2 3 1 1\n";
  myfile << "2 3 2 1\n";
  myfile << "2 3 3 0\n";
  myfile << "3 1 1 1\n";
  myfile << "3 1 2 1\n";
  myfile << "3 1 3 0\n";
  myfile << "3 2 1 1\n";
  myfile << "3 2 2 1\n";
  myfile << "3 2 3 0\n";
  myfile << "3 3 1 1\n";
  myfile << "3 3 2 0\n";
  myfile << "3 3 3 0\n";
  myfile.close();

  myfile.open (fileForDiffusiontensor);
  myfile << "1 1 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 1 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 1 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 2 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 2 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 2 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 3 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 3 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "1 3 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 1 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 1 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 1 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 2 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 2 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 2 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 3 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 3 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "2 3 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 1 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 1 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 1 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 2 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 2 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 2 3 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 3 1 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 3 2 0 0 0 0 0 0 0 0 0\n";
  myfile << "3 3 3 0 0 0 0 0 0 0 0 0\n";
  myfile.close();

  typedef Dune::grid::Multiscale::ProviderFunctionbased<> MultiscaleGridProviderType;
  typedef typename MultiscaleGridProviderType::NonMultiscaleType NonMultiscaleType;
  typedef typename MultiscaleGridProviderType::GridType GridType;
  typedef double      DomainFieldType;
  const int           dimDomain = GridType::dimension;
  typedef double      RangeFieldType;
  typedef double      DiffusionRangeFieldType;
  typedef Dune::FieldMatrix< DiffusionRangeFieldType, 3, 3> RangeType;

  typedef typename Dune::Stuff::FunctionFromFile< DomainFieldType, dimDomain, RangeFieldType ,1 ,1 >
      FunctionFromFileScalarType;
  typedef typename Dune::Stuff::FunctionFromFile< DomainFieldType, dimDomain, DiffusionRangeFieldType ,dimDomain ,dimDomain >
      FunctionFromFileMatrixType;

  typedef Dune::Stuff::GridProviderCube< GridType > CubeProviderType;

  typedef typename Dune::FieldVector<DomainFieldType, dimDomain> DomainType;
  DomainType lowerLeft(0.0);
  DomainType upperRight(3.0);
  std::vector< size_t > numElements = {3, 3, 3};

  // function representing the diffusion tensors
  std::shared_ptr< FunctionFromFileMatrixType > diffusion = std::make_shared< FunctionFromFileMatrixType >(fileForDiffusiontensor,
                                                                                                          lowerLeft,
                                                                                                          upperRight,
                                                                                                          numElements);

  // creating the grid with the help of the data from brainmask.txt
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
  info << "  took " << timer.elapsed()
       << " sec (has " << msGrid->globalGridPart()->grid().size(0) << " elements, "
       << msGrid->size() << " subdomains)" << std::endl;

  info << "visualizing... " << std::flush;
  timer.reset();
  multiscaleProvider.visualize();
  info << " done (took " << timer.elapsed() << " sek)" << std::endl;

  // visualization of the diffusion tensor
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
    // write the values of ret into the nine vectors
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
  info << "visualizing diffusion tensors... " << std::flush;
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
}

int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch(Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch(...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
}
