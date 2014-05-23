// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define HAVE_DUNE_GRID_MULTISCALE 1

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#else
# error This test requires ALUGrid!
#endif

#include <dune/grid/multiscale/provider/cube.hh>

using namespace Dune;
using namespace Dune::grid::Multiscale;

int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);

    typedef ALUConformGrid< 2, 2 > GridType;
    typedef Providers::Cube< GridType > GridProviderType;
    auto config = GridProviderType::default_config();
    config["lower_left"]     = "-1";
    config["upper_right"]    = "1";
    config["num_elements"]   = "4";
    config["num_partitions"] = "3";
    auto grid_provider = GridProviderType::create(config);
    grid_provider->visualize("grid");

    const auto ms_grid = grid_provider->ms_grid();
    const auto global_grid_part = ms_grid->globalGridPart();
    const auto coupling_grid_part = ms_grid->couplingGridPart(0, 1);
    const auto entity_it_end = coupling_grid_part->template end< 0 >();
    for (auto entity_it = coupling_grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const size_t entity_id = global_grid_part->indexSet().index(entity);
      if (entity_id == 6) {
        const auto intersection_it_end = coupling_grid_part->iend(entity);
        for (auto intersection_it = coupling_grid_part->ibegin(entity);
             intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          const size_t local_intersection_id = intersection.indexInInside();
          std::cout << "entity " << entity_id << ", intersection " << local_intersection_id << std::endl;
        }
      }
    }

  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error:\n" << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
