// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH

#include <vector>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/dgfparser.hh>

//#include <dune/stuff/grid/provider/cube.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/grid/multiscale/factory/default.hh>
#include <dune/grid/multiscale/default.hh>

#include "interface.hh"

namespace Dune {
namespace grid {
namespace Multiscale {
namespace Providers {


static inline std::string cube_gridprovider_id()
{
  return "grid.multiscale.provider.cube";
}


static inline XT::Common::Configuration cube_gridprovider_default_config()
{
  XT::Common::Configuration config = XT::Grid::cube_gridprovider_default_config();
  config["num_partitions"] = "[2 2 2]";
  config["oversampling_layers"] = "0";
  config["inner_boundary_segment_index"] = XT::Common::to_string(std::numeric_limits<size_t>::max() - 42);
  return config;
}


#if HAVE_DUNE_FEM

template <class GridType>
class CubeGridProviderFactory
{
  static_assert(XT::Grid::is_grid<GridType>::value, "");

public:
  typedef Default<GridType> MsGridType;

  static std::string static_id()
  {
    return cube_gridprovider_id();
  }

  static XT::Common::Configuration default_config()
  {
    return cube_gridprovider_default_config();
  }

  static DefaultProvider<GridType> create(const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
                                          const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
                                          const std::array<unsigned int, GridType::dimension>& num_elements,
                                          const unsigned int num_refinements,
                                          const std::array<unsigned int, GridType::dimension>& overlap_size,
                                          const std::array<unsigned int, GridType::dimension>& num_partitions,
                                          const size_t num_oversampling_layers,
                                          const size_t inner_boundary_segment_index)
  {
    auto grid = XT::Grid::make_cube_grid<GridType>(lower_left, upper_right, num_elements, num_refinements, overlap_size)
                    .grid_ptr();

    typedef Factory::Default<GridType> MsGridFactoryType;
    const size_t neighbor_recursion_level = Factory::NeighborRecursionLevel<GridType>::compute();
    // prepare
    MsGridFactoryType factory(*grid, inner_boundary_segment_index);
    factory.prepare();
    // global grid part
    const auto global_grid_part = factory.globalGridPart();
    // walk the grid
    const auto entity_it_end = global_grid_part->template end<0>();
    for (auto entity_it = global_grid_part->template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      // get center of entity
      const auto& entity = *entity_it;
      const auto center = entity.geometry().center();
      // decide on the subdomain this entity shall belong to
      std::vector<size_t> whichPartition(GridType::dimension, 0);
      for (size_t dd = 0; dd < GridType::dimension; ++dd)
        whichPartition[dd] =
            (std::min((unsigned int)(std::floor(
                          num_partitions[dd] * ((center[dd] - lower_left[dd]) / (upper_right[dd] - lower_left[dd])))),
                      num_partitions[dd] - 1));
      size_t subdomain = 0;
      if (GridType::dimension == 1)
        subdomain = whichPartition[0];
      else if (GridType::dimension == 2)
        subdomain = whichPartition[0] + whichPartition[1] * num_partitions[0];
      else if (GridType::dimension == 3)
        subdomain = whichPartition[0] + whichPartition[1] * num_partitions[0]
                    + whichPartition[2] * num_partitions[1] * num_partitions[0];
      else
        DUNE_THROW(Dune::NotImplemented,
                   "ERROR in " << static_id() << ": not implemented for grid dimensions other than 1, 2 or 3!");
      // add entity to subdomain
      factory.add(entity, subdomain /*, prefix + "  ", out*/);
    } // walk the grid
    // finalize
    factory.finalize(num_oversampling_layers, neighbor_recursion_level /*, prefix + "  ", out*/);
    // be done with it
    return DefaultProvider<GridType>(grid, factory.createMsGrid());
  } // ... create(...)
}; // class CubeGridProviderFactory


template <class GridType>
typename std::enable_if<XT::Grid::is_grid<GridType>::value, DefaultProvider<GridType>>::type make_cube_grid(
    const FieldVector<typename GridType::ctype, GridType::dimension>& lower_left,
    const FieldVector<typename GridType::ctype, GridType::dimension>& upper_right,
    const std::array<unsigned int, GridType::dimension> num_elements =
        cube_gridprovider_default_config().template get<std::vector<unsigned int>>("num_elements")[0],
    const unsigned int num_refinements =
        cube_gridprovider_default_config().template get<unsigned int>("num_refinements"),
    const std::array<unsigned int, GridType::dimension> overlap_size =
        XT::Common::make_array<unsigned int, GridType::dimension>(
            cube_gridprovider_default_config().template get<std::vector<unsigned int>>("overlap_size")),
    const std::array<unsigned int, GridType::dimension> num_partitions =
        XT::Common::make_array<unsigned int, GridType::dimension>(
            cube_gridprovider_default_config().template get<std::vector<unsigned int>>("num_partitions")),
    const size_t num_oversampling_layers = cube_gridprovider_default_config().template get<size_t>("num_refinements"),
    const size_t inner_boundary_segment_index = cube_gridprovider_default_config().template get<int>("inner_boundary_segment_index"))
{
  return CubeGridProviderFactory<GridType>::create(lower_left,
                                                   upper_right,
                                                   num_elements,
                                                   num_refinements,
                                                   overlap_size,
                                                   num_partitions,
                                                   num_oversampling_layers,
                                                   inner_boundary_segment_index);
}


#else // HAVE_DUNE_FEM

// template <class GridImp>
// class Cube
//{
//  static_assert(AlwaysFalse<GridImp>::value, "Your are missing dune-fem!");
//};

#endif // HAVE_DUNE_FEM

} // namespace Providers
} // namespace Multiscale
} // namespace Grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
