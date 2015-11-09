// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#elif HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#include <dune/grid/multiscale/glued.hh>

using namespace Dune;

template <class LocalGridType>
struct GluedMultiscaleGrid : public ::testing::Test
{
  void works()
  {
    typedef YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> MacroGridType;
    Stuff::Grid::Providers::Cube<MacroGridType> macro_grid({0., 0.}, {1., 1.}, {3, 3});

    grid::Multiscale::Glued<MacroGridType, LocalGridType> multiscale_grid(macro_grid, 2);

    size_t num_couplings        = 0;
    const auto& macro_grid_view = multiscale_grid.macro_grid_view();
    for (auto&& macro_entity : elements(macro_grid_view)) {
      const auto entity_index = macro_grid_view.indexSet().index(macro_entity);
      for (auto&& macro_intersection : intersections(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor = macro_intersection.outside();
          const auto neighbor_index = macro_grid_view.indexSet().index(macro_neighbor);
          const auto& coupling = multiscale_grid.coupling(macro_entity, 1, macro_neighbor, 2);
          EXPECT_EQ(4, coupling.size()) << "entity " << entity_index << ", neighbor " << neighbor_index;
          ++num_couplings;
        }
      }
    }
    // because we have 3x3 macro entities -> 12 coupling intersections (twice, since we are looking from each side)
    EXPECT_EQ(num_couplings, 24);
  } // ... works(...)
};  // struct GluedMultiscaleGrid

typedef ::testing::Types<YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> // clang-format off
#if HAVE_DUNE_ALUGRID
                        , ALUGrid<2, 2, simplex, nonconforming>
#elif HAVE_ALUGRID
                        , ALUGrid<2, 2, simplex, conforming>
#endif
                        > LocalGridTypes; // clang-format on

TYPED_TEST_CASE(GluedMultiscaleGrid, LocalGridTypes);
TYPED_TEST(GluedMultiscaleGrid, works) { this->works(); }
