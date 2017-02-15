// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/common/version.hh>

// alugrid
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#elif HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
// yaspgrid or sgrid
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
#include <dune/grid/yaspgrid.hh>
#else
#include <dune/grid/sgrid.hh>
#endif
// elements or DSC::entityRange
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
#include <dune/grid/common/rangegenerators.hh>
#else
#include <dune/stuff/common/ranges.hh>
#endif

#include <dune/grid/multiscale/glued.hh>

using namespace Dune;


template<class G, bool anything = true>
struct num_refinements
{
  static_assert(AlwaysFalse<G>::value, "");
};

template<int dim, int dimworld, typename _ctype, bool anything>
struct num_refinements<SGrid<dim, dimworld, _ctype>, anything>
{
  int operator()()
  {
    return 2;
  }
};

template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm, bool anything>
struct num_refinements<ALUGrid<dim, dimworld, elType, refineType, Comm>, anything>
{
  int operator()()
  {
    return 4;
  }
};

template <class LocalGridType>
struct GluedMultiscaleGrid : public ::testing::Test
{
  void works()
  {
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
    typedef YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> MacroGridType;
#else
    typedef SGrid<2, 2> MacroGridType;
#endif
    Stuff::Grid::Providers::Cube<MacroGridType> macro_grid({0., 0.}, {1., 1.}, {3, 3});
    const int levels = num_refinements<LocalGridType>()();
    grid::Multiscale::Glued<MacroGridType, LocalGridType> multiscale_grid(macro_grid, levels);

    size_t num_couplings        = 0;
    const auto& macro_grid_view = multiscale_grid.macro_grid_view();
    for (auto&& macro_entity :
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
                               elements
#else
                               DSC::entityRange
#endif
                                               (macro_grid_view)) {
      const auto entity_index = macro_grid_view.indexSet().index(macro_entity);
      for (auto&& macro_intersection :
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
                                       intersections
#else
                                       DSC::intersectionRange
#endif
                                                             (macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
          const auto macro_neighbor = macro_intersection.outside();
#else
          const auto macro_neighbor_ptr = macro_intersection.outside();
          const auto& macro_neighbor = *macro_neighbor_ptr;
#endif
          const auto neighbor_index = macro_grid_view.indexSet().index(macro_neighbor);
          const auto& coupling = multiscale_grid.coupling(macro_entity, 1, macro_neighbor, levels);
          EXPECT_EQ(4, coupling.size()) << "entity " << entity_index << ", neighbor " << neighbor_index;
          ++num_couplings;
        }
      }
    }
    // because we have 3x3 macro entities -> 12 coupling intersections (twice, since we are looking from each side)
    EXPECT_EQ(num_couplings, 24);
  } // ... works(...)
};  // struct GluedMultiscaleGrid

typedef ::testing::Types< // clang-format off
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
                          YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>
#else
                          SGrid<2, 2>
#endif
#if HAVE_DUNE_ALUGRID
                        , ALUGrid<2, 2, simplex, nonconforming>
#elif HAVE_ALUGRID
                        , ALUGrid<2, 2, simplex, conforming>
#endif
                         > LocalGridTypes; // clang-format on

TYPED_TEST_CASE(GluedMultiscaleGrid, LocalGridTypes);
TYPED_TEST(GluedMultiscaleGrid, works) { this->works(); }
