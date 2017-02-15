#ifndef DUNE_GRID_MULTISCALE_TEST_GLUED_HH
#define DUNE_GRID_MULTISCALE_TEST_GLUED_HH

#include <dune/common/version.hh>

// alugrid
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming> Alu2dSimplexType;
#elif HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>  Alu2dSimplexType;
#endif

// yaspgrid or sgrid
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
#include <dune/grid/yaspgrid.hh>
#else
#include <dune/grid/sgrid.hh>
#endif

template <int dim>
struct YaspOrSGrid
{
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
  typedef Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim>> type;
#else
  typedef Dune::SGrid<dim, dim> type;
#endif
};


// elements or DSC::entityRange
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
#include <dune/grid/common/rangegenerators.hh>
#else
#include <dune/stuff/common/ranges.hh>
#endif

template <class... Args>
auto entity_range(Args&& ...args) -> decltype(
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
                                              Dune::elements(std::forward<Args>(args)...)
#else
                                              DSC::entityRange(std::forward<Args>(args)...)
#endif
                                              )
{
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
  return Dune::elements(std::forward<Args>(args)...);
#else
  return DSC::entityRange(std::forward<Args>(args)...);
#endif
}

template <class... Args>
auto intersection_range(Args&& ...args) -> decltype(
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
                                                    Dune::intersections(std::forward<Args>(args)...)
#else
                                                    DSC::intersectionRange(std::forward<Args>(args)...)
#endif
                                                    )
{
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
  return Dune::intersections(std::forward<Args>(args)...);
#else
  return DSC::intersectionRange(std::forward<Args>(args)...);
#endif
}


#include <dune/grid/multiscale/glued.hh>


using namespace Dune;


template <class M, class L, bool aything = true>
struct Switch
{
  static_assert(AlwaysFalse<M>::value, "Please add me for this grid!");
};

template<typename _ctype, bool anything>
struct Switch<SGrid<2, 2, _ctype>, SGrid<2, 2, _ctype>, anything>
{
  static int num_coarse_refinements()
  {
    return 0;
  }

  static int num_local_refinements()
  {
    return 2;
  }

  static std::string id()
  {
    return "2d_sgrid_sgrid";
  }

  static int num_expected_coupling_intersections()
  {
    return 24;
  }
};

#if HAVE_DUNE_ALUGRID || HAVE_ALUGRID

template<typename _ctype, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm, bool anything>
struct Switch<SGrid<2, 2, _ctype>, ALUGrid<2, 2, elType, refineType, Comm>, anything>
{
  static int num_coarse_refinements()
  {
    return 0;
  }

  static int num_local_refinements()
  {
    return 4;
  }

  static std::string id()
  {
    return "2d_sgrid_alugrid";
  }

  static int num_expected_coupling_intersections()
  {
    return 24;
  }
};

template<ALUGridElementType elType, ALUGridRefinementType refineType, class Comm, bool anything>
struct Switch<ALUGrid<2, 2, elType, refineType, Comm>, ALUGrid<2, 2, elType, refineType, Comm>, anything>
{
  static int num_coarse_refinements()
  {
    return 1;
  }

  static int num_local_refinements()
  {
    return 4;
  }

  static std::string id()
  {
    return "2d_alugrid_alugrid";
  }

  static int num_expected_coupling_intersections()
  {
    return 96;
  }
};

#endif // HAVE_DUNE_ALUGRID || HAVE_ALUGRID


template <class GridTuple>
struct GluedMultiscaleGrid : public ::testing::Test
{
  typedef typename std::tuple_element<0, GridTuple>::type MacroGridType;
  typedef typename std::tuple_element<1, GridTuple>::type LocalGridType;
  typedef Switch<MacroGridType, LocalGridType> Switcher;

  void setup()
  {
    if (!macro_grid_)
      macro_grid_ = DSC::make_unique<Stuff::Grid::Providers::Cube<MacroGridType>>(0., 1., 3, Switcher::num_coarse_refinements());
    EXPECT_NE(macro_grid_, nullptr);
    if (!multiscale_grid_)
      multiscale_grid_ = DSC::make_unique<grid::Multiscale::Glued<MacroGridType, LocalGridType>>(*macro_grid_, Switcher::num_local_refinements());
    EXPECT_NE(multiscale_grid_, nullptr);
  } // ... setup()

  void couplings_are_of_correct_size()
  {
    setup();
    if (!macro_grid_ || !multiscale_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!");

    size_t num_couplings        = 0;
    const auto& macro_grid_view = multiscale_grid_->macro_grid_view();
    for (auto&& macro_entity : entity_range(macro_grid_view)) {
      const auto entity_index = macro_grid_view.indexSet().index(macro_entity);
      for (auto&& macro_intersection : intersection_range(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor_ptr = macro_intersection.outside();
          const auto& macro_neighbor = *macro_neighbor_ptr;
          const auto neighbor_index = macro_grid_view.indexSet().index(macro_neighbor);
          const auto& coupling = multiscale_grid_->coupling(macro_entity, 1, macro_neighbor, Switcher::num_local_refinements());
          EXPECT_EQ(4, coupling.size()) << "entity " << entity_index << ", neighbor " << neighbor_index;
          ++num_couplings;
        }
      }
    }
    EXPECT_EQ(num_couplings, Switcher::num_expected_coupling_intersections());
  } // ... couplings_are_of_correct_size(...)

  void visualize_is_callable()
  {
    setup();
    if (!macro_grid_ || !multiscale_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!");

    multiscale_grid_->visualize(Switcher::id());
  } // ... visualize_is_callable(...)

  std::unique_ptr<Stuff::Grid::Providers::Cube<MacroGridType>> macro_grid_;
  std::unique_ptr<grid::Multiscale::Glued<MacroGridType, LocalGridType>> multiscale_grid_;
};  // struct GluedMultiscaleGrid


#endif // DUNE_GRID_MULTISCALE_TEST_GLUED_HH
