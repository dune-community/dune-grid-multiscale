// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include "glued.hh"

namespace Dune {
namespace grid {
namespace Multiscale {
namespace Test {


template <bool anything>
struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
                       YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
                       anything>
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
    return "2d_yaspgrid_yaspgrid";
  }

  static std::set<size_t> num_local_couplings_intersections()
  {
    return {4};
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return false;
  }

  static bool failure_for_higher()
  {
    return false;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
    return {{{0, 0}, 0},
            {{0, 1}, 0},
            {{0, 2}, 0},
            {{1, 0}, 24},
            {{1, 1}, 0},
            {{1, 2}, 0},
            {{2, 0}, 72},
            {{2, 1}, 48},
            {{2, 2}, 0}};
  }
}; // struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, YaspGrid<2,
// EquidistantOffsetCoordinates<double, 2>>, anything>

#if HAVE_DUNE_ALUGRID || HAVE_ALUGRID

template <class Comm, bool anything>
struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
                       ALUGrid<2, 2, simplex, conforming, Comm>,
                       anything>
{
  static int num_coarse_refinements()
  {
#if HAVE_ALUGRID
    return 0;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static int num_local_refinements()
  {
#if HAVE_ALUGRID
    return 3;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static std::string id()
  {
    return "2d_yaspgrid_alugridconforming";
  }

  static std::set<int> num_local_couplings_intersections()
  {
#if HAVE_ALUGRID
    return {2};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return {0};
#endif
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return false;
  }

  static bool failure_for_higher()
  {
    return false;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
#if HAVE_ALUGRID
    return {{{0, 0}, 0},
            {{0, 1}, 0},
            {{0, 2}, 0},
            {{0, 3}, 0},
            {{1, 0}, 0},
            {{1, 1}, 0},
            {{1, 2}, 0},
            {{1, 3}, 0},
            {{2, 0}, 24},
            {{2, 1}, 24},
            {{2, 2}, 0},
            {{2, 3}, 0},
            {{3, 0}, 24},
            {{3, 1}, 24},
            {{3, 2}, 0},
            { {3, 3},
              0 }};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return {};
#endif
  }
}; // struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, ALUGrid<2, 2, simplex, conforming,
// Comm>, anything>

template <class Comm, bool anything>
struct ExpectedResults<ALUGrid<2, 2, simplex, conforming, Comm>, ALUGrid<2, 2, simplex, conforming, Comm>, anything>
{
  static int num_coarse_refinements()
  {
#if HAVE_ALUGRID
    return 1;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static int num_local_refinements()
  {
#if HAVE_ALUGRID
    return 3;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static std::string id()
  {
    return "2d_alugridcoforming_alugridconforming";
  }

  static std::set<size_t> num_local_couplings_intersections()
  {
#if HAVE_ALUGRID
    return {2, 4};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return 0;
#endif
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return false;
  }

  static bool failure_for_higher()
  {
    return false;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
#if HAVE_ALUGRID
    return {{{0, 0}, 0},
            {{0, 1}, 0},
            {{0, 2}, 0},
            {{0, 3}, 0},
            {{1, 0}, 24},
            {{1, 1}, 0},
            {{1, 2}, 0},
            {{1, 3}, 0},
            {{2, 0}, 96},
            {{2, 1}, 72},
            {{2, 2}, 0},
            {{2, 3}, 0},
            {{3, 0}, 144},
            {{3, 1}, 120},
            {{3, 2}, 48},
            { {3, 3},
              0 }};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return {};
#endif
  }
}; // struct ExpectedResults<ALUGrid<2, 2, simplex, conforming, Comm>, ALUGrid<2, 2, simplex, conforming, Comm>,
// anything>

#endif // HAVE_DUNE_ALUGRID || HAVE_ALUGRID
#if HAVE_UG

template <bool anything>
struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, UGGrid<2>, anything>
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
    return "2d_yaspgrid_uggrid";
  }

  static std::set<size_t> num_local_couplings_intersections()
  {
    return {4};
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return false;
  }

  static bool failure_for_higher()
  {
    return false;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
    return {{{0, 0}, 0},
            {{0, 1}, 0},
            {{0, 2}, 0},
            {{1, 0}, 24},
            {{1, 1}, 0},
            {{1, 2}, 0},
            {{2, 0}, 72},
            {{2, 1}, 48},
            {{2, 2}, 0}};
  }
}; // struct ExpectedResults<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, UGGrid<2>, anything>

#endif // HAVE_UG


} // namespace Test
} // namespace Multiscale
} // namespace grid
} // namespace Dune


using namespace Dune;
using namespace Dune::grid::Multiscale::Test;


// clang-format off
typedef ::testing::Types< std::tuple<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
                                     YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>
#if HAVE_ALUGRID
                        , std::tuple<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>,
                                     ALUGrid<2, 2, simplex, conforming>>
                        , std::tuple<ALUGrid<2, 2, simplex, conforming>, ALUGrid<2, 2, simplex, conforming>>
#endif // HAVE_ALUGRID
#if HAVE_UG
                        , std::tuple<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, UGGrid<2>>
#endif
                        > GridTypes; // clang-format on

TYPED_TEST_CASE(GluedMultiscaleGridTest, GridTypes);
TYPED_TEST(GluedMultiscaleGridTest, setup_works)
{
  this->setup();
}
TYPED_TEST(GluedMultiscaleGridTest, visualize_is_callable)
{
  this->visualize_is_callable();
}
TYPED_TEST(GluedMultiscaleGridTest, couplings_are_of_correct_size)
{
  this->couplings_are_of_correct_size();
}
TYPED_TEST(GluedMultiscaleGridTest, intersections_are_correctly_oriented_for_equal_levels)
{
  this->check_intersection_orientation_for_equal_levels();
}
TYPED_TEST(GluedMultiscaleGridTest, intersections_are_correctly_oriented_for_higher_neighbor_levels)
{
  this->check_intersection_orientation_for_higher_neighbor_levels();
}
TYPED_TEST(GluedMultiscaleGridTest,
           __STILL_BROKEN__intersection_orientation_is_wrong_for_lower_or_equal_neighbor_levels)
{
  this->check_intersection_orientation_for_lower_or_equal_neighbor_levels();
}
