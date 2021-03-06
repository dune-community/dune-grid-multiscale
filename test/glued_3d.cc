// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include "glued.hh"

template<typename _ctype, bool anything>
struct ExpectedResults<SGrid<3, 3, _ctype>, SGrid<3, 3, _ctype>, anything>
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
    return "3d_sgrid_sgrid";
  }

  static std::set<int> num_local_couplings_intersections()
  {
    // we expect 16 rectangles, each containing two triangles
    return {32};
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return true;
  }

  static bool failure_for_higher()
  {
    return true;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
    return {{{0, 0}, 216},
            {{0, 1}, 864},
            {{0, 2}, 3456},
            {{1, 0}, 108},
            {{1, 1}, 864},
            {{1, 2}, 3456},
            {{2, 0}, 108},
            {{2, 1}, 108},
            {{2, 2}, 3456}};
  }
}; // struct ExpectedResults<SGrid<3, 3, _ctype>, SGrid<3, 3, _ctype>, anything>

#if HAVE_DUNE_ALUGRID || HAVE_ALUGRID

template<typename _ctype, class Comm, bool anything>
struct ExpectedResults<ALUGrid<3, 3, cube, conforming, Comm>, SGrid<3, 3, _ctype>, anything>
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
    return 2;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static std::string id()
  {
    return "3d_alucubeconforminggrid_sgrid";
  }

  static std::set<int> num_local_couplings_intersections()
  {
#if HAVE_ALUGRID
    return {32};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return true;
  }

  static bool failure_for_higher()
  {
    return true;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
#if HAVE_ALUGRID
    return {{{0, 0}, 216},
            {{0, 1}, 864},
            {{0, 2}, 3456},
            {{1, 0}, 108},
            {{1, 1}, 864},
            {{1, 2}, 3456},
            {{2, 0}, 108},
            {{2, 1}, 108},
            {{2, 2}, 3456}};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return {};
#endif
  }
}; // struct ExpectedResults<SGrid<3, 3, _ctype>, ALUGrid<3, 3, simplex, nonconforming, Comm>, anything>

template<typename _ctype, class Comm, bool anything>
struct ExpectedResults<SGrid<3, 3, _ctype>, ALUGrid<3, 3, cube, nonconforming, Comm>, anything>
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
    return 2;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static std::string id()
  {
    return "3d_alucubenonconforminggrid_sgrid";
  }

  static std::set<int> num_local_couplings_intersections()
  {
#if HAVE_ALUGRID
    return {32};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return true;
  }

  static bool failure_for_higher()
  {
    return true;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
#if HAVE_ALUGRID
    return {{{0, 0}, 216},
            {{0, 1}, 864},
            {{0, 2}, 3456},
            {{1, 0}, 108},
            {{1, 1}, 864},
            {{1, 2}, 3456},
            {{2, 0}, 108},
            {{2, 1}, 108},
            {{2, 2}, 3456}};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return {};
#endif
  }
}; // ExpectedResults<SGrid<3, 3, _ctype>, ALUGrid<3, 3, cube, nonconforming, Comm>, anything>

template<typename _ctype, class Comm, bool anything>
struct ExpectedResults<SGrid<3, 3, _ctype>, ALUGrid<3, 3, simplex, nonconforming, Comm>, anything>
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
    return 2;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static std::string id()
  {
    return "3d_alusimplexnonconforminggrid_sgrid";
  }

  static std::set<int> num_local_couplings_intersections()
  {
#if HAVE_ALUGRID
    return {32};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return true;
  }

  static bool failure_for_higher()
  {
    return true;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
#if HAVE_ALUGRID
    return {{{0, 0}, 216},
            {{0, 1}, 864},
            {{0, 2}, 3456},
            {{1, 0}, 108},
            {{1, 1}, 864},
            {{1, 2}, 3456},
            {{2, 0}, 108},
            {{2, 1}, 108},
            {{2, 2}, 3456}};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return {};
#endif
  }
}; // struct ExpectedResults<SGrid<3, 3, _ctype>, ALUGrid<3, 3, simplex, nonconforming, Comm>, anything>

template<class Comm, bool anything>
struct ExpectedResults<ALUGrid<3, 3, cube, conforming, Comm>, ALUGrid<3, 3, simplex, nonconforming, Comm>, anything>
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
    return 2;
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static std::string id()
  {
    return "3d_alucubeconforminggrid_alusimplexnonconforminggrid";
  }

  static std::set<int> num_local_couplings_intersections()
  {
#if HAVE_ALUGRID
    return {32};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return -1;
#endif
  }

  static bool failure_for_lower_or_equal()
  {
    return true;
  }

  static bool failure_for_equal()
  {
    return true;
  }

  static bool failure_for_higher()
  {
    return true;
  }

  static std::map<std::pair<size_t, size_t>, size_t> results()
  {
#if HAVE_ALUGRID
    return {{{0, 0}, 216},
            {{0, 1}, 864},
            {{0, 2}, 3456},
            {{1, 0}, 108},
            {{1, 1}, 864},
            {{1, 2}, 3456},
            {{2, 0}, 108},
            {{2, 1}, 108},
            {{2, 2}, 3456}};
#else
    DUNE_THROW(InvalidStateException, "Please update these for dune-alugrid!");
    return {};
#endif
  }
}; // struct ExpectedResults<ALUGrid<3, 3, cube, conforming, Comm>, ALUGrid<3, 3, simplex, nonconforming, Comm>, anything>

#endif // HAVE_DUNE_ALUGRID || HAVE_ALUGRID


typedef ::testing::Types<
                          std::tuple<typename YaspOrSGrid<3>::type, typename YaspOrSGrid<3>::type>
#if HAVE_DUNE_ALUGRID || HAVE_ALUGRID
//                      , std::tuple<typename YaspOrSGrid<3>::type, ALUGrid<3, 3, cube, conforming>> // <- knwon to fail completely
                        , std::tuple<typename YaspOrSGrid<3>::type, ALUGrid<3, 3, cube, nonconforming>>
//                      , std::tuple<typename YaspOrSGrid<3>::type, ALUGrid<3, 3, simplex, conforming>> // <- knwon to fail completely
                        , std::tuple<typename YaspOrSGrid<3>::type, ALUGrid<3, 3, simplex, nonconforming>>
                        , std::tuple<ALUGrid<3, 3, cube, conforming>, typename YaspOrSGrid<3>::type>
//                      , std::tuple<ALUGrid<3, 3, simplex, conforming>, ALUGrid<3, 3, simplex, nonconforming>> // <- knwon to fail completely
//                      , std::tuple<ALUGrid<3, 3, simplex, nonconforming>, ALUGrid<3, 3, simplex, nonconforming>> // <- knwon to fail completely
                        , std::tuple<ALUGrid<3, 3, cube, conforming>, ALUGrid<3, 3, simplex, nonconforming>>
#endif // HAVE_DUNE_ALUGRID || HAVE_ALUGRID
                         > GridTypes;

TYPED_TEST_CASE(GluedMultiscaleGridTest, GridTypes);
TYPED_TEST(GluedMultiscaleGridTest, setup_works) { this->setup(); }
TYPED_TEST(GluedMultiscaleGridTest, visualize_is_callable) { this->visualize_is_callable(); }
TYPED_TEST(GluedMultiscaleGridTest, couplings_are_of_correct_size) { this->couplings_are_of_correct_size(); }
TYPED_TEST(GluedMultiscaleGridTest, __STILL_BROKEN__intersections_are_correctly_oriented_for_equal_levels)
{
  this->check_intersection_orientation_for_equal_levels();
}
TYPED_TEST(GluedMultiscaleGridTest, __STILL_BROKEN__intersections_are_correctly_oriented_for_higher_neighbor_levels)
{
  this->check_intersection_orientation_for_higher_neighbor_levels();
}
TYPED_TEST(GluedMultiscaleGridTest, __STILL_BROKEN__intersection_orientation_is_wrong_for_lower_or_equal_neighbor_levels)
{
  this->check_intersection_orientation_for_lower_or_equal_neighbor_levels();
}
