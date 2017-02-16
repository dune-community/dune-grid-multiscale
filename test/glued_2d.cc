// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include "glued.hh"

typedef ::testing::Types<
                          std::tuple<typename YaspOrSGrid<2>::type, typename YaspOrSGrid<2>::type> // <- works with 0a3346b3ac4119857f74140b5ddb807a867715d1
#if HAVE_DUNE_ALUGRID || HAVE_ALUGRID
//                        , std::tuple<typename YaspOrSGrid<2>::type, Alu2dSimplexType> // <- better with 0a3346b3ac4119857f74140b5ddb807a867715d1
//                        , std::tuple<Alu2dSimplexType,              Alu2dSimplexType> // <- better with 0a3346b3ac4119857f74140b5ddb807a867715d1
#endif // HAVE_DUNE_ALUGRID || HAVE_ALUGRID
                         > GridTypes;

TYPED_TEST_CASE(GluedMultiscaleGrid, GridTypes);
TYPED_TEST(GluedMultiscaleGrid, setup_works) { this->setup(); }
TYPED_TEST(GluedMultiscaleGrid, couplings_are_of_correct_size) { this->couplings_are_of_correct_size(); }
TYPED_TEST(GluedMultiscaleGrid, visualize_is_callable) { this->visualize_is_callable(); }
TYPED_TEST(GluedMultiscaleGrid, coupling_intersections_are_correctly_oriented) { this->coupling_intersections_are_correctly_oriented(); }
