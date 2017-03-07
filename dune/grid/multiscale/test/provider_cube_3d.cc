// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/common/filesystem.hh>

#include "provider_cube.hh"

template <bool anything>
struct ExpectedResults<YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>, anything>
{
  static std::string grid_name()
  {
    return "yasp_3d";
  }

  static std::vector<size_t> local_sizes()
  {
    return std::vector<size_t>(27, 27);
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 19},  {1, 15},  {2, 19},  {3, 15}, {4, 9},   {5, 15},  {6, 19},  {7, 15},  {8, 19},
            {9, 15},  {10, 9},  {11, 15}, {12, 9}, {14, 9},  {15, 15}, {16, 9},  {17, 15}, {18, 19},
            {19, 15}, {20, 19}, {21, 15}, {22, 9}, {23, 15}, {24, 19}, {25, 15}, {26, 19}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 9}, {3, 9}, {9, 9}},
            {{0, 9}, {2, 9}, {4, 9}, {10, 9}},
            {{1, 9}, {5, 9}, {11, 9}},
            {{0, 9}, {4, 9}, {6, 9}, {12, 9}},
            {{1, 9}, {3, 9}, {5, 9}, {7, 9}, {13, 9}},
            {{2, 9}, {4, 9}, {8, 9}, {14, 9}},
            {{3, 9}, {7, 9}, {15, 9}},
            {{4, 9}, {6, 9}, {8, 9}, {16, 9}},
            {{5, 9}, {7, 9}, {17, 9}},
            {{0, 9}, {10, 9}, {12, 9}, {18, 9}},
            {{1, 9}, {9, 9}, {11, 9}, {13, 9}, {19, 9}},
            {{2, 9}, {10, 9}, {14, 9}, {20, 9}},
            {{3, 9}, {9, 9}, {13, 9}, {15, 9}, {21, 9}},
            {{4, 9}, {10, 9}, {12, 9}, {14, 9}, {16, 9}, {22, 9}},
            {{5, 9}, {11, 9}, {13, 9}, {17, 9}, {23, 9}},
            {{6, 9}, {12, 9}, {16, 9}, {24, 9}},
            {{7, 9}, {13, 9}, {15, 9}, {17, 9}, {25, 9}},
            {{8, 9}, {14, 9}, {16, 9}, {26, 9}},
            {{9, 9}, {19, 9}, {21, 9}},
            {{10, 9}, {18, 9}, {20, 9}, {22, 9}},
            {{11, 9}, {19, 9}, {23, 9}},
            {{12, 9}, {18, 9}, {22, 9}, {24, 9}},
            {{13, 9}, {19, 9}, {21, 9}, {23, 9}, {25, 9}},
            {{14, 9}, {20, 9}, {22, 9}, {26, 9}},
            {{15, 9}, {21, 9}, {25, 9}},
            {{16, 9}, {22, 9}, {24, 9}, {26, 9}},
            {{17, 9}, {23, 9}, {25, 9}}};
  }
}; // ... YaspGrid<3, EquidistantOffsetCoordinates<double, 3>> ...

#if HAVE_DUNE_ALUGRID

template <bool anything>
struct ExpectedResults<Dune::ALUGrid<3, 3, cube, nonconforming>, anything>
{
  static std::string grid_name()
  {
    return "alu_3d_cube_nonconforming";
  }

  static std::vector<size_t> local_sizes()
  {
    return std::vector<size_t>(27, 27);
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 19},  {1, 15},  {2, 19},  {3, 15}, {4, 9},   {5, 15},  {6, 19},  {7, 15},  {8, 19},
            {9, 15},  {10, 9},  {11, 15}, {12, 9}, {14, 9},  {15, 15}, {16, 9},  {17, 15}, {18, 19},
            {19, 15}, {20, 19}, {21, 15}, {22, 9}, {23, 15}, {24, 19}, {25, 15}, {26, 19}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 9}, {3, 9}, {9, 9}},
            {{0, 9}, {2, 9}, {4, 9}, {10, 9}},
            {{1, 9}, {5, 9}, {11, 9}},
            {{0, 9}, {4, 9}, {6, 9}, {12, 9}},
            {{1, 9}, {3, 9}, {5, 9}, {7, 9}, {13, 9}},
            {{2, 9}, {4, 9}, {8, 9}, {14, 9}},
            {{3, 9}, {7, 9}, {15, 9}},
            {{4, 9}, {6, 9}, {8, 9}, {16, 9}},
            {{5, 9}, {7, 9}, {17, 9}},
            {{0, 9}, {10, 9}, {12, 9}, {18, 9}},
            {{1, 9}, {9, 9}, {11, 9}, {13, 9}, {19, 9}},
            {{2, 9}, {10, 9}, {14, 9}, {20, 9}},
            {{3, 9}, {9, 9}, {13, 9}, {15, 9}, {21, 9}},
            {{4, 9}, {10, 9}, {12, 9}, {14, 9}, {16, 9}, {22, 9}},
            {{5, 9}, {11, 9}, {13, 9}, {17, 9}, {23, 9}},
            {{6, 9}, {12, 9}, {16, 9}, {24, 9}},
            {{7, 9}, {13, 9}, {15, 9}, {17, 9}, {25, 9}},
            {{8, 9}, {14, 9}, {16, 9}, {26, 9}},
            {{9, 9}, {19, 9}, {21, 9}},
            {{10, 9}, {18, 9}, {20, 9}, {22, 9}},
            {{11, 9}, {19, 9}, {23, 9}},
            {{12, 9}, {18, 9}, {22, 9}, {24, 9}},
            {{13, 9}, {19, 9}, {21, 9}, {23, 9}, {25, 9}},
            {{14, 9}, {20, 9}, {22, 9}, {26, 9}},
            {{15, 9}, {21, 9}, {25, 9}},
            {{16, 9}, {22, 9}, {24, 9}, {26, 9}},
            {{17, 9}, {23, 9}, {25, 9}}};
  }
}; // ... Dune::ALUGrid<3, 3, cube, nonconforming> ...

template <Dune::ALUGridRefinementType ref, bool anything>
struct ExpectedResults<Dune::ALUGrid<3, 3, simplex, ref>, anything>
{
  static std::string grid_name()
  {
    return std::string("alu_3d_simplex_") + (ref == Dune::ALUGridRefinementType::conforming ? "" : "non")
           + "conforming";
  }

  static std::vector<size_t> local_sizes()
  {
    return std::vector<size_t>(27, 162);
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 54},  {1, 36},  {2, 48},  {3, 36},  {4, 18},  {5, 33},  {6, 48},  {7, 33},  {8, 48},
            {9, 36},  {10, 18}, {11, 33}, {12, 18}, {14, 18}, {15, 33}, {16, 18}, {17, 36}, {18, 48},
            {19, 33}, {20, 48}, {21, 33}, {22, 18}, {23, 36}, {24, 48}, {25, 36}, {26, 54}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 18}, {3, 18}, {9, 18}},
            {{0, 18}, {2, 18}, {4, 18}, {10, 18}},
            {{1, 18}, {5, 18}, {11, 18}},
            {{0, 18}, {4, 18}, {6, 18}, {12, 18}},
            {{1, 18}, {3, 18}, {5, 18}, {7, 18}, {13, 18}},
            {{2, 18}, {4, 18}, {8, 18}, {14, 18}},
            {{3, 18}, {7, 18}, {15, 18}},
            {{4, 18}, {6, 18}, {8, 18}, {16, 18}},
            {{5, 18}, {7, 18}, {17, 18}},
            {{0, 18}, {10, 18}, {12, 18}, {18, 18}},
            {{1, 18}, {9, 18}, {11, 18}, {13, 18}, {19, 18}},
            {{2, 18}, {10, 18}, {14, 18}, {20, 18}},
            {{3, 18}, {9, 18}, {13, 18}, {15, 18}, {21, 18}},
            {{4, 18}, {10, 18}, {12, 18}, {14, 18}, {16, 18}, {22, 18}},
            {{5, 18}, {11, 18}, {13, 18}, {17, 18}, {23, 18}},
            {{6, 18}, {12, 18}, {16, 18}, {24, 18}},
            {{7, 18}, {13, 18}, {15, 18}, {17, 18}, {25, 18}},
            {{8, 18}, {14, 18}, {16, 18}, {26, 18}},
            {{9, 18}, {19, 18}, {21, 18}},
            {{10, 18}, {18, 18}, {20, 18}, {22, 18}},
            {{11, 18}, {19, 18}, {23, 18}},
            {{12, 18}, {18, 18}, {22, 18}, {24, 18}},
            {{13, 18}, {19, 18}, {21, 18}, {23, 18}, {25, 18}},
            {{14, 18}, {20, 18}, {22, 18}, {26, 18}},
            {{15, 18}, {21, 18}, {25, 18}},
            {{16, 18}, {22, 18}, {24, 18}, {26, 18}},
            {{17, 18}, {23, 18}, {25, 18}}};
  }
}; // ... Dune::ALUGrid<3, 3, simplex, conforming> ...

#endif // HAVE_DUNE_ALUGRID

// clang-format off
typedef ::testing::Types< YaspGrid<3, EquidistantOffsetCoordinates<double, 3>>
#if HAVE_DUNE_ALUGRID
                        , Dune::ALUGrid<3, 3, cube, nonconforming>
                        , Dune::ALUGrid<3, 3, simplex, conforming>
                        , Dune::ALUGrid<3, 3, simplex, nonconforming>
#endif
                        > GridTypes; // clang-format on

TYPED_TEST_CASE(CubeProviderTest, GridTypes);
TYPED_TEST(CubeProviderTest, setup_works)
{
  this->setup();
}
TYPED_TEST(CubeProviderTest, visualize_is_callable)
{
  this->visualize_is_callable(XT::Common::filename_only(::testing::internal::GetInjectableArgvs().at(0)));
}
TYPED_TEST(CubeProviderTest, entity_to_subdomain_mapping_is_correct)
{
  this->entity_to_subdomain_mapping_is_correct();
}

TYPED_TEST(CubeProviderTest, local_parts_are_of_correct_size)
{
  this->local_parts_are_of_correct_size();
}
TYPED_TEST(CubeProviderTest, local_parts_are_indexed_consecutively)
{
  this->local_parts_are_indexed_consecutively();
}
TYPED_TEST(CubeProviderTest, local_parts_report_correct_boundary_id)
{
  this->local_parts_report_correct_boundary_id();
}

TYPED_TEST(CubeProviderTest, boundary_parts_are_of_correct_size)
{
  this->boundary_parts_are_of_correct_size();
}
TYPED_TEST(CubeProviderTest, boundary_parts_are_indexed_consecutively)
{
  this->boundary_parts_are_indexed_consecutively();
}
TYPED_TEST(CubeProviderTest, boundary_parts_contain_only_boundary_entities_and_intersections)
{
  this->boundary_parts_contain_only_boundary_entities_and_intersections();
}
TYPED_TEST(CubeProviderTest, domain_boundary_is_exactly_covered_by_the_sum_of_local_boundaries)
{
  this->domain_boundary_is_exactly_covered_by_the_sum_of_local_boundaries();
}

TYPED_TEST(CubeProviderTest, coupling_parts_are_of_correct_size)
{
  this->coupling_parts_are_of_correct_size();
}
TYPED_TEST(CubeProviderTest, coupling_parts_are_indexed_consecutively)
{
  this->coupling_parts_are_indexed_consecutively();
}
TYPED_TEST(CubeProviderTest, coupling_parts_contain_only_inner_entities_and_intersections)
{
  this->coupling_parts_contain_only_inner_entities_and_intersections();
}
TYPED_TEST(CubeProviderTest, subdomain_connections_are_exactly_covered_by_couplings)
{
  this->subdomain_connections_are_exactly_covered_by_couplings();
}
