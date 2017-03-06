// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/stuff/common/filesystem.hh>

#include "provider_cube.hh"


template <bool anything>
struct ExpectedResults<typename YaspOrSGrid<3>::type, anything>
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
}; // ... typename YaspOrSGrid<3>::type ...

#if HAVE_ALUGRID

template <Dune::ALUGridRefinementType ref, bool anything>
struct ExpectedResults<ALUGrid<3, 3, cube, ref>, anything>
{
  static std::string grid_name()
  {
    if (ref == Dune::ALUGridRefinementType::conforming)
      return "alu_3d_cube_conforming";
    else
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
}; // ... ALUGrid<3, 3, cube, ref> ...

#endif // HAVE_ALUGRID


// clang-format off
typedef ::testing::Types< typename YaspOrSGrid<3>::type
#if HAVE_ALUGRID
                        , ALUGrid<3, 3, cube, conforming>
                        , ALUGrid<3, 3, cube, nonconforming>
#endif
                        >
    GridTypes; // clang-format on

TYPED_TEST_CASE(CubeProviderTest, GridTypes);
TYPED_TEST(CubeProviderTest, setup_works)
{
  this->setup();
}
TYPED_TEST(CubeProviderTest, visualize_is_callable)
{
  this->visualize_is_callable(Stuff::Common::filenameOnly(::testing::internal::GetArgvs()[0]));
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
