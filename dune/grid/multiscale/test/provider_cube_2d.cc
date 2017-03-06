// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/stuff/common/filesystem.hh>

#include "provider_cube.hh"


template <bool anything>
struct ExpectedResults<typename YaspOrSGrid<2>::type, anything>
{
  static std::vector<size_t> local_sizes()
  {
    return {9, 9, 9, 9, 9, 9, 9, 9, 9};
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 5}, {1, 3}, {2, 5}, {3, 3}, {5, 3}, {6, 5}, {7, 3}, {8, 5}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 3}, {3, 3}},
            {{0, 3}, {2, 3}, {4, 3}},
            {{1, 3}, {5, 3}},
            {{0, 3}, {4, 3}, {6, 3}},
            {{1, 3}, {3, 3}, {5, 3}, {7, 3}},
            {{2, 3}, {4, 3}, {8, 3}},
            {{3, 3}, {7, 3}},
            {{4, 3}, {6, 3}, {8, 3}},
            {{5, 3}, {7, 3}}};
  }
};

#if HAVE_ALUGRID

template <bool anything>
struct ExpectedResults<ALUGrid<2, 2, simplex, conforming>, anything>
{
  static std::vector<size_t> local_sizes()
  {
    return {36, 36, 36, 36, 36, 36, 36, 36, 36};
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 6}, {1, 3}, {2, 6}, {3, 3}, {5, 3}, {6, 6}, {7, 3}, {8, 6}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 3}, {3, 3}},
            {{0, 3}, {2, 3}, {4, 3}},
            {{1, 3}, {5, 3}},
            {{0, 3}, {4, 3}, {6, 3}},
            {{1, 3}, {3, 3}, {5, 3}, {7, 3}},
            {{2, 3}, {4, 3}, {8, 3}},
            {{3, 3}, {7, 3}},
            {{4, 3}, {6, 3}, {8, 3}},
            {{5, 3}, {7, 3}}};
  }
};

#endif // HAVE_ALUGRID


// clang-format off
typedef ::testing::Types< typename YaspOrSGrid<2>::type
#if HAVE_ALUGRID
                        , ALUGrid<2, 2, simplex, conforming>
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
