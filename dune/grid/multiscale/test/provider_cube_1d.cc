// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/common/filesystem.hh>

#include "provider_cube.hh"


struct Expected1dResults
{
  static std::vector<size_t> local_sizes()
  {
    return {3, 3, 3};
  }

  static std::map<size_t, size_t> boundary_sizes()
  {
    return {{0, 1}, {2, 1}};
  }

  static std::vector<std::map<size_t, size_t>> coupling_sizes()
  {
    return {{{1, 1}}, {{0, 1}, {2, 1}}, {{1, 1}}};
  }
};


template <bool anything>
struct ExpectedResults<YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>, anything> : public Expected1dResults
{
  static std::string grid_name()
  {
    return "yasp_1d";
  }
};

template <bool anything>
struct ExpectedResults<OneDGrid, anything> : public Expected1dResults
{
  static std::string grid_name()
  {
    return "oned_1d";
  }
};

#if HAVE_ALBERTA

template <bool anything>
struct ExpectedResults<AlbertaGrid<1, 1>, anything> : public Expected1dResults
{
  static std::string grid_name()
  {
    return "alberta_1d";
  }
};

#endif // HAVE_ALBERTA


// clang-format off
typedef ::testing::Types< YaspGrid<1, EquidistantOffsetCoordinates<double, 1>>
                        , OneDGrid
#if HAVE_ALBERTA
                        , AlbertaGrid<1, 1>
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
