#ifndef DUNE_GRID_MULTISCALE_TEST_PROVIDER_CUBE_HH
#define DUNE_GRID_MULTISCALE_TEST_PROVIDER_CUBE_HH

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/grid/grids.hh>

#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/grid/multiscale/provider/cube.hh>

using namespace Dune;

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& results)
{
  if (results.size() == 0)
    out << "{}";
  else if (results.size() == 1)
    out << "{" << *results.begin() << "}";
  else {
    auto iterator = results.begin();
    out << "{" << *iterator;
    ++iterator;
    for (; iterator != results.end(); ++iterator) {
      out << ", " << *iterator;
    }
    out << "}";
  }
  return out;
}

template <class T>
std::ostream& operator<<(std::ostream& out, const std::set<T>& results)
{
  if (results.size() == 0)
    out << "{}";
  else if (results.size() == 1)
    out << "{" << *results.begin() << "}";
  else {
    auto iterator = results.begin();
    out << "{" << *iterator;
    ++iterator;
    for (; iterator != results.end(); ++iterator) {
      out << ", " << *iterator;
    }
    out << "}";
  }
  return out;
}

template <class L, class R>
std::ostream& operator<<(std::ostream& out, const std::pair<L, R>& results)
{
  out << "{" << results.first << ", " << results.second << "}";
  return out;
}

template <class F, class S>
std::ostream& operator<<(std::ostream& out, const std::map<F, S>& results)
{
  if (results.size() == 0)
    out << "{}" << std::endl;
  else if (results.size() == 1)
    out << "{{" << results.begin()->first << ", " << results.begin()->second << "}}";
  else {
    auto iterator = results.begin();
    out << "{{" << iterator->first << ", " << iterator->second << "}";
    ++iterator;
    for (; iterator != results.end(); ++iterator) {
      out << ",\n {" << iterator->first << ", " << iterator->second << "}";
    }
    out << "}";
  }
  return out;
}

template <class G, bool anything = true>
struct ExpectedResults
{
  static_assert(AlwaysFalse<G>::value, "Please add me for this grid!");
};

template <class G>
struct CubeProviderTest : public ::testing::Test
{
  static const constexpr size_t d = G::dimension;
  typedef grid::Multiscale::DefaultProvider<G> ProviderType;
  typedef typename ProviderType::MsGridType MsGridType;

  typedef ExpectedResults<G> Expected;

  void setup()
  {
    // grid definition (see below)
    //   {"lower_left", "upper_right", "num_elements", "num_partitions"},
    //   {"[0 0 0]", "[1 1 1]", "[9 9 9]", "[3 3 3]"});
    // should ensure at least one completely inner subdomain and in each
    // subdomain at least one completely inner entity
    if (!ms_grid_provider_)
      ms_grid_provider_ =
          std::make_shared<ProviderType>(grid::Multiscale::Providers::make_cube_grid<G>(lower_left(),
                                                                                        upper_right(),
                                                                                        num_elements(),
                                                                                        num_refinements(),
                                                                                        overlap_size(),
                                                                                        num_partitions(),
                                                                                        /*oversampling_layers=*/0,
                                                                                        local_boundary_id()));
    if (!ms_grid_provider_w_oversampling_)
      ms_grid_provider_w_oversampling_ =
          std::make_shared<ProviderType>(grid::Multiscale::Providers::make_cube_grid<G>(lower_left(),
                                                                                        upper_right(),
                                                                                        num_elements(),
                                                                                        num_refinements(),
                                                                                        overlap_size(),
                                                                                        num_partitions(),
                                                                                        /*oversampling_layers=*/2,
                                                                                        local_boundary_id()));

    ASSERT_EQ(num_subdomains(), ms_grid_provider_->ms_grid()->size());
    ASSERT_EQ(num_subdomains(), ms_grid_provider_w_oversampling_->ms_grid()->size());
    //    ASSERT_FALSE(ms_grid_provider_->oversampling_available());
    //    ASSERT_TRUE(ms_grid_provider_w_oversampling_->oversampling_available());
  } // ... setup(...)

  void visualize_is_callable(const std::string& filename_prefix)
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    ms_grid_provider_->visualize(/*filename=*/filename_prefix + "__" + Expected::grid_name()
                                     + "__without_oversampling_with_coupling",
                                 /*with_coupling=*/true);
    ms_grid_provider_->visualize(/*filename=*/filename_prefix + "__" + Expected::grid_name()
                                     + "__without_oversampling_without_coupling",
                                 /*with_coupling=*/false);
    //    ms_grid_provider_w_oversampling_->visualize(/*filename=*/filename_prefix
    //    + "_with_oversampling_with_coupling", /*with_coupling=*/true);
    //    ms_grid_provider_w_oversampling_->visualize(/*filename=*/filename_prefix
    //    + "_with_oversampling_without_coupling", /*with_coupling=*/false);
  } // ... visualize_is_callable(...)

  void entity_to_subdomain_mapping_is_correct()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto global_grid_part = ms_grid_provider_->ms_grid()->globalGridPart();
    const auto& entity_to_subdomain_map = *ms_grid_provider_->ms_grid()->entityToSubdomainMap();
    ASSERT_EQ(global_grid_part.indexSet().size(0), entity_to_subdomain_map.size());
    for (auto&& entity : elements(global_grid_part)) {
      const auto entity_index = global_grid_part.indexSet().index(entity);
      auto expected_subdomain = compute_subdomain(entity);
      ASSERT_NE(entity_to_subdomain_map.end(), entity_to_subdomain_map.find(entity_index));
      EXPECT_EQ(expected_subdomain, entity_to_subdomain_map.at(entity_index));
      EXPECT_EQ(expected_subdomain, ms_grid_provider_->ms_grid()->subdomainOf(entity));
      EXPECT_EQ(expected_subdomain, ms_grid_provider_->ms_grid()->subdomainOf(entity_index));
    }
  } // ... entity_to_subdomain_mapping_is_correct(...)

  void local_parts_are_of_correct_size()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto expected_local_sizes = Expected::local_sizes();
    ASSERT_EQ(ms_grid_provider_->ms_grid()->size(), expected_local_sizes.size())
        << "Please update the expected results!"
        << "\n"
        << "expected_local_sizes: " << expected_local_sizes << "\n"
        << "actual local sizes:   " << compute_local_sizes(*ms_grid_provider_);
    size_t total_size = 0;

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      total_size += expected_local_sizes[ss];
      auto local_grid_part = ms_grid_provider_->ms_grid()->localGridPart(ss, false);
      EXPECT_EQ(expected_local_sizes[ss], local_grid_part.indexSet().size(0))
          << "ss: " << ss << "\n"
          << "expected_local_sizes: " << expected_local_sizes << "\n"
          << "actual local sizes:   " << compute_local_sizes(*ms_grid_provider_);
    }
    EXPECT_EQ(ms_grid_provider_->ms_grid()->globalGridPart().indexSet().size(0), total_size);
  } // ... local_parts_are_of_correct_size(...)

  void local_parts_are_indexed_consecutively()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      auto local_grid_part = ms_grid_provider_->ms_grid()->localGridPart(ss, false);
      auto local_indices = compute_local_indices(local_grid_part);
      EXPECT_EQ(1, local_indices.count(0)) << "local indices have to start with 0!\n"
                                           << "ss: " << ss << "\n"
                                           << "local_indices: " << local_indices;
      EXPECT_EQ(Expected::local_sizes()[ss] - 1, *local_indices.rbegin())
          << "local indices have to be numbered consecutively!\n"
          << "ss: " << ss << "\n"
          << "local_indices: " << local_indices;
    }
  } // ... local_parts_are_indexed_consecutively(...)

  void local_parts_report_correct_boundary_id()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto global_grid_part = ms_grid_provider_->ms_grid()->globalGridPart();
    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      auto local_grid_part = ms_grid_provider_->ms_grid()->localGridPart(ss, false);
      for (auto&& entity : elements(local_grid_part)) {
        // we cannot use entity.hasBoundaryIntersections() here!
        for (auto local_intersection_it = local_grid_part.ibegin(entity); // <- DO NOT USE intersection_range HERE!
             local_intersection_it != local_grid_part.iend(entity);
             ++local_intersection_it) {
          const auto& local_intersection = *local_intersection_it;
          const auto local_intersection_index = local_intersection.indexInInside();
          if (local_intersection.boundary() && !local_intersection.neighbor()) {
            size_t global_equals_local = 0;
            // this entity lies on the boundary of the subdomain, so lets see
            // what it looks like globally
            for (auto global_intersection_it =
                     global_grid_part.ibegin(entity); // <- DO NOT USE intersection_range HERE!
                 global_intersection_it != global_grid_part.iend(entity);
                 ++global_intersection_it) {
              const auto& global_intersection = *global_intersection_it;
              if (global_intersection.indexInInside() == local_intersection_index) {
                // this should be the corresponding global intersection
                ++global_equals_local;
                if (global_intersection.boundary() && !global_intersection.neighbor()) {
                  // and this is also on the domain boundary
                  EXPECT_EQ(global_intersection.boundarySegmentIndex(), local_intersection.boundarySegmentIndex());
                } else if (global_intersection.neighbor() && !global_intersection.boundary()) {
                  // and this is an inner intersection globally
                  EXPECT_EQ(local_boundary_id(), local_intersection.boundarySegmentIndex())
                      << "The wrapped intersections of the local grid part "
                         "should report a predefined boundary id on "
                         "subdomain boundaries!";
                } else
                  DUNE_THROW(Dune::InvalidStateException,
                             "This should only happen in parallel runs, which "
                             "are not yet considered here!");
              }
            }
            EXPECT_EQ(1, global_equals_local) << "This must not happen!";
          }
        }
      }
    }
  } // ... local_parts_report_correct_boundary_id(...)

  void boundary_parts_are_of_correct_size()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto expected_boundary_sizes = Expected::boundary_sizes();
    ASSERT_EQ(num_subdomains() - 1, expected_boundary_sizes.size())
        << "The grid is designed to have only 1 subdomain which does not touch "
           "the domain boundary!"
        << "\n"
        << "actual boundary sizes: " << compute_boundary_sizes(*ms_grid_provider_);

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      if (ms_grid_provider_->ms_grid()->boundary(ss)) {
        ASSERT_NE(expected_boundary_sizes.end(), expected_boundary_sizes.find(ss))
            << "This subdomain is not supposed to have a boundary grid part!\n"
            << "ss: " << ss << "\n"
            << "expected_boundary_sizes: " << expected_boundary_sizes << "\n"
            << "actual boundary sizes: " << compute_boundary_sizes(*ms_grid_provider_);
        ASSERT_NE(expected_boundary_sizes.end(), expected_boundary_sizes.find(ss))
            << "Please record the expected results!\n"
            << "expected_boundary_sizes: " << expected_boundary_sizes << "\n"
            << "actual boundary sizes: " << compute_boundary_sizes(*ms_grid_provider_);
        auto boundary_grid_part = ms_grid_provider_->ms_grid()->boundaryGridPart(ss);
        EXPECT_EQ(expected_boundary_sizes[ss], boundary_grid_part.indexSet().size(0))
            << "ss: " << ss << "\n"
            << "expected_boundary_sizes: " << expected_boundary_sizes << "\n"
            << "actual boundary sizes:   " << compute_boundary_sizes(*ms_grid_provider_);
      } else {
        EXPECT_EQ(expected_boundary_sizes.end(), expected_boundary_sizes.find(ss))
            << "This subdomain is supposed to have a boundary grid part!\n"
            << "ss: " << ss << "expected_boundary_sizes: " << expected_boundary_sizes << "\n"
            << "actual boundary sizes: " << compute_boundary_sizes(*ms_grid_provider_);
      }
    }
  } // ... boundary_parts_are_of_correct_size(...)

  void boundary_parts_are_indexed_consecutively()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto expected_boundary_sizes = Expected::boundary_sizes();
    ASSERT_EQ(num_subdomains() - 1, expected_boundary_sizes.size())
        << "The grid is designed to have only 1 subdomain which does not touch "
           "the domain boundary!"
        << "\n"
        << "actual boundary sizes: " << compute_boundary_sizes(*ms_grid_provider_);

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      if (ms_grid_provider_->ms_grid()->boundary(ss)) {
        auto boundary_grid_part = ms_grid_provider_->ms_grid()->boundaryGridPart(ss);
        auto boundary_indices = compute_local_indices(boundary_grid_part);
        EXPECT_EQ(1, boundary_indices.count(0)) << "boundary indices have to start with 0!\n"
                                                << "ss: " << ss << "\n"
                                                << "boundary_indices: " << boundary_indices;
        EXPECT_EQ(expected_boundary_sizes[ss] - 1, *boundary_indices.rbegin())
            << "boundary indices have to be numbered consecutively!\n"
            << "ss: " << ss << "\n"
            << "boundary_indices: " << boundary_indices;
      }
    }
  } // ... boundary_parts_are_indexed_consecutively(...)

  void boundary_parts_contain_only_boundary_entities_and_intersections()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      if (ms_grid_provider_->ms_grid()->boundary(ss)) {
        auto boundary_grid_part = ms_grid_provider_->ms_grid()->boundaryGridPart(ss);
        for (auto&& entity : elements(boundary_grid_part)) {
          EXPECT_TRUE(entity.hasBoundaryIntersections()) << "ss: " << ss;
          for (auto intersection_it = boundary_grid_part.ibegin(entity); // <- DO NOT USE intersection_range HERE!
               intersection_it != boundary_grid_part.iend(entity);
               ++intersection_it) {
            const auto& intersection = *intersection_it;
            EXPECT_TRUE(intersection.boundary() && !intersection.neighbor())
                << "ss: " << ss << "\n"
                << "intersection.boundary(): " << intersection.boundary() << "\n"
                << "intersection.neighbor():" << intersection.neighbor();
          }
        }
      }
    }
  } // ... boundary_parts_contain_only_boundary_entities_and_intersections(...)

  void domain_boundary_is_exactly_covered_by_the_sum_of_local_boundaries()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto global_grid_part = ms_grid_provider_->ms_grid()->globalGridPart();
    auto global_boundary_entities = compute_boundary_indices(global_grid_part);
    std::map<size_t, std::set<size_t>> sum_of_local_boundary_entities;

    size_t min_num_boundary_intersections = 0;
    if (d == 1)
      min_num_boundary_intersections = 2;
    else {
      // we have at least 9 entities in each dimension
      for (size_t ii = 0; ii < d; ++ii)
        min_num_boundary_intersections += 9;
    }
    ASSERT_GE(global_boundary_entities.size(), min_num_boundary_intersections);

    // test that each local boundary entity and intersection is also a global
    // one
    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      if (ms_grid_provider_->ms_grid()->boundary(ss)) {
        auto boundary_grid_part = ms_grid_provider_->ms_grid()->boundaryGridPart(ss);
        auto boundary_entities = compute_boundary_indices(boundary_grid_part);
        for (auto&& entity : elements(boundary_grid_part)) {
          auto boundary_entity_index = boundary_grid_part.indexSet().index(entity);
          ASSERT_NE(boundary_entities.end(), boundary_entities.find(boundary_entity_index))
              << "A boundary grid part should contain only boundary entities!";
          auto boundary_entity_intersections = boundary_entities[boundary_entity_index];
          auto global_entity_index = global_grid_part.indexSet().index(entity);
          ASSERT_NE(global_boundary_entities.end(), global_boundary_entities.find(global_entity_index))
              << "Each entity of a boundary grid part has to be a boundary "
                 "entity of the global grid part!";
          auto global_entity_intersections = global_boundary_entities[global_entity_index];
          EXPECT_EQ(boundary_entity_intersections, global_entity_intersections)
              << "Each boundary intersection of an entity of a boundary grid "
                 "part has to be a boundary intersection of the same entity in "
                 "the global grid part! ";
          EXPECT_EQ(sum_of_local_boundary_entities.find(global_entity_index), sum_of_local_boundary_entities.end())
              << "Each boundary entity of the global grid part should only be "
                 "contained in one boundary grid part!";
          sum_of_local_boundary_entities[global_entity_index] = boundary_entity_intersections;
        }
      }
    }

    // test that there were no more
    EXPECT_EQ(sum_of_local_boundary_entities, global_boundary_entities);
  } // ...
  // domain_boundary_is_exactly_covered_by_the_sum_of_local_boundaries(...)

  void coupling_parts_are_of_correct_size()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto expected_coupling_sizes = Expected::coupling_sizes();
    ASSERT_EQ(expected_coupling_sizes.size(), ms_grid_provider_->ms_grid()->size())
        << "actual coupling sizes: " << compute_coupling_sizes(*ms_grid_provider_->ms_grid());

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      std::set<size_t> expected_neighbors;
      for (const auto& element : expected_coupling_sizes[ss])
        expected_neighbors.insert(element.first);
      auto actual_neighbors = ms_grid_provider_->ms_grid()->neighborsOf(ss);
      ASSERT_EQ(expected_neighbors, actual_neighbors) << "ss: " << ss;
      for (const auto& nn : actual_neighbors) {
        EXPECT_EQ(expected_coupling_sizes[ss][nn],
                  ms_grid_provider_->ms_grid()->couplingGridPart(ss, nn).indexSet().size(0))
            << "ss: " << ss << "\n"
            << "nn: " << nn << "\n"
            << "expected_coupling_sizes: " << expected_coupling_sizes << "\n"
            << "actual coupling sizes: " << compute_coupling_sizes(*ms_grid_provider_->ms_grid());
      }
    }
  } // ... coupling_parts_are_of_correct_size(...)

  void coupling_parts_are_indexed_consecutively()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      for (const auto& nn : ms_grid_provider_->ms_grid()->neighborsOf(ss)) {
        auto coupling_indices = compute_local_indices(ms_grid_provider_->ms_grid()->couplingGridPart(ss, nn));
        EXPECT_EQ(1, coupling_indices.count(0)) << "coupling indices have to start with 0!\n"
                                                << "ss: " << ss << "\n"
                                                << "nn: " << nn << "\n"
                                                << "coupling_indices: " << coupling_indices;
        EXPECT_EQ(Expected::coupling_sizes()[ss][nn] - 1, *coupling_indices.rbegin())
            << "coupling indices have to be numbered consecutively!\n"
            << "ss: " << ss << "\n"
            << "nn: " << nn << "\n"
            << "coupling_indices: " << coupling_indices;
      }
    }
  } // ... coupling_parts_are_indexed_consecutively(...)

  void coupling_parts_contain_only_inner_entities_and_intersections()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      for (const auto& nn : ms_grid_provider_->ms_grid()->neighborsOf(ss)) {
        auto coupling_grid_part = ms_grid_provider_->ms_grid()->couplingGridPart(ss, nn);
        for (auto&& entity : elements(coupling_grid_part)) {
          for (auto intersection_it = coupling_grid_part.ibegin(entity); // <- DO NOT USE intersection_range HERE!
               intersection_it != coupling_grid_part.iend(entity);
               ++intersection_it) {
            const auto& intersection = *intersection_it;
            EXPECT_TRUE(!intersection.boundary() && intersection.neighbor())
                << "ss: " << ss << "\n"
                << "intersection.boundary(): " << intersection.boundary() << "\n"
                << "intersection.neighbor():" << intersection.neighbor();
          }
        }
      }
    }
  } // ... coupling_parts_contain_only_inner_entities_and_intersections(...)

  void subdomain_connections_are_exactly_covered_by_couplings()
  {
    setup();
    ASSERT_NE(nullptr, ms_grid_provider_);
    ASSERT_NE(nullptr, ms_grid_provider_w_oversampling_);

    auto global_grid_part = ms_grid_provider_->ms_grid()->globalGridPart();

    // compute expectations using geometrical information
    auto expected_coupling_indices = compute_coupling_indices(global_grid_part);

    // check information from ms grid
    std::vector<std::map<size_t, std::map<size_t, std::set<size_t>>>> actual_coupling_indices(
        ms_grid_provider_->ms_grid()->size());
    for (size_t ss = 0; ss < ms_grid_provider_->ms_grid()->size(); ++ss) {
      for (const auto& nn : ms_grid_provider_->ms_grid()->neighborsOf(ss)) {
        auto coupling_grid_part = ms_grid_provider_->ms_grid()->couplingGridPart(ss, nn);
        for (auto&& entity : elements(coupling_grid_part)) {
          EXPECT_EQ(compute_subdomain(entity), ss);
          size_t num_coupling_intersections = 0;
          for (auto intersection_it = coupling_grid_part.ibegin(entity); // <- DO NOT USE intersection_range HERE!
               intersection_it != coupling_grid_part.iend(entity);
               ++intersection_it) {
            const auto& intersection = *intersection_it;
            const auto neighbor = intersection.outside(); // has to exist, see test above
            EXPECT_EQ(compute_subdomain(neighbor), nn);
            actual_coupling_indices[ss][nn][global_grid_part.indexSet().index(entity)].insert(
                intersection.indexInInside());
            ++num_coupling_intersections;
          }
          EXPECT_GT(num_coupling_intersections, 0);
        }
        EXPECT_EQ(expected_coupling_indices[ss][nn], actual_coupling_indices[ss][nn]) << "ss: " << ss << "\n"
                                                                                      << "nn: " << nn;
      }
    }
  } // ... subdomain_connections_are_exactly_covered_by_couplings(...)

  // some helper functions

  template <class E>
  size_t compute_subdomain(const E& entity)
  {
    auto center = entity.geometry().center();
    // this code is copied from Providers::Cube
    std::vector<size_t> whichPartition(d, 0);
    for (size_t dd = 0; dd < d; ++dd)
      whichPartition[dd] =
          (std::min((unsigned int)(std::floor(num_partitions()[dd] * ((center[dd] - lower_left()[dd])
                                                                      / (upper_right()[dd] - lower_left()[dd])))),
                    num_partitions()[dd] - 1));
    size_t subdomain = 0;
    if (d == 1)
      subdomain = whichPartition[0];
    else if (d == 2)
      subdomain = whichPartition[0] + whichPartition[1] * num_partitions()[0];
    else if (d == 3)
      subdomain = whichPartition[0] + whichPartition[1] * num_partitions()[0]
                  + whichPartition[2] * num_partitions()[1] * num_partitions()[0];
    else
      DUNE_THROW(Dune::NotImplemented, "for grid dimensions other than 1, 2 or 3!");
    return subdomain;
  } // ... compute_subdomain(...)

  template <class GV>
  std::set<size_t> compute_local_indices(const GV& grid_view)
  {
    std::set<size_t> ret;
    for (auto&& entity : elements(grid_view))
      ret.insert(grid_view.indexSet().index(entity));
    return ret;
  } // ... compute_local_indices(...)

  std::vector<size_t> compute_local_sizes(const ProviderType& provider)
  {
    std::vector<size_t> ret;

    for (size_t ss = 0; ss < provider.ms_grid()->size(); ++ss) {
      auto local_grid_part = provider.ms_grid()->localGridPart(ss, false);
      ret.push_back(local_grid_part.indexSet().size(0));
    }

    return ret;
  } // ... compute_local_sizes(...)

  std::map<size_t, size_t> compute_boundary_sizes(const ProviderType& provider)
  {
    std::map<size_t, size_t> ret;

    for (size_t ss = 0; ss < provider.ms_grid()->size(); ++ss) {
      if (provider.ms_grid()->boundary(ss)) {
        auto boundary_grid_part = provider.ms_grid()->boundaryGridPart(ss);
        ret[ss] = boundary_grid_part.indexSet().size(0);
      }
    }
    return ret;
  } // ... compute_boundary_sizes(...)

  template <class GV>
  std::map<size_t, std::set<size_t>> compute_boundary_indices(const GV& grid_view)
  {
    std::map<size_t, std::set<size_t>> ret;
    for (auto&& entity : elements(grid_view))
      if (entity.hasBoundaryIntersections()) {
        std::set<size_t> intersection_indices;
        for (auto intersection_it = grid_view.ibegin(entity); // <- DO NOT USE intersection_range HERE!
             intersection_it != grid_view.iend(entity);
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          if (intersection.boundary() && !intersection.neighbor())
            intersection_indices.insert(intersection.indexInInside());
        }
        ret[grid_view.indexSet().index(entity)] = intersection_indices;
      }
    return ret;
  } // ... compute_boundary_indices(...)

  std::vector<std::map<size_t, size_t>> compute_coupling_sizes(const MsGridType& ms_grid)
  {
    std::vector<std::map<size_t, size_t>> ret(ms_grid.size());
    for (size_t ss = 0; ss < ms_grid.size(); ++ss)
      for (const auto& nn : ms_grid.neighborsOf(ss))
        ret[ss][nn] = ms_grid.couplingGridPart(ss, nn).indexSet().size(0);
    return ret;
  } // ... compute_coupling_sizes(...)

  template <class GV>
  std::map<size_t, std::map<size_t, std::map<size_t, std::set<size_t>>>>
  compute_coupling_indices(const GV& global_grid_view)
  {
    std::map<size_t, std::map<size_t, std::map<size_t, std::set<size_t>>>> ret;
    for (auto&& entity : elements(global_grid_view)) {
      auto entity_index = global_grid_view.indexSet().index(entity);
      auto ss = compute_subdomain(entity);
      for (auto intersection_it = global_grid_view.ibegin(entity); intersection_it != global_grid_view.iend(entity);
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbor = intersection.outside();
          auto nn = compute_subdomain(neighbor);
          if (nn != ss) { // this is a coupling intersection
            ret[ss][nn][entity_index].insert(intersection.indexInInside());
          }
        }
      }
    }
    return ret;
  } // ... compute_coupling_indices(...)

  // definition of the grid
  //   {"lower_left", "upper_right", "num_elements", "num_partitions"},
  //   {"[0 0 0]", "[1 1 1]", "[9 9 9]", "[3 3 3]"});

  static Dune::FieldVector<double, d> lower_left()
  {
    return Dune::FieldVector<double, d>(0.);
  }

  static Dune::FieldVector<double, d> upper_right()
  {
    return Dune::FieldVector<double, d>(1.);
  }

  static std::array<unsigned int, d> num_elements()
  {
    return XT::Common::make_array<unsigned int, d>(9);
  }

  template <class Grid, bool anything = true>
  struct get_num_refinements
  {
    unsigned int operator()()
    {
      return 0;
    }
  };

  template <class Comm, bool anything>
  struct get_num_refinements<Dune::ALUGrid<2, 2, simplex, conforming, Comm>, anything>
  {
    unsigned int operator()()
    {
      return 1;
    }
  };

  static unsigned int num_refinements()
  {
    return get_num_refinements<G>()();
  }

  static std::array<unsigned int, d> overlap_size()
  {
    return XT::Common::make_array<unsigned int, d>(0);
  }

  static std::array<unsigned int, d> num_partitions()
  {
    return XT::Common::make_array<unsigned int, d>(3);
  }

  static size_t num_subdomains()
  {
    size_t ret = 1;
    for (size_t dd = 0; dd < d; ++dd)
      ret *= num_partitions()[dd];
    return ret;
  }

  static size_t local_boundary_id()
  {
    return std::numeric_limits<size_t>::max() - 17;
  }

  std::shared_ptr<ProviderType> ms_grid_provider_;
  std::shared_ptr<ProviderType> ms_grid_provider_w_oversampling_;
}; // struct CubeProviderTest

#endif // DUNE_GRID_MULTISCALE_TEST_PROVIDER_CUBE_HH
