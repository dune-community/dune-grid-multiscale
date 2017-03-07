#ifndef DUNE_GRID_MULTISCALE_TEST_GLUED_HH
#define DUNE_GRID_MULTISCALE_TEST_GLUED_HH

#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/grid/multiscale/glued.hh>


template <class T>
std::string convert_to_initializer_list_str(const std::set<T>& results)
{
  std::stringstream out;
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
  return out.str();
}

template <class L, class R>
std::string convert_to_initializer_list_str(const std::pair<L, R>& results)
{
  std::stringstream out;
  out << "{" << results.first << ", " << results.second << "}";
  return out.str();
}

template <class F, class S>
std::string convert_to_initializer_list_str(const std::map<F, S>& results)
{
  std::stringstream out;
  if (results.size() == 0)
    out << "{}" << std::endl;
  else if (results.size() == 1)
    out << "{{" << convert_to_initializer_list_str(results.begin()->first) << ", " << results.begin()->second << "}}";
  else {
    auto iterator = results.begin();
    out << "{{" << convert_to_initializer_list_str(iterator->first) << ", " << iterator->second << "}";
    ++iterator;
    for (; iterator != results.end(); ++iterator) {
      out << ",\n {" << convert_to_initializer_list_str(iterator->first) << ", " << iterator->second << "}";
    }
    out << "}";
  }
  return out.str();
}


namespace Dune {
namespace grid {
namespace Multiscale {
namespace Test {


template <class M, class L, bool aything = true>
struct ExpectedResults
{
  static_assert(AlwaysFalse<M>::value, "Please add me for this grid!");
};


/// \note assumes that all macro entities contain local grids of same refiment levels
template <class GridTuple>
struct GluedMultiscaleGridTest : public ::testing::Test
{
  typedef typename std::tuple_element<0, GridTuple>::type MacroGridType;
  typedef typename std::tuple_element<1, GridTuple>::type LocalGridType;
  typedef ExpectedResults<MacroGridType, LocalGridType> Expectations;

  void setup()
  {
    if (!macro_grid_)
      macro_grid_ = std::make_unique<XT::Grid::GridProvider<MacroGridType>>(
          XT::Grid::make_cube_grid<MacroGridType>(0., 1., 3, Expectations::num_coarse_refinements()));
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    if (!multiscale_grid_)
      multiscale_grid_ = std::make_unique<grid::Multiscale::Glued<MacroGridType, LocalGridType>>(
          *macro_grid_,
          Expectations::num_local_refinements(),
          /*prepare_glues=*/false,
          /*allow_for_broken_orientation_of_coupling_intersections=*/true);
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";
    for (auto&& macro_entity : Dune::elements(multiscale_grid_->macro_grid_view())) {
      EXPECT_EQ(multiscale_grid_->max_local_level(macro_entity), Expectations::num_local_refinements());
    }
  } // ... setup()

  void couplings_are_of_correct_size()
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";

    const auto& macro_grid_view = multiscale_grid_->macro_grid_view();
    for (auto&& macro_entity : Dune::elements(macro_grid_view)) {
      const auto entity_index = macro_grid_view.indexSet().index(macro_entity);
      for (auto&& macro_intersection : Dune::intersections(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor = macro_intersection.outside();
          const auto neighbor_index = macro_grid_view.indexSet().index(macro_neighbor);
          const auto& coupling =
              multiscale_grid_->coupling(macro_entity,
                                         1,
                                         macro_neighbor,
                                         Expectations::num_local_refinements(),
                                         /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          EXPECT_EQ(Expectations::num_local_couplings_intersections().count(coupling.size()), 1)
              << "entity: " << entity_index << "\n"
              << "neighbor: " << neighbor_index << "\n"
              << "expected num_local_couplings_intersections: "
              << convert_to_initializer_list_str(Expectations::num_local_couplings_intersections())
              << "\nactual num_local_couplings_intersections: " << coupling.size();
        }
      }
    }
  } // ... couplings_are_of_correct_size(...)

  void visualize_is_callable()
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";

    multiscale_grid_->visualize(Expectations::id());
  } // ... visualize_is_callable(...)

  void check_intersection_orientation_for_equal_levels(const bool expect_failure = Expectations::failure_for_equal())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    for (ssize_t level = 0; level <= multiscale_grid_->max_local_level(0); ++level)
      failure += check_intersections_for_levels(level, level, expect_failure);

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n"
                << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()) << std::endl;
  } // ... check_intersection_orientation_for_equal_levels(...)

  void check_intersection_orientation_for_higher_neighbor_levels(
      const bool expect_failure = Expectations::failure_for_higher())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    for (ssize_t entity_level = 0; entity_level <= multiscale_grid_->max_local_level(0); ++entity_level)
      for (ssize_t neighbor_level = entity_level; neighbor_level <= multiscale_grid_->max_local_level(0);
           ++neighbor_level)
        failure += check_intersections_for_levels(entity_level, neighbor_level, expect_failure);

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n"
                << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()) << std::endl;
  } // ... check_intersection_orientation_for_higher_neighbor_levels(...)

  void check_intersection_orientation_for_lower_or_equal_neighbor_levels(
      const bool expect_failure = Expectations::failure_for_lower_or_equal())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    for (ssize_t entity_level = 0; entity_level <= multiscale_grid_->max_local_level(0); ++entity_level)
      for (ssize_t neighbor_level = 0; neighbor_level <= entity_level; ++neighbor_level)
        failure += check_intersections_for_levels(entity_level, neighbor_level, expect_failure);

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n"
                << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()) << std::endl;
  } // ... check_intersection_orientation_for_lower_or_equal_neighbor_levels(...)

  size_t
  check_intersections_for_levels(const ssize_t entity_level, const ssize_t neighbor_level, const bool expect_failure)
  {
    size_t failure = 0;
    const auto actual_num_wrongly_oriented_intersections = count_wrong_intersections(entity_level, neighbor_level);
    if (expect_failure) {
      const auto expected_results = Expectations::results();
      const auto search_for_levels_in_expected_results =
          expected_results.find(std::make_pair(entity_level, neighbor_level));
      EXPECT_NE(search_for_levels_in_expected_results, expected_results.end())
          << "missing expected results for entity and neighbor level "
          << convert_to_initializer_list_str(std::make_pair(entity_level, neighbor_level)) << ".\n"
          << "-> PLEASE ADD A RECORD TO\n"
          << "   ExpectedResults<" << XT::Common::Typename<MacroGridType>::value() << ", "
          << XT::Common::Typename<LocalGridType>::value() << ">!";
      if (search_for_levels_in_expected_results == expected_results.end()) {
        DUNE_THROW(InvalidStateException,
                   "Cannot use ASSERT_EQ above, so we need to exit this way.\n\n"
                       << "The actual numbers of broken intersections are\n"
                       << convert_to_initializer_list_str(count_wrong_intersections_on_all_levels()));
      }
      const auto expected_num_wrongly_oriented_intersections = search_for_levels_in_expected_results->second;
      if (expected_num_wrongly_oriented_intersections != actual_num_wrongly_oriented_intersections)
        ++failure;
      EXPECT_EQ(expected_num_wrongly_oriented_intersections, actual_num_wrongly_oriented_intersections)
          << "\nTHIS IS A GOOD THING, AN ACTUAL NUMBER OF FAILURES WHICH IS LOWER THAN THE EXPECTED NUMBER OF "
          << "FAILURES IS AN IMPROVEMENT!\n"
          << "-> BE HAPPY AND UPDATE THE RECORD IN!\n"
          << "   ExpectedResults<" << XT::Common::Typename<MacroGridType>::value() << ", "
          << XT::Common::Typename<LocalGridType>::value() << ">!\n"
          << "IF THE UPDATED RECORDS DO NOT INDICATE FAILURES ANYMORE, UPDATE THE TESTS!\n";
    } else {
      if (actual_num_wrongly_oriented_intersections != 0)
        ++failure;
      EXPECT_EQ(actual_num_wrongly_oriented_intersections, 0);
    }
    return failure;
  } // ... check_intersections(...)

  std::map<std::pair<size_t, size_t>, size_t> count_wrong_intersections_on_all_levels()
  {
    setup();
    if (!macro_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- cannot use ASSERT_NE here, non void return
    if (!multiscale_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- s.a.

    std::map<std::pair<size_t, size_t>, size_t> results;
    for (ssize_t entity_level = 0; entity_level <= multiscale_grid_->max_local_level(0); ++entity_level)
      for (ssize_t neighbor_level = 0; neighbor_level <= multiscale_grid_->max_local_level(0); ++neighbor_level)
        results[std::make_pair(entity_level, neighbor_level)] = count_wrong_intersections(entity_level, neighbor_level);
    return results;
  }

  size_t count_wrong_intersections(const size_t entity_level, const size_t neighbor_level)
  {
    setup();
    if (!macro_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- cannot use ASSERT_NE here, non void return
    if (!multiscale_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- s.a.

    const auto& macro_grid_view = multiscale_grid_->macro_grid_view();
    size_t failures = 0;
    for (auto&& macro_entity : Dune::elements(macro_grid_view)) {
      for (auto&& macro_intersection : Dune::intersections(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor = macro_intersection.outside();
          //          const auto local_grid_view = multiscale_grid_->local_grid(macro_entity).level_view(entity_level);
          const auto& coupling_glue =
              multiscale_grid_->coupling(macro_entity,
                                         entity_level,
                                         macro_neighbor,
                                         neighbor_level,
                                         /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          failures += Dune::grid::Multiscale::check_for_broken_coupling_intersections(coupling_glue);
        }
      }
    }
    return failures;
  } // ... count_wrong_intersections(...)

  std::unique_ptr<XT::Grid::GridProvider<MacroGridType>> macro_grid_;
  std::unique_ptr<grid::Multiscale::Glued<MacroGridType, LocalGridType>> multiscale_grid_;
}; // struct GluedMultiscaleGridTest


} // namespace Test
} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_TEST_GLUED_HH
