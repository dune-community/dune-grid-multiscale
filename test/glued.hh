#ifndef DUNE_GRID_MULTISCALE_TEST_GLUED_HH
#define DUNE_GRID_MULTISCALE_TEST_GLUED_HH

#include <dune/common/version.hh>

// alugrid
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#elif HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
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

#include <dune/stuff/test/gtest/gtest.h>

#include <dune/grid/multiscale/glued.hh>


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


using namespace Dune;


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
      macro_grid_ = DSC::make_unique<Stuff::Grid::Providers::Cube<MacroGridType>>(0., 1., 3, Expectations::num_coarse_refinements());
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    if (!multiscale_grid_)
      multiscale_grid_ = DSC::make_unique<grid::Multiscale::Glued<MacroGridType, LocalGridType>>(*macro_grid_,
                                                                                                 Expectations::num_local_refinements(),
                                                                                                 /*prepare_glues=*/false,
                                                                                                 /*allow_for_broken_orientation_of_coupling_intersections=*/true);
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";
    for (auto&& macro_entity : entity_range(multiscale_grid_->macro_grid_view())) {
      EXPECT_EQ(multiscale_grid_->max_local_level(macro_entity), Expectations::num_local_refinements());
    }
  } // ... setup()

  void couplings_are_of_correct_size()
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";

    const auto& macro_grid_view = multiscale_grid_->macro_grid_view();
    for (auto&& macro_entity : entity_range(macro_grid_view)) {
      const auto entity_index = macro_grid_view.indexSet().index(macro_entity);
      for (auto&& macro_intersection : intersection_range(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor_ptr = macro_intersection.outside();
          const auto& macro_neighbor = *macro_neighbor_ptr;
          const auto neighbor_index = macro_grid_view.indexSet().index(macro_neighbor);
          const auto& coupling = multiscale_grid_->coupling(macro_entity, 1,
                                                            macro_neighbor, Expectations::num_local_refinements(),
                                                            /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          EXPECT_EQ(Expectations::num_local_couplings_intersections().count(coupling.size()), 1)
              << "entity: " << entity_index << "\n"
              << "neighbor: " << neighbor_index << "\n"
              << "expected num_local_couplings_intersections: " << Expectations::num_local_couplings_intersections()
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
      std::cout << "The actual numbers of broken intersections are\n" << count_wrong_intersections_on_all_levels()
                << std::endl;
  } // ... check_intersection_orientation_for_equal_levels(...)

  void check_intersection_orientation_for_higher_neighbor_levels(const bool expect_failure = Expectations::failure_for_higher())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    for (ssize_t entity_level = 0; entity_level <= multiscale_grid_->max_local_level(0); ++entity_level)
      for (ssize_t neighbor_level = entity_level; neighbor_level <= multiscale_grid_->max_local_level(0); ++neighbor_level)
        failure += check_intersections_for_levels(entity_level, neighbor_level, expect_failure);

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n" << count_wrong_intersections_on_all_levels()
                << std::endl;
  } // ... check_intersection_orientation_for_higher_neighbor_levels(...)

  void check_intersection_orientation_for_lower_or_equal_neighbor_levels(const bool expect_failure = Expectations::failure_for_lower_or_equal())
  {
    setup();
    ASSERT_NE(macro_grid_, nullptr) << "This should not happen!";
    ASSERT_NE(multiscale_grid_, nullptr) << "This should not happen!";
    size_t failure = 0;

    for (ssize_t entity_level = 0; entity_level <= multiscale_grid_->max_local_level(0); ++entity_level)
      for (ssize_t neighbor_level = 0; neighbor_level <= entity_level; ++neighbor_level)
        failure += check_intersections_for_levels(entity_level, neighbor_level, expect_failure);

    if (failure)
      std::cout << "The actual numbers of broken intersections are\n" << count_wrong_intersections_on_all_levels()
                << std::endl;
  } // ... check_intersection_orientation_for_lower_or_equal_neighbor_levels(...)

  size_t check_intersections_for_levels(const ssize_t entity_level, const ssize_t neighbor_level, const bool expect_failure)
  {
    size_t failure = 0;
    const auto actual_num_wrongly_oriented_intersections = count_wrong_intersections(entity_level, neighbor_level);
    if (expect_failure) {
      const auto expected_results = Expectations::results();
      const auto search_for_levels_in_expected_results = expected_results.find(std::make_pair(entity_level,
                                                                                              neighbor_level));
      EXPECT_NE(search_for_levels_in_expected_results, expected_results.end())
          << "missing expected results for entity and neighbor level " << std::make_pair(entity_level, neighbor_level)
          << ".\n"
          << "-> PLEASE ADD A RECORD TO\n"
          << "   ExpectedResults<" << DSC::Typename<MacroGridType>::value() << ", "
          << DSC::Typename<LocalGridType>::value() << ">!";
      if (search_for_levels_in_expected_results == expected_results.end()) {
        DUNE_THROW(InvalidStateException,
                   "Cannot use ASSERT_EQ above, so we need to exit this way.\n\n"
                   << "The actual numbers of broken intersections are\n"
                   << count_wrong_intersections_on_all_levels());
      }
      const auto expected_num_wrongly_oriented_intersections = search_for_levels_in_expected_results->second;
      if (expected_num_wrongly_oriented_intersections != actual_num_wrongly_oriented_intersections)
        ++failure;
      EXPECT_EQ(expected_num_wrongly_oriented_intersections, actual_num_wrongly_oriented_intersections)
          << "\nTHIS IS A GOOD THING, AN ACTUAL NUMBER OF FAILURES WHICH IS LOWER THAN THE EXPECTED NUMBER OF "
          << "FAILURES IS AN IMPROVEMENT!\n"
          << "-> BE HAPPY AND UPDATE THE RECORD IN!\n"
          << "   ExpectedResults<" << DSC::Typename<MacroGridType>::value() << ", "
          << DSC::Typename<LocalGridType>::value() << ">!\n"
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
    if(!macro_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- cannot use ASSERT_NE here, non void return
    if(!multiscale_grid_)
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
    if(!macro_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- cannot use ASSERT_NE here, non void return
    if(!multiscale_grid_)
      DUNE_THROW(InvalidStateException, "This should not happen!"); // <- s.a.

    const auto& macro_grid_view = multiscale_grid_->macro_grid_view();
    size_t failures = 0;
    for (auto&& macro_entity : entity_range(macro_grid_view)) {
      for (auto&& macro_intersection : intersection_range(macro_grid_view, macro_entity)) {
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor_ptr = macro_intersection.outside();
          const auto& macro_neighbor = *macro_neighbor_ptr;
//          const auto local_grid_view = multiscale_grid_->local_grid(macro_entity).level_view(entity_level);
          const auto& coupling_glue = multiscale_grid_->coupling(macro_entity, entity_level,
                                                                 macro_neighbor, neighbor_level,
                                                                 /*allow_for_broken_orientation_of_coupling_intersections=*/true);
          failures += Dune::grid::Multiscale::check_for_broken_coupling_intersections(coupling_glue);
//          // walk the coupling
//          const auto coupling_intersection_it_end = coupling_glue.template iend<0>();
//          for (auto coupling_intersection_it = coupling_glue.template ibegin<0>();
//               coupling_intersection_it != coupling_intersection_it_end;
//               ++coupling_intersection_it) {
//            const auto& coupling_intersection = *coupling_intersection_it;
//            const auto coupling_intersection_normal = coupling_intersection.centerUnitOuterNormal();
//            const auto local_entity_ptr = coupling_intersection.inside();
//            const auto& local_entity = *local_entity_ptr;
//            typename std::remove_const<decltype(coupling_intersection_normal)>::type local_intersection_normal(0.);
//            // find the intersection of the local inside entity that corresponds to the coupling intersection
//            size_t found = 0;
//            for (auto&& local_intersection : intersection_range(local_grid_view, local_entity)) {
//              // the coupling intersection may be smaller than the local intersection
//              int corners_inside = 0;
//              for (auto ii : DSC::valueRange(coupling_intersection.geometry().corners()))
//                if (DSG::contains(local_intersection, coupling_intersection.geometry().corner(ii)))
//                  ++corners_inside;
//              if (corners_inside == coupling_intersection.geometry().corners()) {
//                // this is the one
//                ++found;
//                local_intersection_normal = local_intersection.centerUnitOuterNormal();
//              }
//            }
//            EXPECT_EQ(1, found) << "This should not happen!\n"
//                                << "  macro_entity:   " << macro_grid_view.indexSet().index(macro_entity) << "\n"
//                                << "  macro_neighbor: " << macro_grid_view.indexSet().index(macro_neighbor) << "\n"
//                                << "  entity_level:   " << entity_level << "\n"
//                                << "  neighbor_level: " << neighbor_level << "\n"
//                                << "  coupling_intersection: " << coupling_glue.indexSet().index(coupling_intersection);
//            // now the expected normal is local_intersection_normal
//            // and we would like coupling_intersection_normal to point in the same direction
//            // since they have unit length, they should be identical
//            if ((local_intersection_normal - coupling_intersection_normal).infinity_norm() > 1e-15)
//              ++failures;
//          }
        }
      }
    }
    return failures;
  } // ... count_wrong_intersections(...)

std::unique_ptr<Stuff::Grid::Providers::Cube<MacroGridType>> macro_grid_;
std::unique_ptr<grid::Multiscale::Glued<MacroGridType, LocalGridType>> multiscale_grid_;
};  // struct GluedMultiscaleGridTest



#endif // DUNE_GRID_MULTISCALE_TEST_GLUED_HH
