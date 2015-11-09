// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_GLUED_HH
#define DUNE_GRID_MULTISCALE_GLUED_HH

#include <memory>
#include <map>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/merging/contactmerge.hh>

#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/provider/interface.hh>

namespace Dune {
namespace grid {
namespace Multiscale {

template <class MacroGridImp, class LocalGridImp>
class Glued
{
public:
  typedef MacroGridImp MacroGridType;
  typedef Stuff::Grid::ProviderInterface<MacroGridType> MacroGridProviderType;
  typedef typename MacroGridProviderType::LeafGridViewType MacroGridViewType;
  typedef LocalGridImp LocalGridType;
  typedef Stuff::Grid::ProviderInterface<LocalGridType> LocalGridProviderType;

  typedef typename MacroGridType::ctype ctype;
  static const size_t dimDomain = MacroGridType::dimension;
  static const size_t dimWorld  = MacroGridType::dimensionworld;

  typedef typename MacroGridType::template Codim<0>::Entity MacroEntityType;

private:
  typedef typename LocalGridProviderType::LevelGridViewType LocalViewType;
  typedef GridGlue::Codim1Extractor<LocalViewType> LocalExtractorType;
  typedef Dune::GridGlue::GridGlue<LocalExtractorType, LocalExtractorType> GlueType;

  template <class GridView, class MacroIntersectionType>
  class CouplingFaceDescriptor : public GridGlue::ExtractorPredicate<GridView, 1>
  {
    typedef typename GridView::Traits::template Codim<0>::Entity LocalEntityType;
    typedef typename GridView::ctype ctype;

  public:
    CouplingFaceDescriptor(const MacroIntersectionType& macro_intersection) : macro_intersection_(macro_intersection) {}

    virtual bool contains(const LocalEntityType& element, unsigned int face) const override final
    {
      const auto local_intersection           = element.template subEntity<1>(face);
      const auto& local_intersection_geometry = local_intersection.geometry();
      // Check if all corners of the local intersection lie within the macro intersection.
      for (auto ii : DSC::valueRange(local_intersection_geometry.corners()))
        if (!DSG::contains(macro_intersection_, local_intersection_geometry.corner(ii)))
          return false;
      return true;
    } // ... contains(...)

  private:
    const MacroIntersectionType& macro_intersection_;
  }; // CouplingFaceDescriptor

  template <class GV, class MI>
  static CouplingFaceDescriptor<GV, MI> create_descriptor(const GV& /*gv*/, const MI& mi)
  {
    return CouplingFaceDescriptor<GV, MI>(mi);
  }

public:
  Glued(MacroGridProviderType& macro_grid_provider, const size_t num_local_refinements = 0,
        const bool prepare_glues     = false,
        const ctype& allowed_overlap = 10 * Stuff::Common::FloatCmp::DefaultEpsilon<ctype>::value())
    : macro_grid_(macro_grid_provider)
    , allowed_overlap_(allowed_overlap)
    , macro_leaf_view_(macro_grid_.leaf_view())
    , local_grids_(macro_leaf_view_.indexSet().size(0), nullptr)
    , glues_(macro_leaf_view_.indexSet().size(0))
  {
    setup_local_grids();
    if (num_local_refinements > 0)
      for (auto& local_grid_provider : local_grids_) {
        assert(local_grid_provider);
        local_grid_provider->grid().globalRefine(boost::numeric_cast<int>(num_local_refinements));
      }
    if (prepare_glues)
      setup_glues();
  } // Glued(...)

  const MacroGridViewType& macro_grid_view() const { return macro_leaf_view_; }

  const LocalGridProviderType& local_grid(const MacroEntityType& macro_entity) const
  {
    assert(macro_leaf_view_.indexSet().contains(macro_entity));
    const auto macro_index = macro_leaf_view_.indexSet().index(macro_entity);
    return *(local_grids_.at(macro_index));
  }

  LocalGridProviderType& local_grid(const MacroEntityType& macro_entity)
  {
    assert(macro_leaf_view_.indexSet().contains(macro_entity));
    const auto macro_index = macro_leaf_view_.indexSet().index(macro_entity);
    return *(local_grids_.at(macro_index));
  }

  /**
   * \brief Returns (and creates, if it does not exist) the coupling glue between the local grid view of level
   *        macro_entity_level on macro_entity and the local grid view of level macro_neighbor_level on macro_neighbor.
   * \note  Access is not implemented efficiently. This could be improved by using std::map::find instead of
   *        std::map::operator[].
   */
  const GlueType& coupling(const MacroEntityType& macro_entity, const int macro_entity_level,
                           const MacroEntityType& macro_neighbor, const int macro_neighbor_level) const
  {
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    assert(macro_index_set.contains(macro_entity));
    assert(macro_index_set.contains(macro_neighbor));
    const auto entity_index   = macro_index_set.index(macro_entity);
    const auto neighbor_index = macro_index_set.index(macro_neighbor);
    if (glues_[entity_index][neighbor_index][macro_entity_level][macro_neighbor_level] == nullptr) {
      // find the corresponding macro intersection ...
      const auto macro_intersection_it_end = macro_leaf_view_.iend(macro_entity);
      for (auto macro_intersection_it = macro_leaf_view_.ibegin(macro_entity);
           macro_intersection_it != macro_intersection_it_end;
           ++macro_intersection_it) {
        const auto& macro_intersection = *macro_intersection_it;
        if (macro_intersection.neighbor() && !macro_intersection.boundary())
          if (macro_index_set.index(macro_intersection.outside()) == neighbor_index)
            glues_[entity_index][neighbor_index][macro_entity_level][macro_neighbor_level] =
                create_glue(macro_entity, macro_neighbor, macro_intersection, macro_entity_level, macro_neighbor_level);
      } // ... find the corresponding macro intersection
    }
    return *(glues_[entity_index][neighbor_index][macro_entity_level][macro_neighbor_level]);
  } // ... coupling(...)

private:
  template <class MacroEntityType>
  static std::shared_ptr<LocalGridProviderType> create_grid_of_simplex(const MacroEntityType& macro_entity)
  {
    try {
      GridFactory<LocalGridType> subdomain_factory;
      const auto num_vertices = macro_entity.subEntities(dimDomain);
      std::vector<unsigned int> vertex_ids(num_vertices, 0);
      for (unsigned int local_vertex_id = 0; local_vertex_id < num_vertices; ++local_vertex_id) {
        const auto vertex = macro_entity.template subEntity<dimDomain>(local_vertex_id).geometry().center();
        subdomain_factory.insertVertex(vertex);
        vertex_ids[local_vertex_id] = local_vertex_id;
      }
      subdomain_factory.insertElement(macro_entity.geometry().type(), vertex_ids);
      return std::make_shared<Stuff::Grid::Providers::Default<LocalGridType>>(subdomain_factory.createGrid());
    } catch (GridError& ee) {
      DUNE_THROW(GridError,
                 "It was not possible to create a grid for this simplex with the given GridType!\n\n"
                     << "GridType: "
                     << Stuff::Common::Typename<LocalGridType>::value()
                     << "\n\n"
                     << "This was the original error: "
                     << ee.what());
    }
  } // ... create_grid_of_simplex(...)

  template <class MacroEntityType>
  static std::shared_ptr<LocalGridProviderType> create_grid_of_cube(const MacroEntityType& macro_entity)
  {
    const auto num_vertices = macro_entity.subEntities(dimDomain);
    FieldVector<ctype, dimDomain> lower_left(std::numeric_limits<ctype>::max());
    FieldVector<ctype, dimDomain> upper_right(std::numeric_limits<ctype>::min());
    for (unsigned int local_vertex_id = 0; local_vertex_id < num_vertices; ++local_vertex_id) {
      const auto vertex = macro_entity.template subEntity<dimDomain>(local_vertex_id).geometry().center();
      for (size_t dd = 0; dd < dimDomain; ++dd) {
        lower_left[dd]  = std::min(lower_left[dd], vertex[dd]);
        upper_right[dd] = std::max(upper_right[dd], vertex[dd]);
      }
    }
    return std::make_shared<Stuff::Grid::Providers::Cube<LocalGridType>>(lower_left, upper_right, 1);
  } // ... create_grid_of_cube(...)

  void setup_local_grids()
  {
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      auto macro_entity_index = macro_index_set.index(macro_entity);
      if (macro_entity.type().isSimplex())
        local_grids_[macro_entity_index] = create_grid_of_simplex(macro_entity);
      else if (macro_entity.type().isCube())
        local_grids_[macro_entity_index] = create_grid_of_cube(macro_entity);
      else
        DUNE_THROW(GridError, "Unknown entity.type() encountered: " << macro_entity.type());
    }
  } // ... setup_local_grids()

  template <class MacroIntersectionType>
  std::shared_ptr<GlueType> create_glue(const MacroEntityType& macro_entity, const MacroEntityType& macro_neighbor,
                                        const MacroIntersectionType& macro_intersection, const int local_entity_level,
                                        const int local_neighbor_level) const
  {
    assert(local_entity_level >= 0);
    assert(local_neighbor_level >= 0);
    const auto& local_entity_grid   = local_grid(macro_entity);
    const auto& local_neighbor_grid = local_grid(macro_neighbor);
    assert(local_entity_level <= local_entity_grid.grid().maxLevel());
    assert(local_neighbor_level <= local_neighbor_grid.grid().maxLevel());
    auto local_entity_view   = local_entity_grid.level_view(local_entity_level);
    auto local_neighbor_view = local_neighbor_grid.level_view(local_neighbor_level);
    // create descriptors, these can be discarded after creating the extractors
    auto entity_descriptor   = create_descriptor(local_entity_view, macro_intersection);
    auto neighbor_descriptor = create_descriptor(local_neighbor_view, macro_intersection);
    // create extractors and merger as shared_ptr, so glue will handle memory
    auto entity_extractor   = std::make_shared<LocalExtractorType>(local_entity_view, entity_descriptor);
    auto neighbor_extractor = std::make_shared<LocalExtractorType>(local_neighbor_view, neighbor_descriptor);
    auto contact_merger     = std::make_shared<GridGlue::ContactMerge<dimWorld, ctype>>(allowed_overlap_);
    // create glue
    auto glue = std::make_shared<GlueType>(entity_extractor, neighbor_extractor, contact_merger);
    glue->build();
    if (glue->size() == 0)
      DUNE_THROW(GridError, // clang-format off
                 "Something went wrong, the coupling glue is empty!\n"
                     << "   macro_entity "         << macro_leaf_view_.indexSet().index(macro_entity) << "\n"
                     << "   local_entity_level "   << local_entity_level << "\n"
                     << "   macro_neighbor "       << macro_leaf_view_.indexSet().index(macro_neighbor) << "\n"
                     << "   local_neighbor_level " << local_neighbor_level << "\n"); // clang-format on
    return glue;
  } // ... create_glue(...)

  void setup_glues()
  {
    const auto& macro_index_set = macro_leaf_view_.indexSet();
    for (auto&& macro_entity : elements(macro_leaf_view_)) {
      const auto macro_entity_index = macro_index_set.index(macro_entity);
      auto& entity_glues            = glues_[macro_entity_index];
      // walk the neighbors ...
      const auto macro_intersection_it_end = macro_leaf_view_.iend(macro_entity);
      for (auto macro_intersection_it = macro_leaf_view_.ibegin(macro_entity);
           macro_intersection_it != macro_intersection_it_end;
           ++macro_intersection_it) {
        const auto& macro_intersection = *macro_intersection_it;
        if (macro_intersection.neighbor() && !macro_intersection.boundary()) {
          const auto macro_neighbor       = macro_intersection.outside();
          const auto macro_neighbor_index = macro_index_set.index(macro_neighbor);
          for (auto local_entity_level : DSC::valueRange(local_grids_[macro_entity_index]->grid().maxLevel() + 1))
            for (auto local_neighbor_level : DSC::valueRange(local_grids_[macro_neighbor_index]->grid().maxLevel() + 1))
              entity_glues[macro_neighbor_index][local_entity_level][local_neighbor_level] = create_glue(
                  macro_entity, macro_neighbor, macro_intersection, local_entity_level, local_neighbor_level);
        }
      } // ... walk the neighbors
    }
  } // ... setup_glues(...)

  MacroGridProviderType& macro_grid_;
  const ctype allowed_overlap_;
  MacroGridViewType macro_leaf_view_;
  std::vector<std::shared_ptr<LocalGridProviderType>> local_grids_;
  mutable std::vector<std::map<size_t, std::map<int, std::map<int, std::shared_ptr<GlueType>>>>> glues_;
}; // class Glued

} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_GLUED_HH
