// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH

#include <memory>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/grid/provider/interface.hh>

#include "../default.hh"

namespace Dune {
namespace grid {
namespace Multiscale {

template <class GridImp>
class ProviderInterface : public Stuff::Grid::ConstProviderInterface<GridImp>
{
  typedef Stuff::Grid::ConstProviderInterface<GridImp> BaseType;

  template <class MSG, Stuff::Grid::ChoosePartView type>
  struct ChooseGlobalPartView
  {
    typedef typename MSG::GlobalGridViewType Type;

    static Type create(const MSG& msg) { return msg.globalGridView(); }
  };

  template <class MSG>
  struct ChooseGlobalPartView<MSG, Stuff::Grid::ChoosePartView::part>
  {
    typedef typename MSG::GlobalGridPartType Type;

    static Type create(const MSG& msg) { return msg.globalGridPart(); }
  };

  template <class MSG, Stuff::Grid::ChoosePartView type>
  struct ChooseLocalPartView
  {
    typedef typename MSG::LocalGridViewType Type;

    static Type create(const MSG& msg, const size_t ss, const bool over)
    {
      if (over)
        DUNE_THROW(NotImplemented,
                   "Please add the corresponding method 'localGridView(..., true)' to the multiscale"
                       << " grid first!");
      return msg.localGridView(ss);
    }
  };

  template <class MSG>
  struct ChooseLocalPartView<MSG, Stuff::Grid::ChoosePartView::part>
  {
    typedef typename MSG::LocalGridPartType Type;

    static Type create(const MSG& msg, const size_t ss, const bool over) { return msg.localGridPart(ss, over); }
  };

  template <class MSG, Stuff::Grid::ChoosePartView type>
  struct ChooseBoundaryPartView
  {
    typedef typename MSG::BoundaryGridViewType Type;

    static Type create(const MSG& msg, const size_t ss) { return msg.boundaryGridView(ss); }
  };

  template <class MSG>
  struct ChooseBoundaryPartView<MSG, Stuff::Grid::ChoosePartView::part>
  {
    typedef typename MSG::BoundaryGridPartType Type;

    static Type create(const MSG& msg, const size_t ss) { return msg.boundaryGridPart(ss); }
  };

  template <class MSG, Stuff::Grid::ChoosePartView type>
  struct ChooseCouplingPartView
  {
    typedef typename MSG::CouplingGridViewType Type;

    static Type create(const MSG& msg, const size_t ss, const size_t nn) { return msg.couplingGridView(ss, nn); }
  };

  template <class MSG>
  struct ChooseCouplingPartView<MSG, Stuff::Grid::ChoosePartView::part>
  {
    typedef typename MSG::CouplingGridPartType Type;

    static Type create(const MSG& msg, const size_t ss, const size_t nn) { return msg.couplingGridPart(ss, nn); }
  };

  template <class MSG, Stuff::Grid::ChooseLayer layer, Stuff::Grid::ChoosePartView part_view>
  struct LayerChooser
  {
    typedef typename Stuff::Grid::Layer<typename MSG::GridType, layer, part_view>::Type Type;

    static Type create(const MSG& ms_grid, const int level)
    {
      typename MSG::GridType& non_const_grid = const_cast<typename MSG::GridType&>(*(ms_grid.grid()));
      return Stuff::Grid::Layer<typename MSG::GridType, layer, part_view>::create(non_const_grid, level);
    }
  };

  template <class MSG, Stuff::Grid::ChoosePartView part_view>
  struct LayerChooser<MSG, Stuff::Grid::ChooseLayer::local, part_view>
  {
    typedef typename ChooseLocalPartView<MSG, part_view>::Type Type;

    static Type create(const MSG& ms_grid, const int level)
    {
      return ChooseLocalPartView<MSG, part_view>::create(ms_grid, level, false);
    }
  };

  template <class MSG, Stuff::Grid::ChoosePartView part_view>
  struct LayerChooser<MSG, Stuff::Grid::ChooseLayer::local_oversampled, part_view>
  {
    typedef typename ChooseLocalPartView<MSG, part_view>::Type Type;

    static Type create(const MSG& ms_grid, const int level)
    {
      return ChooseLocalPartView<MSG, part_view>::create(ms_grid, level, true);
    }
  };

public:
  typedef typename BaseType::GridType GridType;
  typedef Default<GridType> MsGridType;

  static const unsigned int dimDomain = GridImp::dimension;
  typedef typename GridType::ctype DomainFieldType;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

  template <Stuff::Grid::ChoosePartView type>
  struct Global
  {
    typedef typename ChooseGlobalPartView<MsGridType, type>::Type Type;
  };

  template <Stuff::Grid::ChoosePartView type>
  struct Local
  {
    typedef typename ChooseLocalPartView<MsGridType, type>::Type Type;
  };

  template <Stuff::Grid::ChoosePartView type>
  struct Boundary
  {
    typedef typename ChooseBoundaryPartView<MsGridType, type>::Type Type;
  };

  template <Stuff::Grid::ChoosePartView type>
  struct Coupling
  {
    typedef typename ChooseCouplingPartView<MsGridType, type>::Type Type;
  };

  static const std::string static_id() { return "grid.multiscale.provider"; }

  virtual ~ProviderInterface() {}

  virtual const std::shared_ptr<const MsGridType>& ms_grid() const = 0;

  template <Stuff::Grid::ChoosePartView type>
  typename Global<type>::Type global() const
  {
    return ChooseGlobalPartView<MsGridType, type>::create(*ms_grid());
  }

  size_t num_subdomains() const { return ms_grid()->size(); }

  bool oversampling_available() const { return ms_grid()->oversampling(); }

  template <Stuff::Grid::ChoosePartView type>
  typename Local<type>::Type local(const size_t subdomain, const bool oversampling = false) const
  {
    if (subdomain >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "You requsted subdomain number " << subdomain << " of a multiscale grid with only " << num_subdomains()
                                                  << " subdomains!\n"
                                                  << "Check 'num_subdomains()' first!");
    return ChooseLocalPartView<MsGridType, type>::create(ms_grid(), subdomain, oversampling);
  } // ... local(...)

  bool is_boundary(const size_t subdomain) const
  {
    if (subdomain >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "You requsted subdomain number " << subdomain << " of a multiscale grid with only " << num_subdomains()
                                                  << " subdomains!\n"
                                                  << "Check 'num_subdomains()' first!");
    return ms_grid()->boundary(subdomain);
  } // ... is_boundary(...)

  template <Stuff::Grid::ChoosePartView type>
  typename Boundary<type>::Type boundary(const size_t subdomain) const
  {
    if (!is_boundary(subdomain))
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Subdomain " << subdomain << " does not lie at the boundary!"
                              << "Check 'is_boundary(subdomain)' first!");
    return ChooseBoundaryPartView<MsGridType, type>::create(ms_grid(), subdomain);
  } // ... local(...)

  std::set<size_t> neighbors(const size_t subdomain) const { return ms_grid()->neighborsOf(subdomain); }

  template <Stuff::Grid::ChoosePartView type>
  typename Coupling<type>::Type coupling(const size_t subdomain, const size_t neighbor) const
  {
    if (subdomain >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "You requsted subdomain number " << subdomain << " of a multiscale grid with only " << num_subdomains()
                                                  << " subdomains!\n"
                                                  << "Check 'num_subdomains()' first!");
    if (neighbor >= num_subdomains())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "You requsted subdomain number " << neighbor << " of a multiscale grid with only " << num_subdomains()
                                                  << " subdomains!\n"
                                                  << "Check 'num_subdomains()' first!");
    if (ms_grid()->neighborsOf(subdomain).count(neighbor) == 0)
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "Subdomains " << subdomain << " and " << neighbor << " are not neighbors!"
                               << "Check 'neighbors(subdomain)' first!");
    return ChooseCouplingPartView<MsGridType, type>::create(ms_grid(), subdomain, neighbor);
  } // ... local(...)

  template <Stuff::Grid::ChooseLayer layer_type, Stuff::Grid::ChoosePartView part_view_type>
  typename LayerChooser<MsGridType, layer_type, part_view_type>::Type layer(const int level) const
  {
    return LayerChooser<MsGridType, layer_type, part_view_type>::create(*ms_grid(), level);
  }

  using BaseType::visualize;

  virtual void visualize(const std::string filename, const bool coupling) const
  {
    // vtk writer
    typedef typename MsGridType::GlobalGridViewType GlobalGridViewType;
    const GlobalGridViewType& globalGridView = ms_grid()->globalGridView();
    Dune::VTKWriter<GlobalGridViewType> vtkwriter(globalGridView);
    // data
    std::map<std::string, std::vector<double>> data;
    data["subdomain"]          = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
    data["global_entity_id"]   = std::vector<double>(globalGridView.indexSet().size(0), -1.);
    data["global_boundary_id"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
    data["local_boundary_id"]  = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
    if (coupling) {
      for (unsigned int s = 0; s < ms_grid()->size(); ++s) {
        data[DSC::toString(s) + " __local_entity_id"]   = std::vector<double>(globalGridView.indexSet().size(0), -1.);
        data[DSC::toString(s) + "__boundary_entity_id"] = std::vector<double>(globalGridView.indexSet().size(0), -1.);
        for (auto nn : ms_grid()->neighborsOf(s)) {
          data[DSC::toString(s) + "_" + DSC::toString(nn) + "__coupling"]
              = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
          data[DSC::toString(s) + "_" + DSC::toString(nn) + "__coupling_entity_id"]
              = std::vector<double>(globalGridView.indexSet().size(0), -1.);
        }
      }
    }
    // walk the global grid view
    for (auto it = globalGridView.template begin<0>(); it != globalGridView.template end<0>(); ++it) {
      const auto& entity = *it;
      const size_t index = globalGridView.indexSet().index(entity);
      data["subdomain"][index]        = ms_grid()->subdomainOf(index);
      data["global_entity_id"][index] = double(index);
      // compute global boundary id
      data["global_boundary_id"][index] = 0.0;
      int numberOfBoundarySegments      = 0;
      bool isOnBoundary = false;
      for (auto intersectionIt = globalGridView.ibegin(entity); intersectionIt != globalGridView.iend(entity); ++intersectionIt) {
        if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
          isOnBoundary = true;
          numberOfBoundarySegments += 1;
          data["global_boundary_id"][index] += double(intersectionIt->boundaryId());
        }
      }
      if (isOnBoundary) {
        data["global_boundary_id"][index] /= double(numberOfBoundarySegments);
      } // compute global boundary id
    }   // walk the global grid view
    // walk the subdomains
    for (unsigned int s = 0; s < ms_grid()->size(); ++s) {
      // walk the local grid view
      const auto localGridView = ms_grid()->localGridPart(s);
      for (auto it = localGridView.template begin<0>(); it != localGridView.template end<0>(); ++it) {
        const auto& entity       = *it;
        const unsigned int index = globalGridView.indexSet().index(entity);
        data[DSC::toString(s) + " __local_entity_id"][index] = localGridView.indexSet().index(entity);
        // compute local boundary id
        unsigned int numberOfBoundarySegments = 0u;
        for (auto intersectionIt = localGridView.ibegin(entity); intersectionIt != localGridView.iend(entity);
             ++intersectionIt) {
          if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
            numberOfBoundarySegments += 1u;
            data["local_boundary_id"][index] += double(intersectionIt->boundaryId());
          }
        }
        if (numberOfBoundarySegments > 0)
          data["local_boundary_id"][index] /= double(numberOfBoundarySegments);
        // visualize coupling
        if (coupling) {
          for (auto nn : ms_grid()->neighborsOf(s)) {
            const auto coupling_grid_view  = ms_grid()->couplingGridPart(s, nn);
            const std::string coupling_str = DSC::toString(s) + "_" + DSC::toString(nn) + "__coupling";
            const auto entity_it_end = coupling_grid_view.template end<0>();
            for (auto entity_it = coupling_grid_view.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
              const auto& coupling_entity          = *entity_it;
              const size_t global_entity_id        = globalGridView.indexSet().index(coupling_entity);
              data[DSC::toString(s) + "__boundary_entity_id"][global_entity_id]
                                                   = coupling_grid_view.indexSet().index(coupling_entity);
              data[coupling_str][global_entity_id] = 1.0;
              for (auto intersection_it = coupling_grid_view.ibegin(coupling_entity);
                   intersection_it != coupling_grid_view.iend(coupling_entity);
                   ++intersection_it) {
                const auto& intersection = *intersection_it;
                if (intersection.neighbor() && !intersection.boundary()) {
                  const auto neighbor_ptr                = intersection.outside();
                  const auto& neighbor                   = *neighbor_ptr;
                  const size_t global_neighbor_id        = globalGridView.indexSet().index(neighbor);
                  data[coupling_str][global_neighbor_id] = 0.5;
                }
              }
            }
          }
          // visualize boundary
          if (ms_grid()->boundary(s)) {
            const auto coupling_grid_view  = ms_grid()->boundaryGridPart(s);
            const auto entity_it_end = coupling_grid_view.template end<0>();
            for (auto entity_it = coupling_grid_view.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
              const auto& entity       = *it;
              const unsigned int index = globalGridView.indexSet().index(entity);
              data[DSC::toString(s) + "__boundary_entity_id"][index] = localGridView.indexSet().index(entity);
            }
          }
        } // if (coupling)
      }   // walk the local grid view
    }     // walk the subdomains
    if (ms_grid()->oversampling()) {
      // walk the oversampled grid parts
      for (size_t ss = 0; ss < ms_grid()->size(); ++ss) {
        const std::string string_id = Stuff::Common::toString(ss) + "__oversampled";
        data[string_id]             = std::vector<double>(globalGridView.indexSet().size(0), -1.0);
        typedef typename MsGridType::LocalGridPartType LocalGridPartType;
        const LocalGridPartType oversampledGridPart = ms_grid()->localGridPart(ss, true);
        for (typename LocalGridPartType::template Codim<0>::IteratorType it = oversampledGridPart.template begin<0>();
             it != oversampledGridPart.template end<0>();
             ++it) {
          const auto& entity = *it;
          const size_t index = globalGridView.indexSet().index(entity);
          data[string_id][index] = 0.0;
          // compute local boundary id
          unsigned int numberOfBoundarySegments = 0u;
          for (auto intersectionIt = oversampledGridPart.ibegin(entity);
               intersectionIt != oversampledGridPart.iend(entity);
               ++intersectionIt) {
            if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
              numberOfBoundarySegments += 1u;
              data[string_id][index]   += double(intersectionIt->boundaryId());
            }
          }
          if (numberOfBoundarySegments > 0) {
            data[string_id][index] /= double(numberOfBoundarySegments);
          } // compute global boundary id
        }
      }
    } // if (ms_grid()->oversampling())

    // write
    for (const auto& data_pair : data)
      vtkwriter.addCellData(data_pair.second, data_pair.first);
    vtkwriter.write(filename, Dune::VTK::appendedraw);
  } // void visualize(const std::string filename = id()) const
};  // class ProviderInterface

} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
