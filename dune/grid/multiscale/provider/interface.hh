// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH

#include <memory>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/gridpart/common/gridpart2gridview.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include "../default.hh"

namespace Dune {
namespace grid {
namespace Multiscale {


template <class GridImp>
class DefaultProvider : public XT::Grid::GridProvider<GridImp>
{
  typedef XT::Grid::GridProvider<GridImp> BaseType;

  //  template <class MSG, XT::Grid::Backends type>
  //  struct ChooseGlobalPartView
  //  {
  //    typedef typename MSG::GlobalGridViewType Type;

  //    static Type create(const MSG& msg)
  //    {
  //      return msg.globalGridView();
  //    }
  //  };

  //  template <class MSG>
  //  struct ChooseGlobalPartView<MSG, XT::Grid::Backends::part>
  //  {
  //    typedef typename MSG::GlobalGridPartType Type;

  //    static Type create(const MSG& msg)
  //    {
  //      return msg.globalGridPart();
  //    }
  //  };

  //  template <class MSG, XT::Grid::Backends type>
  //  struct ChooseLocalPartView
  //  {
  //    typedef typename MSG::LocalGridViewType Type;

  //    static Type create(const MSG& msg, const size_t ss, const bool over)
  //    {
  //      if (over)
  //        DUNE_THROW(NotImplemented,
  //                   "Please add the corresponding method 'localGridView(..., true)' to the multiscale"
  //                       << " grid first!");
  //      return msg.localGridView(ss);
  //    }
  //  };

  //  template <class MSG>
  //  struct ChooseLocalPartView<MSG, XT::Grid::Backends::part>
  //  {
  //    typedef typename MSG::LocalGridPartType Type;

  //    static Type create(const MSG& msg, const size_t ss, const bool over)
  //    {
  //      return msg.localGridPart(ss, over);
  //    }
  //  };

  //  template <class MSG, XT::Grid::Backends type>
  //  struct ChooseBoundaryPartView
  //  {
  //    typedef typename MSG::BoundaryGridViewType Type;

  //    static Type create(const MSG& msg, const size_t ss)
  //    {
  //      return msg.boundaryGridView(ss);
  //    }
  //  };

  //  template <class MSG>
  //  struct ChooseBoundaryPartView<MSG, XT::Grid::Backends::part>
  //  {
  //    typedef typename MSG::BoundaryGridPartType Type;

  //    static Type create(const MSG& msg, const size_t ss)
  //    {
  //      return msg.boundaryGridPart(ss);
  //    }
  //  };

  //  template <class MSG, XT::Grid::Backends type>
  //  struct ChooseCouplingPartView
  //  {
  //    typedef typename MSG::CouplingGridViewType Type;

  //    static Type create(const MSG& msg, const size_t ss, const size_t nn)
  //    {
  //      return msg.couplingGridView(ss, nn);
  //    }
  //  };

  //  template <class MSG>
  //  struct ChooseCouplingPartView<MSG, XT::Grid::Backends::part>
  //  {
  //    typedef typename MSG::CouplingGridPartType Type;

  //    static Type create(const MSG& msg, const size_t ss, const size_t nn)
  //    {
  //      return msg.couplingGridPart(ss, nn);
  //    }
  //  };

  //  template <class MSG, XT::Grid::Layers layer, XT::Grid::Backends part_view>
  //  struct LayerChooser
  //  {
  //    typedef typename XT::Grid::Layer<typename MSG::GridType, layer, part_view>::Type Type;

  //    static Type create(const MSG& ms_grid, const int level)
  //    {
  //      typename MSG::GridType& non_const_grid = const_cast<typename MSG::GridType&>(*(ms_grid.grid()));
  //      return XT::Grid::Layer<typename MSG::GridType, layer, part_view>::create(non_const_grid, level);
  //    }
  //  };

  //  template <class MSG, XT::Grid::Backends part_view>
  //  struct LayerChooser<MSG, XT::Grid::Layers::local, part_view>
  //  {
  //    typedef typename ChooseLocalPartView<MSG, part_view>::Type Type;

  //    static Type create(const MSG& ms_grid, const int level)
  //    {
  //      return ChooseLocalPartView<MSG, part_view>::create(ms_grid, level, false);
  //    }
  //  };

  //  template <class MSG, XT::Grid::Backends part_view>
  //  struct LayerChooser<MSG, XT::Grid::Layers::local_oversampled, part_view>
  //  {
  //    typedef typename ChooseLocalPartView<MSG, part_view>::Type Type;

  //    static Type create(const MSG& ms_grid, const int level)
  //    {
  //      return ChooseLocalPartView<MSG, part_view>::create(ms_grid, level, true);
  //    }
  //  };

public:
  typedef typename BaseType::GridType GridType;
  typedef Default<GridType> MsGridType;

  //  static const unsigned int dimDomain = GridImp::dimension;
  //  typedef typename GridType::ctype DomainFieldType;
  //  typedef FieldVector<DomainFieldType, dimDomain> DomainType;
  //  typedef typename GridType::template Codim<0>::Entity EntityType;

  //  template <XT::Grid::Backends type>
  //  struct Global
  //  {
  //    typedef typename ChooseGlobalPartView<MsGridType, type>::Type Type;
  //  };

  //  template <XT::Grid::Backends type>
  //  struct Local
  //  {
  //    typedef typename ChooseLocalPartView<MsGridType, type>::Type Type;
  //  };

  //  template <XT::Grid::Backends type>
  //  struct Boundary
  //  {
  //    typedef typename ChooseBoundaryPartView<MsGridType, type>::Type Type;
  //  };

  //  template <XT::Grid::Backends type>
  //  struct Coupling
  //  {
  //    typedef typename ChooseCouplingPartView<MsGridType, type>::Type Type;
  //  };

  //  static const std::string static_id()
  //  {
  //    return "grid.multiscale.provider";
  //  }

  DefaultProvider(std::shared_ptr<GridType> grd, std::shared_ptr<MsGridType> ms_grd)
    : BaseType(grd)
    , ms_grid_(ms_grd)
  {
  }

  //  virtual ~ProviderInterface()
  //  {
  //  }

  const std::shared_ptr<const MsGridType>& ms_grid() const
  {
    return ms_grid_;
  }

  //  template <XT::Grid::Backends type>
  //  typename Global<type>::Type global() const
  //  {
  //    return ChooseGlobalPartView<MsGridType, type>::create(*ms_grid());
  //  }

  //  size_t num_subdomains() const
  //  {
  //    return ms_grid()->size();
  //  }

  //  bool oversampling_available() const
  //  {
  //    return ms_grid()->oversampling();
  //  }

  //  template <XT::Grid::Backends type>
  //  typename Local<type>::Type local(const size_t subdomain, const bool oversampling = false) const
  //  {
  //    if (subdomain >= num_subdomains())
  //      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
  //                 "You requsted subdomain number " << subdomain << " of a multiscale grid with only " <<
  //                 num_subdomains()
  //                                                  << " subdomains!\n"
  //                                                  << "Check 'num_subdomains()' first!");
  //    return ChooseLocalPartView<MsGridType, type>::create(*ms_grid(), subdomain, oversampling);
  //  } // ... local(...)

  //  bool is_boundary(const size_t subdomain) const
  //  {
  //    if (subdomain >= num_subdomains())
  //      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
  //                 "You requsted subdomain number " << subdomain << " of a multiscale grid with only " <<
  //                 num_subdomains()
  //                                                  << " subdomains!\n"
  //                                                  << "Check 'num_subdomains()' first!");
  //    return ms_grid()->boundary(subdomain);
  //  } // ... is_boundary(...)

  //  template <XT::Grid::Backends type>
  //  typename Boundary<type>::Type boundary(const size_t subdomain) const
  //  {
  //    if (!is_boundary(subdomain))
  //      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
  //                 "Subdomain " << subdomain << " does not lie at the boundary!"
  //                              << "Check 'is_boundary(subdomain)' first!");
  //    return ChooseBoundaryPartView<MsGridType, type>::create(ms_grid(), subdomain);
  //  } // ... local(...)

  //  std::set<size_t> neighbors(const size_t subdomain) const
  //  {
  //    return ms_grid()->neighborsOf(subdomain);
  //  }

  //  template <XT::Grid::Backends type>
  //  typename Coupling<type>::Type coupling(const size_t subdomain, const size_t neighbor) const
  //  {
  //    if (subdomain >= num_subdomains())
  //      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
  //                 "You requsted subdomain number " << subdomain << " of a multiscale grid with only " <<
  //                 num_subdomains()
  //                                                  << " subdomains!\n"
  //                                                  << "Check 'num_subdomains()' first!");
  //    if (neighbor >= num_subdomains())
  //      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
  //                 "You requsted subdomain number " << neighbor << " of a multiscale grid with only " <<
  //                 num_subdomains()
  //                                                  << " subdomains!\n"
  //                                                  << "Check 'num_subdomains()' first!");
  //    if (ms_grid()->neighborsOf(subdomain).count(neighbor) == 0)
  //      DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
  //                 "Subdomains " << subdomain << " and " << neighbor << " are not neighbors!"
  //                               << "Check 'neighbors(subdomain)' first!");
  //    return ChooseCouplingPartView<MsGridType, type>::create(ms_grid(), subdomain, neighbor);
  //  } // ... local(...)

  //  template <XT::Grid::Layers layer_type, XT::Grid::Backends part_view_type>
  //  typename LayerChooser<MsGridType, layer_type, part_view_type>::Type layer(const int level) const
  //  {
  //    return LayerChooser<MsGridType, layer_type, part_view_type>::create(*ms_grid(), level);
  //  }

  using BaseType::visualize;

  void visualize(const std::string filename, const bool with_coupling) const
  {
    using XT::Common::to_string;
    // vtk writer
    typedef typename MsGridType::GlobalGridPartType GlobalGridPartType;
    const auto& globalGridPart = ms_grid()->globalGridPart();
    typedef Dune::Fem::GridPart2GridView<GlobalGridPartType> GVT;
    GVT globalGridView(globalGridPart);
    Dune::VTKWriter<GVT> vtkwriter(globalGridView);
    // data
    std::map<std::string, std::vector<double>> data;
    data["subdomain"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
    data["global entity id"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
    data["global boundary id"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
    data["local boundary id"] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
    // walk the global grid view
    for (auto it = globalGridView.template begin<0>(); it != globalGridView.template end<0>(); ++it) {
      const auto& entity = *it;
      const auto index = globalGridView.indexSet().index(entity);
      data["subdomain"][index] = ms_grid()->subdomainOf(index);
      data["global entity id"][index] = double(index);
      // compute global boundary id
      data["global boundary id"][index] = 0.0;
      int numberOfBoundarySegments = 0;
      bool isOnBoundary = false;
      for (auto intersectionIt = globalGridView.ibegin(entity); intersectionIt != globalGridView.iend(entity);
           ++intersectionIt) {
        if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
          isOnBoundary = true;
          numberOfBoundarySegments += 1;
          data["global boundary id"][index] += double(intersectionIt->boundarySegmentIndex());
        }
      }
      if (isOnBoundary) {
        data["global boundary id"][index] /= double(numberOfBoundarySegments);
      } // compute global boundary id
    } // walk the global grid view
    // walk the subdomains
    for (unsigned int s = 0; s < ms_grid()->size(); ++s) {
      // walk the local grid view
      const auto localGridView = ms_grid()->localGridPart(s);
      for (auto it = localGridView.template begin<0>(); it != localGridView.template end<0>(); ++it) {
        const auto& entity = *it;
        const unsigned int index = globalGridView.indexSet().index(entity);
        // compute local boundary id
        unsigned int numberOfBoundarySegments = 0u;
        for (auto intersectionIt = localGridView.ibegin(entity); intersectionIt != localGridView.iend(entity);
             ++intersectionIt) {
          if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
            numberOfBoundarySegments += 1u;
            data["local boundary id"][index] += double(intersectionIt->boundarySegmentIndex());
          }
        }
        if (numberOfBoundarySegments > 0)
          data["local boundary id"][index] /= double(numberOfBoundarySegments);
        // visualize coupling
        if (with_coupling) {
          for (auto nn : ms_grid()->neighborsOf(s)) {
            const auto coupling_grid_view = ms_grid()->couplingGridPart(s, nn);
            const std::string coupling_str = "coupling (" + to_string(s) + ", " + to_string(nn) + ")";
            data[coupling_str] = std::vector<double>(globalGridView.indexSet().size(0), 0.0);
            const auto entity_it_end = coupling_grid_view.template end<0>();
            for (auto entity_it = coupling_grid_view.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
              const auto& coupling_entity = *entity_it;
              const size_t global_entity_id = globalGridView.indexSet().index(coupling_entity);
              data[coupling_str][global_entity_id] = 1.0;
              for (auto intersection_it = coupling_grid_view.ibegin(coupling_entity);
                   intersection_it != coupling_grid_view.iend(coupling_entity);
                   ++intersection_it) {
                const auto& intersection = *intersection_it;
                if (intersection.neighbor() && !intersection.boundary()) {
                  const auto neighbor = intersection.outside();
                  const size_t global_neighbor_id = globalGridView.indexSet().index(neighbor);
                  data[coupling_str][global_neighbor_id] = 0.5;
                }
              }
            }
          }
        } // if (with_coupling)
      } // walk the local grid view
    } // walk the subdomains
    if (ms_grid()->oversampling()) {
      // walk the oversampled grid parts
      for (size_t ss = 0; ss < ms_grid()->size(); ++ss) {
        const std::string string_id = "oversampled subdomain " + to_string(ss);
        data[string_id] = std::vector<double>(globalGridView.indexSet().size(0), -1.0);
        typedef typename MsGridType::LocalGridPartType LocalGridPartType;
        const LocalGridPartType oversampledGridPart = ms_grid()->localGridPart(ss, true);
        for (auto it = oversampledGridPart.template begin<0>(); it != oversampledGridPart.template end<0>(); ++it) {
          const auto& entity = *it;
          const auto index = globalGridView.indexSet().index(entity);
          data[string_id][index] = 0.0;
          // compute local boundary id
          unsigned int numberOfBoundarySegments = 0u;
          for (auto intersectionIt = oversampledGridPart.ibegin(entity);
               intersectionIt != oversampledGridPart.iend(entity);
               ++intersectionIt) {
            if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
              numberOfBoundarySegments += 1u;
              data[string_id][index] += double(intersectionIt->boundarySegmentIndex());
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

private:
  std::shared_ptr<const MsGridType> ms_grid_;
}; // class Provider


} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
