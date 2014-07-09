// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH

#include <vector>
#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/sgrid.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/grid/multiscale/factory/default.hh>

#include "interface.hh"

namespace Dune {
namespace grid {
namespace Multiscale {
namespace Providers {


template< class GridImp >
class Cube
  : public ProviderInterface< GridImp >
{
  typedef ProviderInterface< GridImp >  BaseType;
  typedef Cube< GridImp >               ThisType;
public:
  typedef typename BaseType::GridType   GridType;
  typedef typename BaseType::MsGridType MsGridType;

  static const unsigned int             dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".cube";
  }

  static Stuff::Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Stuff::Common::ConfigTree config;
    config["lower_left"] = "[0.0 0.0 0.0]";
    config["upper_right"] = "[1.0 1.0 1.0]";
    config["num_elements"] = "[8 8 8]";
    config["num_partitions"] = "[2 2 2]";
    config["oversampling_layers"] = "0";
    if (sub_name.empty())
      return config;
    else {
      Stuff::Common::ConfigTree tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... createSampleDescription(...)

  static std::unique_ptr< ThisType > create(const Stuff::Common::ConfigTree config = default_config(),
                                            const std::string sub_name = static_id())
  {
    // get correct config
    const Stuff::Common::ConfigTree cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Stuff::Common::ConfigTree default_cfg = default_config();
    return Stuff::Common::make_unique< ThisType >(
        cfg.get("lower_left",          default_cfg.get< DomainType >("lower_left"), dimDomain),
        cfg.get("upper_right",         default_cfg.get< DomainType >("upper_right"), dimDomain),
        cfg.get("num_elements",        default_cfg.get< std::vector< unsigned int > >(   "num_elements"), dimDomain),
        cfg.get("num_partitions",      default_cfg.get< std::vector< size_t > >(         "num_partitions"), dimDomain),
        cfg.get("oversampling_layers", default_cfg.get< size_t >(                        "oversampling_layers")));
  } // ... create(...)

  Cube(const DomainType lower_left                    = default_config().get< DomainType >("lower_left"),
       const DomainType upper_right                   = default_config().get< DomainType >("upper_right"),
       const std::vector< unsigned int > num_elements = default_config().get< std::vector< unsigned int >  >("num_elements"),
       const std::vector< size_t > num_partittions    = default_config().get< std::vector< size_t >  >("num_partitions", dimDomain),
       const size_t num_oversampling_layers           = default_config().get< size_t >("oversampling_layers"),
       std::ostream& out                              = DSC_LOG.devnull(),
       const std::string prefix                       = "")
  {
    if (num_partittions.size() < dimDomain)
      DUNE_THROW_COLORFULLY(Dune::RangeError,
                            "num_partittions has to be at least of size " << dimDomain << " (is "
                            << num_partittions.size() << ")!");
#ifndef DUNE_GRID_MULTISCALE_PROVIDER_CUBE_DISABLE_CHECKS
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (num_partittions[ii] > num_elements[ii])
        DUNE_THROW_COLORFULLY(Dune::RangeError,
                              num_partittions[ii] << " = num_partittions[" << ii << "] has to be smaller than "
                              << "num_elements[" << ii << "] = " << num_elements[ii] << "!)");
    }
#endif // DUNE_GRID_MULTISCALE_PROVIDER_CUBE_DISABLE_CHECKS
    typedef Dune::Stuff::Grid::Providers::Cube< GridType > CubeGridProvider;
    auto grd_ptr = CubeGridProvider(lower_left, upper_right, num_elements).grid();
    if (std::is_same< GridType, ALUConformGrid< 2, 2 > >::value
        || std::is_same< GridType, ALUGrid< 2, 2, simplex, conforming > >::value)
      grd_ptr->globalRefine(1);
    grid_ = grd_ptr;
    setup(lower_left, upper_right, num_partittions, num_oversampling_layers, out, prefix);
  }

  Cube(const std::shared_ptr< const GridType > grd,
       const DomainType lower_left                 = default_config().get< DomainType >("lower_left"),
       const DomainType upper_right                = default_config().get< DomainType >("upper_right"),
       const std::vector< size_t > num_partittions = default_config().get< std::vector< size_t >  >("num_partitions", dimDomain),
       const size_t num_oversampling_layers        = default_config().get< size_t >("oversampling_layers"),
       std::ostream& out                           = DSC_LOG.devnull(),
       const std::string prefix                    = "")
    : grid_(grd)
  {
    if (num_partittions.size() < dimDomain)
      DUNE_THROW_COLORFULLY(Dune::RangeError,
                 "num_partittions has to be at least of size " << dimDomain << " (is " << num_partittions.size() << ")!");
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (lower_left[ii] >= upper_right[ii])
        DUNE_THROW_COLORFULLY(Dune::RangeError,
                   lower_left[ii] << " = lower_left[" << ii << "] has to be smaller than upper_right[" << ii
                   << "] = " << upper_right[ii] << "!)");
    }
    setup(lower_left, upper_right, num_partittions, num_oversampling_layers, out, prefix);
  }

  virtual std::shared_ptr< const GridType > grid() const DS_OVERRIDE
  {
    return grid_;
  }

  virtual std::shared_ptr< const MsGridType > ms_grid() const DS_OVERRIDE
  {
    return ms_grid_;
  }

private:
  void setup(const DomainType& lower_left,
             const DomainType& upper_right,
             const std::vector< size_t >& num_partitions,
             const size_t num_oversampling_layers,
             std::ostream& out = DSC_LOG.devnull(), const std::string prefix = "")
  {
    typedef Dune::grid::Multiscale::Factory::Default< GridType > MsGridFactoryType;

    const size_t neighbor_recursion_level = Factory::NeighborRecursionLevel< GridType >::compute();
    // prepare
    MsGridFactoryType factory(grid_);
    factory.prepare();
#ifndef NDEBUG
    // debug output
    out << prefix << static_id()<< ":" << std::endl;
    Stuff::Common::print(lower_left, "lower_left", out, prefix);
    Stuff::Common::print(upper_right, "upper_right", out, prefix);
    Stuff::Common::print(num_partitions, "num_partitions", out, prefix);
#endif // NDEBUG
    // global grid part
    typedef typename MsGridType::GlobalGridPartType GridPartType;
    const auto global_grid_part = factory.globalGridPart();
    // walk the grid
    const auto entity_it_end = global_grid_part->template end< 0 >();
    for (auto entity_it = global_grid_part->template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      // get center of entity
      const auto& entity = *entity_it;
      const auto center = entity.geometry().center();
#ifndef NDEBUG
      const size_t entity_index = global_grid_part->indexSet().index(entity);
      Stuff::Common::print(center, "entity (" + Stuff::Common::toString(entity_index) + ")", out, prefix);
#endif // NDEBUG
      // decide on the subdomain this entity shall belong to
      std::vector< size_t > whichPartition(dimDomain, 0);
      for (size_t dd = 0; dd < dimDomain; ++dd)
        whichPartition[dd] = (std::min((size_t)(std::floor(num_partitions[dd]*((center[dd] - lower_left[dd])/(upper_right[dd] - lower_left[dd])))),
                                       num_partitions[dd] - 1));
      size_t subdomain = 0;
      if (dimDomain == 1)
        subdomain = whichPartition[0];
      else if (dimDomain == 2)
        subdomain = whichPartition[0] + whichPartition[1]*num_partitions[0];
      else if (dimDomain == 3)
        subdomain = whichPartition[0] + whichPartition[1]*num_partitions[0] + whichPartition[2]*num_partitions[1]*num_partitions[0];
      else
        DUNE_THROW_COLORFULLY(Dune::NotImplemented,
                   "ERROR in " << static_id() << ": not implemented for grid dimDomains other than 1, 2 or 3!");
      // add entity to subdomain
      factory.add(entity, subdomain, prefix + "  ", out);
    } // walk the grid
    // finalize
    factory.finalize(num_oversampling_layers, neighbor_recursion_level, prefix + "  ", out);
//    debug << std::flush;
    // be done with it
    ms_grid_ = factory.createMsGrid();
  } // void setup(const Dune::ParameterTree& paramTree)

  std::shared_ptr< const GridType > grid_;
  std::shared_ptr< const MsGridType > ms_grid_;
}; // class Cube


} // namespace Providers
} // namespace Multiscale
} // namespace Grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_CUBE_HH
