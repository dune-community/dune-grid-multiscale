// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_HH

#include <memory>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"

namespace Dune {
namespace grid {
namespace Multiscale {

template <class GridType>
class MsGridProviders
{
public:
  typedef ProviderInterface<GridType> InterfaceType;

protected:
  template <class GridProviderType>
  static std::unique_ptr<InterfaceType> call_create(const Stuff::Common::Configuration& config)
  {
    if (config.empty())
      return GridProviderType::create();
    else
      return GridProviderType::create(config);
  } // ... call_create(...)

public:
  static std::vector<std::string> available() { return {Providers::Cube<GridType>::static_id()}; } // ... available(...)

  static Stuff::Common::Configuration default_config(const std::string type, const std::string sub_name = "")
  {
    if (type == Providers::Cube<GridType>::static_id())
      return Providers::Cube<GridType>::default_config(sub_name);
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType>
      create(const std::string& type = available()[0],
             const Stuff::Common::Configuration config = Stuff::Common::Configuration())
  {
    if (type == Providers::Cube<GridType>::static_id())
      return call_create<Providers::Cube<GridType>>(config);
    else
      DUNE_THROW(Stuff::Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... create(...)
};  // class MsGridProviders

} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_HH
