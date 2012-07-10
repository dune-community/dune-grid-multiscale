
#ifndef DUNE_GRID_MULTISCALE_MULTIDOMAIN_HH
#define DUNE_GRID_MULTISCALE_MULTIDOMAIN_HH

#ifdef HAVE_DUNE_MULTIDOMAINGRID

// system
#include <set>
#include <sstream>

// dune-common
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

// dune-multidomaingrid
#include <dune/grid/multidomaingrid.hh>

// local
#include "interface.hh"

namespace Dune
{

namespace RB
{

namespace Grid
{

namespace Multiscale
{

template< class GridImp, int maxSubdomains = 50 >
class Multidomain
{
public:
  typedef GridImp HostGridType;

  typedef Multidomain< HostGridType > ThisType;

  typedef Dune::mdgrid::MultiDomainGrid< HostGridType, Dune::mdgrid::FewSubDomainsTraits< HostGridType::dimension, maxSubdomains > > GridType;

  typedef typename GridType::LeafGridView GridViewType;

  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;

  typedef typename Traits::LocalGridViewType LocalGridViewType;

  static const std::string id;

  Multidomain(HostGridType& hostGrid)
    : BaseType(),
      hostGrid_(hostGrid),
      mdGrid_(Dune::shared_ptr< GridType >(new GridType(hostGrid_))),
      finalized_(false),
      size_(0),
      subdomains_()
  {
  }

  GridType& grid()
  {
    return *mdGrid_;
  }

  int size() const
  {
    assert(finalized_);
    return size_;
  }

  bool finalized() const
  {
    return finalized_;
  }

  void prepare()
  {
    mdGrid_->startSubDomainMarking();
  }

  void add(const EntityType& entity, int subdomain)
  {
    // assert
    assert(!finalized_);
    assert(0 <= subdomain);
    // add subdomain id
    subdomains_.insert(subdomain);
    // add entity to subdomain
    mdGrid_->addToSubDomain(subdomain, entity);
  } // void add(const EntityType& entity, int subdomain)

  void finalize()
  {
    // check for consecutive numbering
    size_ = subdomains_.size();
    // set flag
    finalized_ = true;
    for (int i = 0; i < size(); ++i) {
      if (subdomains_.count(i) == 0) {
        finalized_ = false;
        std::stringstream msg;
        msg << "Error in " << id << ": numbering of subdomains not consecutive in the following set (after ordering):" << std::endl;
        typedef std::set< int >::iterator IteratorType;
        for (IteratorType it = subdomains_.begin(); it != subdomains_.end(); ++it) {
          msg << *it << " ";
        }
        DUNE_THROW(Dune::InvalidStateException, msg.str());
      }
    } // check for consecutive numbering
    // finalize multidomaingrid
    mdGrid_->preUpdateSubDomains();
    mdGrid_->updateSubDomains();
    mdGrid_->postUpdateSubDomains();
    // set flag
    finalized_ = true;
  } // void finalize()

  GridViewType gridView()
  {
    return mdGrid_->leafView();
  }

  LocalGridViewType localGridView(int subdomain)
  {
    // assert
    assert(0 <= subdomain);
    assert(subdomain < size());
    assert(finalized_);
    // return
    return mdGrid_->subDomain(subdomain).leafView();
  } // LocalGridViewType localGridView(int subdomain)

private:
  HostGridType& hostGrid_;
  Dune::shared_ptr< GridType > mdGrid_;
  bool finalized_;
  int size_;
  std::set< int > subdomains_;
}; // class Multidomain

template< class GridType, int maxSubdomains >
const std::string Multidomain< GridImp, maxSubdomains >::id = "grid.multiscale.multidomain";

} // namespace Multiscale

} // namespace Grid

} // namespace RB

} // namespace Dune

#endif // HAVE_DUNE_MULTIDOMAINGRID

#endif // DUNE_RB_GRID_MULTISCALE_MULTIDOMAIN_HH
