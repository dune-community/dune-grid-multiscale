// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
#define DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH

#include <memory>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/common/string.hh>

#include "../default.hh"

namespace Dune {
namespace grid {
namespace Multiscale {


template< class GridImp >
class ProviderInterface
{
public:
  typedef GridImp GridType;
  typedef Dune::grid::Multiscale::Default< GridType > MsGridType;
  static const unsigned int dimension = GridType::dimension;
  typedef typename GridType::ctype ctype;
  typedef Dune::FieldVector< ctype, dimension > CoordinateType;

//  typedef ProviderInterface< GridType > ThisType;

  static std::string static_id()
  {
    return "dune.grid.multiscale.provider";
  }

  virtual std::string id() const
  {
    return "dune.grid.multiscale.provider";
  }

  virtual ~ProviderInterface(){}

  virtual std::shared_ptr< const GridType > grid() const = 0;

  virtual std::shared_ptr< const MsGridType > msGrid() const = 0;

  virtual void visualize(const std::string filename = static_id()) const
  {
    // vtk writer
    typedef typename MsGridType::GlobalGridViewType GlobalGridViewType;
    const GlobalGridViewType& globalGridView = *(msGrid()->globalGridView());
    Dune::VTKWriter< GlobalGridViewType > vtkwriter(globalGridView);
    // data
    std::vector< double > subdomainId(globalGridView.indexSet().size(0), 0.0);
    std::vector< double > entityId(globalGridView.indexSet().size(0), 0.0);
    std::vector< double > globalBoundaryId(globalGridView.indexSet().size(0), 0.0);
    std::vector< double > localBoundaryId(globalGridView.indexSet().size(0), 0.0);
    std::vector< std::vector< double > > oversampledSubdomains(msGrid()->size());
    // walk the global grid view
    for (typename GlobalGridViewType::template Codim< 0 >::Iterator it = globalGridView.template begin< 0 >();
         it != globalGridView.template end< 0 >();
         ++it) {
      const typename GlobalGridViewType::template Codim< 0 >::Entity& entity = *it;
      const typename GlobalGridViewType::IndexSet::IndexType index = globalGridView.indexSet().index(entity);
      subdomainId[index] = msGrid()->subdomainOf(index);
      entityId[index] = double(index);
      // compute global boundary id
      globalBoundaryId[index] = 0.0;
      int numberOfBoundarySegments = 0;
      bool isOnBoundary = false;
      for (typename GlobalGridViewType::IntersectionIterator intersectionIt = globalGridView.ibegin(entity);
           intersectionIt != globalGridView.iend(entity);
           ++intersectionIt) {
        if (!intersectionIt->neighbor() && intersectionIt->boundary()) {
          isOnBoundary = true;
          numberOfBoundarySegments += 1;
          globalBoundaryId[index] += double(intersectionIt->boundaryId());
        }
      }
      if (isOnBoundary) {
        globalBoundaryId[index] /= double(numberOfBoundarySegments);
      } // compute global boundary id
    } // walk the global grid view
    // walk the subdomains
    for (unsigned int s = 0; s < msGrid()->size(); ++s) {
      // walk the local grid view
      typedef typename MsGridType::LocalGridViewType LocalGridViewType;
      const LocalGridViewType& localGridView = *(msGrid()->localGridView(s));
      for (typename LocalGridViewType::template Codim< 0 >::Iterator it = localGridView.template begin< 0 >();
           it != localGridView.template end< 0 >();
           ++it) {
        const typename LocalGridViewType::template Codim< 0 >::Entity& entity = *it;
        const typename GlobalGridViewType::IndexSet::IndexType index = globalGridView.indexSet().index(entity);
        // compute local boundary id
        unsigned int numberOfBoundarySegments = 0u;
        for (typename LocalGridViewType::IntersectionIterator intersectionIt = localGridView.ibegin(entity);
             intersectionIt != localGridView.iend(entity);
             ++intersectionIt) {
          if (!intersectionIt->neighbor() && intersectionIt->boundary()){
            numberOfBoundarySegments += 1u;
            localBoundaryId[index] += double(intersectionIt->boundaryId());
          }
        }
        if (numberOfBoundarySegments > 0) {
          localBoundaryId[index] /= double(numberOfBoundarySegments);
        } // compute global boundary id
      } // walk the local grid view
    } // walk the subdomains
    if (msGrid()->oversampling()) {
      // walk the oversampled grid parts
      for (size_t ss = 0; ss < msGrid()->size(); ++ss) {
        oversampledSubdomains[ss] = std::vector< double >(globalGridView.indexSet().size(0), -1.0);
        typedef typename MsGridType::LocalGridPartType LocalGridPartType;
        const LocalGridPartType& oversampledGridPart = *(msGrid()->localGridPart(ss, true));
        for (typename LocalGridPartType::template Codim< 0 >::IteratorType it = oversampledGridPart.template begin< 0 >();
             it != oversampledGridPart.template end< 0 >();
             ++it) {
          const typename LocalGridPartType::template Codim< 0 >::EntityType& entity = *it;
          const typename GlobalGridViewType::IndexSet::IndexType index = globalGridView.indexSet().index(entity);
          oversampledSubdomains[ss][index] = 0.0;
          // compute local boundary id
          unsigned int numberOfBoundarySegments = 0u;
          for (typename LocalGridPartType::IntersectionIteratorType intersectionIt = oversampledGridPart.ibegin(entity);
               intersectionIt != oversampledGridPart.iend(entity);
               ++intersectionIt) {
            if (!intersectionIt->neighbor() && intersectionIt->boundary()){
              numberOfBoundarySegments += 1u;
              oversampledSubdomains[ss][index] += double(intersectionIt->boundaryId());
            }
          }
          if (numberOfBoundarySegments > 0) {
            oversampledSubdomains[ss][index] /= double(numberOfBoundarySegments);
          } // compute global boundary id
        }
      }
    } // if (msGrid()->oversampling())

    // write
    vtkwriter.addCellData(subdomainId, "subdomain Id");
    vtkwriter.addCellData(globalBoundaryId, "global boundary id");
    vtkwriter.addCellData(localBoundaryId, "local boundary id");
    vtkwriter.addCellData(entityId, "entity id");
    if (msGrid()->oversampling()) {
      for (size_t ss = 0; ss < msGrid()->size(); ++ss)
        vtkwriter.addCellData(oversampledSubdomains[ss], "oversampled subdomain " + Dune::Stuff::Common::toString(ss));
    }
    vtkwriter.write(filename, Dune::VTK::ascii);
  } // void visualize(const std::string filename = id()) const
}; // class ProviderInterface


} // namespace Multiscale
} // namespace grid
} // namespace Dune

#endif // DUNE_GRID_MULTISCALE_PROVIDER_INTERFACE_HH
