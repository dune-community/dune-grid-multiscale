// This file is part of the dune-grid-multiscale project:
//   http://users.dune-project.org/projects/dune-grid-multiscale
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hh> // <- has to come first

#include <dune/grid/sgrid.hh>

#include <dune/grid/multiscale/provider/cube.hh>


using namespace Dune;
typedef SGrid< 2, 2 > GridType;


template< class GridPartType >
void measureTiming(const GridPartType& gridPart, Dune::Stuff::Common::LogStream& out, const std::string name = "")
{
  out << "  walking " << name << " grid part... " << std::flush;
  Dune::Timer timer;
  unsigned int elements = 0;
  for (typename GridPartType::template Codim< 0 >::IteratorType it = gridPart.template begin< 0 >();
       it != gridPart.template end< 0 >();
       ++it)
    ++elements;
  out << " done (has " << elements << " elements, took " << timer.elapsed() << " sek)" << std::endl;
} // void measureTiming()


template< class GlobalGridPartType, class LocalGridPartType, class OutStreamType >
struct Inspect
{
  template< int dim, int codim >
  struct Codim
  {
    static void entities(const GlobalGridPartType& globalGridPart, const LocalGridPartType& localGridPart, const std::string prefix, OutStreamType& out)
    {
      // walk all codim entitites
      for (auto entityIt = localGridPart.template begin< codim >(); entityIt != localGridPart.template end< codim >(); ++entityIt) {
        const auto& entity = *entityIt;
        const auto& geometryType = entity.type();
        const auto globalIndex = globalGridPart.indexSet().index(entity);
        const auto localIndex = localGridPart.indexSet().index(entity);
        out << prefix << geometryType << ", global index " << globalIndex << ", local index " << localIndex << std::endl;
        out << prefix << "  corners: ";
        auto geometry = entity.geometry();
        for (int i = 0; i < (geometry.corners() - 1); ++i)
          out << "(" << geometry.corner(i) << "), ";
        out << "(" << geometry.corner(geometry.corners() - 1) << ")" << std::endl;
      } // walk all codim entitites
      // increase Codim
      Inspect< GlobalGridPartType, LocalGridPartType, OutStreamType >::Codim< dim, codim + 1 >::entities(globalGridPart, localGridPart, prefix, out);
    }
  }; // struct Codim
}; // struct Inspect


template< class GlobalGridPartType, class LocalGridPartType, class OutStreamType >
template< int codim >
struct Inspect< GlobalGridPartType, LocalGridPartType, OutStreamType >::Codim< codim, codim >
{
  static void entities(const GlobalGridPartType& globalGridPart, const LocalGridPartType& localGridPart, const std::string prefix, OutStreamType& out)
  {
    // walk all codim entitites
    for (auto entityIt = localGridPart.template begin< codim >(); entityIt != localGridPart.template end< codim >(); ++entityIt) {
      const auto& entity = *entityIt;
      const auto& geometryType = entity.type();
      const auto globalIndex = globalGridPart.indexSet().index(entity);
      const auto localIndex = localGridPart.indexSet().index(entity);
      out << prefix << geometryType << ", global index " << globalIndex << ", local index " << localIndex << std::endl;
      out << prefix << "  corners: (" << entity.geometry().center() << ")" << std::endl;
    } // walk all codim entitites
  }
};


class CubeProvider
  : public ::testing::Test
{
protected:
  typedef grid::Multiscale::Providers::Cube< GridType > ProviderType;
  typedef typename ProviderType::MsGridType MsGridType;

  CubeProvider()
    : ms_grid_(ProviderType::create())
  {}

  std::shared_ptr< MsGridType > ms_grid_;
}; // class CubeProvider


TEST_F(CubeProvider, global_gridpart)
{
  DSC_LOG_INFO << "inspecting global grid part:" << std::endl;
  const auto& globalGridPart = *(ms_grid_->globalGridPart());
  Inspect< MsGridType::GlobalGridPartType, MsGridType::GlobalGridPartType, Dune::Stuff::Common::LogStream >
      ::Codim< MultiscaleGridProviderType::dim, 0 >
      ::entities(globalGridPart, globalGridPart, "  ", DSC_LOG_DEBUG);
}

TEST_F(CubeProvider, local_grid_parts)
{
  DSC_LOG_INFO << "inspecting local grid parts:" << std::endl;
  for (unsigned int subdomain = 0; subdomain < ms_grid_->size(); ++subdomain) {
    DSC_LOG_INFO << "subdomain " << subdomain << std::endl;
    const auto& localGridPart = *(ms_grid_->localGridPart(subdomain));
    Inspect< MsGridType::GlobalGridPartType, MsGridType::LocalGridPartType, Dune::Stuff::Common::LogStream >
        ::Codim< MultiscaleGridProviderType::dim, 0 >
        ::entities(globalGridPart, localGridPart, "  ", DSC_LOG_DEBUG);
  }
}

TEST_F(CubeProvider, boundary_grid_parts)
{
    DSC_LOG_INFO << "inspecting boundary grid parts:" << std::endl;
    for (unsigned int subdomain = 0; subdomain < ms_grid_->size(); ++subdomain) {
      if (ms_grid_->boundary(subdomain)) {
        DSC_LOG_INFO << "subdomain " << subdomain << std::endl;
        const auto& boundaryGridPart = *(ms_grid_->boundaryGridPart(subdomain));
        Inspect< MsGridType::GlobalGridPartType, MsGridType::BoundaryGridPartType, Dune::Stuff::Common::LogStream >
            ::Codim< MultiscaleGridProviderType::dim, 0 >
            ::entities(globalGridPart, boundaryGridPart, "  ", DSC_LOG_DEBUG);
      }
    }
}

TEST_F(CubeProvider, coupling_grid_parts)
{
  for (size_t ss = 0; ss < ms_grid_->size(); ++ss) {
    for (auto nn : ms_grid_->neighborsOf(ss)) {
      const auto coupling_grid_part = ms_grid_->couplingGridPart(ss, nn);
      DSC_LOG_INFO << "testing coupling grid part: " << ss << ", " << nn << std::endl;
      // walk the grid part
      for (typename CouplingGridPartType::template Codim< 0 >::IteratorType entityIterator = coupling_grid_part.template begin< 0 >();
           entityIterator != coupling_grid_part.template end< 0 >();
           ++entityIterator) {
        const typename CouplingGridPartType::template Codim< 0 >::EntityType& entity = *entityIterator;
        const unsigned int entityIndex = coupling_grid_part.indexSet().index(entity);
        DSC_LOG_INFO << "entity " << entityIndex << ", neighbors " << std::flush;
        // walk the intersections
        for (typename CouplingGridPartType::IntersectionIteratorType intersectionIterator = couplingGridPart.ibegin(entity);
             intersectionIterator != couplingGridPart.iend(entity);
             ++intersectionIterator) {
          const typename CouplingGridPartType::IntersectionIteratorType::Intersection& intersection = *intersectionIterator;
          const typename CouplingGridPartType::IntersectionIteratorType::Intersection::EntityPointer neighborPtr = intersection.outside();
          const typename CouplingGridPartType::template Codim< 0 >::EntityType& neighbor = *neighborPtr;
          const unsigned int neighborIndex = globalGridPart.indexSet().index(neighbor);
          DSC_LOG_INFO << neighborIndex << " ";
        } // walk the intersections
        DSC_LOG_INFO << std::endl;
      } // walk the grid part
    }
  }
}

TEST_F(CubeProvider, timings)
{
  // time grid parts
  DSC_LOG_INFO << "timing grid parts:" << std::endl;
  typedef MsGridType::GlobalGridPartType GlobalGridPartType;
  const Dune::shared_ptr< const GlobalGridPartType > globalGridPart = ms_grid_->globalGridPart();
  measureTiming(*globalGridPart, DSC_LOG_INFO, "global");
  typedef MsGridType::CouplingGridPartType CouplingGridPartType;
  const unsigned int neighbor = *(ms_grid_->neighborsOf(0)->begin());
  const Dune::shared_ptr< const CouplingGridPartType > couplingGridPart = ms_grid_->couplingGridPart(0, neighbor);
  measureTiming(*couplingGridPart, DSC_LOG_INFO, "coupling (subdomain 0 with " + Dune::Stuff::Common::String::convertTo(neighbor) + ")");
  typedef MsGridType::LocalGridPartType LocalGridPartType;
  const Dune::shared_ptr< const LocalGridPartType > firstLocalGridPart = ms_grid_->localGridPart(0);
  measureTiming(*firstLocalGridPart, DSC_LOG_INFO, "local (subdomain 0)");
  for (unsigned int subdomain = 1; subdomain < ms_grid_->size(); ++subdomain) {
    const Dune::shared_ptr< const LocalGridPartType > localGridPart = ms_grid_->localGridPart(subdomain);
    measureTiming(*localGridPart, DSC_LOG_DEBUG, "local (subdomain " + Dune::Stuff::Common::String::convertTo(subdomain) + ")");
  }
}

