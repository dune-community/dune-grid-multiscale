#ifndef DUNE_STUFF_TEST_TOOLS_HH
#define DUNE_STUFF_TEST_TOOLS_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/mpihelper.hh>
#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/logging.hh>
#ifdef HAVE_DUNE_FEM
  #include <dune/fem/misc/mpimanager.hh>
#endif
#include <gtest.h>

#include <random>
#include <fstream>
#include <sys/time.h>

#include <dune/stuff/fem/namespace.hh>

template < template <class> class Test >
struct TestRunner {
    struct Visitor {
        template <class T>
        void visit(const T&) {
            Test<T>().run();
        }
    };

    template < class Tuple >
    static void run() {
        Tuple t;
        Dune::ForEachValue<Tuple> fe(t);
        Visitor v;
        fe.apply(v);
    }
};

template < int i >
struct Int {
  static const int value = i;
};

//! where sleep only counts toward wall time, this wastes actual cpu time
void busywait(const int ms)  {
  // "round" up to next full 10 ms to align with native timer res
  const int milliseconds = (ms/10)*10 + 10;
  timeval start, end;
  gettimeofday(&start, NULL);
  do  {
   gettimeofday(&end, NULL);
  } while( ((end.tv_sec - start.tv_sec )*1e6) + ((end.tv_usec - start.tv_usec)) < milliseconds * 1000 );
}


typedef Dune::tuple<double, float, //Dune::bigunsignedint,
  int, unsigned int, unsigned long, long long, char> BasicTypes;


void test_init(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  DSC_CONFIG.readOptions(argc, argv);
#ifdef HAVE_DUNE_FEM
#if DUNE_FEM_IS_MULTISCALE_COMPATIBLE
  Dune::MPIManager::initialize(argc, argv);
#elif DUNE_FEM_IS_LOCALFUNCTIONS_COMPATIBLE
  Dune::Fem::MPIManager::initialize(argc, argv);
#else
  Dune::Fem::MPIManager::initialize(argc, argv);
#endif
#else
  Dune::MPIHelper::instance(argc, argv);
#endif
  DSC::Logger().create(DSC::LOG_CONSOLE | DSC::LOG_ERROR);
}

#endif // DUNE_STUFF_TEST_TOOLS_HH
