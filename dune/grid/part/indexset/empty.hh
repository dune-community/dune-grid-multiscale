#ifndef DUNE_GRID_PART_INDEXSET_EMPTY_HH
#define DUNE_GRID_PART_INDEXSET_EMPTY_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// system
#include <iostream>
#include <string>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cassert>

// dune-common
#include <dune/common/misc.hh>

namespace Dune {

namespace grid {

namespace Part {

namespace IndexSet {

/*!
  The EmptyIndexSet implements all additional method of a DUNE fem index set with
  an empty default implementation.
*/
class Empty
{
  // dummy value
  enum { myType = -1 };
public:
  //! return false mean the no memory has to be allocated
  //! and no compress of date has to be done
  bool compress () {
    return false;
  }

  //! return true if the index set is consecutive
  bool consecutive () const
  {
    return false;
  }

  //! return true if the index set is persistent
  bool persistent () const
  {
    return false;
  }

  //! do nothing here, because fathers index should already exist
  template< class EntityType >
  void insertEntity ( const EntityType &entity )
  {}

  //! do nothing here, because fathers index should already exist
  template< class EntityType >
  void removeEntity ( const EntityType &entity )
  {}

  //! nothing to do here
  void resize ()
  {}

  //! no extra memory for restriction is needed
  int additionalSizeEstimate () const
  {
    return 0;
  }

  static int type ()
  {
    return myType;
  }

  //! we have no old size
  int numberOfHoles ( const int /*codim*/ ) const
  {
    return 0;
  }

  //! return old index, for dof manager only
  int oldIndex ( const int /*hole*/, const int /*codim*/ ) const
  {
    return 0;
  }

  //! return new index, for dof manager only
  int newIndex ( const int /*hole*/, const int /*codim*/ ) const
  {
    return 0;
  }

  //! write index set to xdr file
  bool write_xdr ( const std::string &filename, int timestep )
  {
    FILE  *file;
    XDR   xdrs;
    const char *path = "";

    std::string fnstr  = genFilename(path,filename, timestep);
    const char * fn = fnstr.c_str();
    file = fopen(fn, "wb");
    if (!file)
    {
        std::cerr << "\aERROR in DefaultGridIndexSet::write_xdr(..): could not open <"
                  << filename << ">!" << std::endl;
      return false;
    }

    xdrstdio_create(&xdrs, file, XDR_ENCODE);
    this->processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);

    return true;
  }

  //! read index set to xdr file
  bool read_xdr ( const std::string &filename, int timestep )
  {
    FILE   *file;
    XDR     xdrs;
    const char *path = "";

    std::string fnstr = genFilename(path,filename, timestep);
    const char * fn  = fnstr.c_str();
    std::cout << "Reading <" << fn << "> \n";
    file = fopen(fn, "rb");
    if(!file)
    {
      std::cerr << "\aERROR in DefaultGridIndexSet::read_xdr(..): could not open <"
                << filename << ">!" << std::endl;
      return(false);
    }

    // read xdr
    xdrstdio_create(&xdrs, file, XDR_DECODE);
    this->processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);
    return true;
  }

protected:
  // read/write from/to xdr stream
  bool processXdr( XDR *xdrs )
  {
    int type = myType;
    xdr_int ( xdrs, &type);
    if(type != myType)
    {
      std::cerr << "\nERROR: DefaultGridIndex: wrong type choosen! \n\n";
      assert(type == myType);
    }
    return true;
  }
}; // class Empty

} // namespace IndexSet

} // namespace Part

} // namespace grid

} // namespace Dune

#endif // DUNE_GRID_PART_INDEXSET_EMPTY_HH
