  # -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid
# this implies checking for [dune-common], [dune-geometry], [dune-grid]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/grid/Makefile
  dune/grid/multiscale/Makefile
  dune/grid/multiscale/provider/Makefile
  examples/Makefile
  examples/provider/Makefile
  m4/Makefile
  dune-grid-multiscale.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL