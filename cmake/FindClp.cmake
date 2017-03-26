# vim: expandtab:ts=2:sw=2
# Find coin-or linear programming solver
#
# Once done this will define
#   Clp_FOUND               - System has CLP
#   Clp_INCLUDE_DIRS        - The lemon include directories
#   Clp_LIBRARIES           - The libraries needed to use CLP
#   Clp_DEFINITIONS         - Compiler switches required for CLP
#
# You may define
#   Clp_PREFIX
# to define a custom installation prefix.
find_package(PkgConfig)
pkg_check_modules(PC_Clp QUIET clp)
set(Lemon_DEFINITIONS ${PC_Clp_CFLAGS_OTHER})

find_path(Clp_INCLUDE_DIR ClpSimplex.hpp HINTS
  ${Clp_PREFIX}/include/coin
  ${PC_Clp_INCLUDE_DIR}
  ${PC_Clp_INCLUDE_DIRS}
  /usr/include/coin
  )

set(Clp_LIBRARY_SEARCH_DIRS
  ${Clp_PREFIX}/lib
  ${PC_Clp_LIBDIR}
  ${PC_Clp_LIBRARY_DIRS}
  /usr/lib
  /usr/lib64
  /usr/lib32
  )

find_library(Clp_LIBRARY NAMES Clp HINTS
  ${Clp_LIBRARY_SEARCH_DIRS}
  )

find_library(ClpSolver_LIBRARY NAMES ClpSolver HINTS
  ${Clp_LIBRARY_SEARCH_DIRS}
  )

find_library(OsiClp_LIBRARY NAMES OsiClp HINTS
  ${Clp_LIBRARY_SEARCH_DIRS}
  )

find_library(CoinUtils_LIBRARY NAMES CoinUtils HINTS
  ${Clp_LIBRARY_SEARCH_DIRS}
  )

set(Clp_INCLUDE_DIRS
  ${Clp_INCLUDE_DIR}
  )

set(Clp_LIBRARIES
  ${Clp_LIBRARY}
  ${ClpSolver_LIBRARY}
  ${OsiClp_LIBRARY}
  ${CoinUtils_LIBRARY}
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Clp DEFAULT_MSG
  Clp_INCLUDE_DIR
  Clp_LIBRARY
  ClpSolver_LIBRARY
  CoinUtils_LIBRARY
  OsiClp_LIBRARY
  )

mark_as_advanced(
  Clp_INCLUDE_DIR
  Clp_LIBRARY_SEARCH_DIRS
  Clp_LIBRARY
  ClpSolver_LIBRARY
  OsiClp_LIBRARY
  CoinUtils_LIBRARY
  )

