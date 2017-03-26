# mcf

## Introduction

This library provides data structures and solvers for multi-frame data
association, suitable for the min-cost flow formulation of multiple object
tracking [1]. The library provides:

* Solvers for the multiple object tracking problem based on COIN-OR Clp and
  COIN-OR Lemon,
* an implementation of an efficient k-shortest path solver [2],
* a convenience class to optimize the sequence in batches over a fixed-length
  optimization window,
* Python bindings.

## Dependencies

Software development is carried out in Linux. The library should compile and
run on other systems with few changes to the build system, but this has
not been tested. The core library has no dependencies other than a recent
version of CMake (>= 3.2.0) and a C++14 compliant compiler. However, depending
on the specific feature set, optional dependencies may be required:

* [Boost](http://www.boost.org/) (for an alternative Dijkstra implementation)
* [pybind11](https://github.com/pybind/pybind11) (for Python bindings)
* [COIN-OR linear programming](https://projects.coin-or.org/Clp) (for Clp solver)
* [COIN-OR LEMON](http://lemon.cs.elte.hu/trac/lemon) (for Lemon solver)

The library ships with a Makefile to install these dependencies locally.

## Installation

First, clone the repository and create a build directory:
```
git clone https://github.com/nwojke/mcf.git
mkdir mcf/build
cd mcf/build
```
Then, install dependencies using the provided Makefile (these commands must be
called from within the build directory):
```
# The following command installs pybind11 to enable Python bindings. This
# should usually be sufficient.
make -f ../Makefile.external pybind11

# Alternatively, you may choose to install all dependencies, but this may take
# a while.
make -f ../Makefile.external
```
Then, configure the project and build the library (again, we are inside the
build directory):
```
# -DMCF_BUILD_PYTHON=ON to explicitly enable Python binds
# -DMCF_BUILD_STATIC=ON to compile a static library
cmake -DCMAKE_BUILD_TYPE=RELEASE -DMCF_BUILD_PYTHON=ON -DMCF_BUILD_STATIC=ON ..
make
```
On completion, the static library should reside in `build/lib`, the compiled
Python module in `build/python_lib`, and library header files in
`build/include`. Installation scripts are currently not provided; simply
copy the files over to your development environment or modify the include
and library paths accordingly.

## Getting started

The main data structure is the tracking graph defined in `graph.hpp`.
An efficient k-shortest path solver is defined in `k_shortest_path_solver.hpp`.
The batch solver is defined in `batch_processing.hpp`.
All of the classes have accompanying documentation that should help you get
started. For a minimalistic example check out the `examples` directory.

## Literature

1. Zhang, Li, Nevatia: Global data association for multi-object tracking
using network flows, CVPR (2008).
2. Berclaz, Fleuret, Tueretken, Fua: Multiple Object Tracking using
K-shortes Paths, PAMI, 33(9), 2011.

