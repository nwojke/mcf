// vim: expandtab:ts=2:sw=2
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <mcf/batch_processing.hpp>
#include <mcf/graph.hpp>
#include <mcf/k_shortest_path_solver.hpp>
#include <mcf/logging.hpp>

#ifdef MCF_USE_Lemon
#include <mcf/lemon_solver.hpp>
#endif

#ifdef MCF_USE_Clp
#include <mcf/clp_solver.hpp>
#endif

namespace py = pybind11;
using namespace py::literals;

//! A list of available solver names that can be passed on to the solver class.
static std::vector<std::string> kAvailableSolvers = {
    "ssp",
#ifdef MCF_USE_Lemon
    "lemon",
#endif
#ifdef MCF_USE_Clp
    "clp",
#endif
};

// A small wrapper class around mcf::Graph with py::dict node attribute storage.
class PyGraph : public mcf::Graph {
 public:
  PyGraph() = default;

  void Reserve(const int num_edges) { mcf::Graph::Reserve(num_edges); }

  //! Add a location to the graph and set location specific attributes.
  int Add(const double cost, const py::dict& attributes) {
    const int location_handle = mcf::Graph::Add(cost);
    location_attributes_[location_handle] = attributes;
    return location_handle;
  }

  void Link(const int src, const int dst, const double cost) {
    mcf::Graph::Link(src, dst, cost);
  }

  //! Get location specific attributes.
  py::dict operator[](const int location_handle) {
    return location_attributes_[location_handle];
  }

  //! Return internal data as a pickable object.
  py::tuple GetPickableState() const {
    std::vector<py::tuple> edges(edges_.size());
    for (std::size_t i = 0; i < edges_.size(); ++i) {
      const mcf::Edge& edge = edges_[i];
      edges[i] =
          py::make_tuple(edge.source_index, edge.target_index, edge.cost);
    }

    return py::make_tuple(edges, next_id_, location_attributes_);
  }

  //! Set internal data from pickle object.
  void SetPickledState(const py::tuple& data) {
    if (data.size() != 3) {
      throw std::runtime_error("Invalid state.");
    }

    py::list edge_list = data[0].cast<py::list>();
    edges_.resize(edge_list.size());
    for (std::size_t i = 0; i < edge_list.size(); ++i) {
      std::tuple<int, int, double> edge_element =
          edge_list[i].cast<std::tuple<int, int, double>>();
      edges_[i].source_index = std::get<0>(edge_element);
      edges_[i].target_index = std::get<1>(edge_element);
      edges_[i].cost = std::get<2>(edge_element);
    }

    next_id_ = data[1].cast<int>();
    location_attributes_ = data[2].cast<std::unordered_map<int, py::dict>>();
  }

 private:
  std::unordered_map<int, py::dict> location_attributes_;
};

//! A small Python wrapper around the solver class.
class PySolver {
 public:
  //! Create solver with appropriate method (see kAvailableSolvers).
  PySolver(const std::string& method) {
    if (method == "ssp") {
      solver_ = std::make_unique<mcf::ShortestPathSolver>();
    }
#ifdef MCF_USE_Lemon
    else if (method == "lemon") {
      solver_ = std::make_unique<mcf::LemonSolver>();
    }
#endif
#ifdef MCF_USE_Clp
    else if (method == "clp") {
      solver_ = std::make_unique<mcf::ClpSolver>();
    }
#endif
    else {
      std::stringstream msg;
      msg << "Unknown method '" << method << "', must be one of the following:";
      for (const auto& name : kAvailableSolvers) {
        msg << " '" << name << "'";
      }
      throw std::invalid_argument(msg.str().c_str());
    }
  }

  //! Create solver with appropriate method and initialize with given graph.
  PySolver(const mcf::Graph& graph, const std::string& method)
      : PySolver(method) {
    Build(graph);
  }

  //! Wrapper function that releases the GIL and calls the solver.
  void Build(const mcf::Graph& graph) {
    py::gil_scoped_release release;
    (void)release;
    solver_->Build(graph);
  }

  //! Wrapper function that releases the GIL and calls the solver.
  std::vector<std::vector<int>> Run(const int flow) {
    py::gil_scoped_release release;
    (void)release;

    std::vector<std::vector<int>> trajectories;
    solver_->Run(flow, trajectories);
    return trajectories;
  }

  //! Wrapper function that releases the GIL and calls the solver.
  std::vector<std::vector<int>> RunSearch(const int min_flow,
                                          const int max_flow) {
    py::gil_scoped_release release;
    (void)release;

    std::vector<std::vector<int>> trajectories;
    solver_->RunSearch(min_flow, max_flow, trajectories);
    return trajectories;
  }

 private:
  std::unique_ptr<mcf::Solver> solver_;
};

//! A small wrapper class around mcf::BatchProcessing with py::Dict node
// attribute storage.
class PyBatchProcessing {
 public:
  PyBatchProcessing(int window_len)
      : processor_(window_len), min_active_location_(0) {}

  void Reserve(const int num_edges) { processor_.Reserve(num_edges); }

  //! Add a location to the graph and set location specific attributes.
  int Add(const double cost, const py::dict& attributes) {
    const int location_handle = processor_.Add(cost);
    location_attributes_[location_handle] = attributes;
    return location_handle;
  }

  void Link(const int src, const int dst, const double cost) {
    processor_.Link(src, dst, cost);
  }

  void FinalizeTimeStep() { processor_.FinalizeTimeStep(); }

  std::vector<mcf::BatchProcessing::Trajectory> RunSearch(
      bool ignore_last_exit_cost) {
    std::vector<mcf::BatchProcessing::Trajectory> trajectories;
    processor_.RunSearch(trajectories, ignore_last_exit_cost);
    for (const auto& trajectory : trajectories) {
      if (trajectory.empty()) {
        continue;
      }
      min_active_location_ = std::min(min_active_location_, trajectory.front());
    }
    return trajectories;
  }

  mcf::BatchProcessing::TrajectoryMap ComputeTrajectories(
      bool ignore_last_exit_cost) {
    mcf::BatchProcessing::TrajectoryMap trajectories =
        processor_.ComputeTrajectories(ignore_last_exit_cost);
    for (const auto& index_and_trajectory : trajectories) {
      assert(!index_and_trajectory.second.empty());
      min_active_location_ =
          std::min(min_active_location_, index_and_trajectory.second.front());
    }
    return trajectories;
  }

  void RemoveInactiveTracks() {
    // Remove unusued tracks.
    processor_.RemoveInactiveTracks();

    // Remove unused attributes.
    AttributeMap new_attribute_map;
    for (const auto& location_and_attributes : location_attributes_) {
      if (location_and_attributes.first >= min_active_location_) {
        continue;
      }
      new_attribute_map.insert(location_and_attributes);
    }

    location_attributes_ = std::move(new_attribute_map);
  }

  //! Get location specific attributes.
  py::dict operator[](const int location_handle) {
    return location_attributes_[location_handle];
  }

 private:
  using AttributeMap =
      std::unordered_map<mcf::BatchProcessing::Index, py::dict>;

  mcf::BatchProcessing processor_;
  AttributeMap location_attributes_;

  // Location with smallest index on an active trajectory.
  mcf::BatchProcessing::Index min_active_location_;
};

PYBIND11_PLUGIN(mcf) {
  py::module module("mcf", "A small library to solve min-cost flow networks");

  module.def("set_verbose",
             [](const bool verbose) { mcf::printer().set_verbose(verbose); },
             "Set verbosity. If True, prints debug information to standard "
             "output.",
             "verbose"_a);

  module.def(
      "is_verbose", []() { return mcf::printer().is_verbose(); },
      "Check verbosity. If True, prints debug information to standard output.");

  py::class_<PyGraph>(module, "Graph")
      .def(py::init<>())
      .def("reserve", &PyGraph::Reserve, py::arg("num_edges"))
      .def("add", &PyGraph::Add, "cost"_a, "attributes"_a = py::dict())
      .def("link", &PyGraph::Link, "src"_a, "dst"_a, "cost"_a)
      .def("__getitem__", &PyGraph::operator[], "location"_a)
      .def("__getstate__",
           [](const PyGraph& graph) { return graph.GetPickableState(); })
      .def("__setstate__",
           [](PyGraph& graph, py::tuple data) {
             new (&graph) PyGraph();
             graph.SetPickledState(data);
           })
      .def_property_readonly("ST",
                             [](const py::object&) { return mcf::Graph::ST; });

  py::class_<PySolver>(module, "Solver")
      .def(py::init<const std::string&>(), "method"_a = std::string("ssp"))
      .def(py::init<const PyGraph&, const std::string&>(), "graph"_a,
           "method"_a = std::string("ssp"))
      .def("build", &PySolver::Build, "graph"_a)
      .def("run", &PySolver::Run, "flow"_a)
      .def("run_search", &PySolver::RunSearch, "min_flow"_a, "max_flow"_a);

  py::class_<PyBatchProcessing>(module, "BatchProcessing")
      .def(py::init<int>(), "window_len"_a)
      .def("reserve", &PyBatchProcessing::Reserve, py::arg("num_edges"))
      .def("add", &PyBatchProcessing::Add, "cost"_a,
           "attributes"_a = py::dict())
      .def("link", &PyBatchProcessing::Link, "src"_a, "dst"_a, "cost"_a)
      .def("__getitem__", &PyBatchProcessing::operator[], "location"_a)
      .def("finalize_timestep", &PyBatchProcessing::FinalizeTimeStep)
      .def("run_search", &PyBatchProcessing::RunSearch,
           "ignore_last_exit_cost"_a = true)
      .def("compute_trajectories", &PyBatchProcessing::ComputeTrajectories,
           "ignore_last_exit_cost"_a)
      .def("remove_inactive_tracks", &PyBatchProcessing::RemoveInactiveTracks)
      .def_property_readonly(
          "ST", [](const py::object&) { return mcf::BatchProcessing::ST; });

  return module.ptr();
}
