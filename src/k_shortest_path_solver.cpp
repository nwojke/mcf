// vim: expandtab:ts=2:sw=2
#include <mcf/k_shortest_path_solver.hpp>

#include <mcf/internal/k_shortest_path.hpp>
#include <mcf/internal/util.hpp>
#include <mcf/logging.hpp>

namespace mcf {

ShortestPathSolver::ShortestPathSolver(ShortestPathSolverType solver_type)
    : solver_type_(solver_type) {}

ShortestPathSolver::ShortestPathSolver(const Graph& graph,
                                       ShortestPathSolverType solver_type)
    : solver_type_(solver_type) {
  Build(graph);
}

void ShortestPathSolver::Build(const Graph& graph) {
  graph_ = graph;
  internal::BuildEdgeMap(graph_, node_index_to_incoming_edges_,
                         node_index_to_outgoing_edges_);
}

double ShortestPathSolver::Run(
    const int flow, std::vector<std::vector<int>>& trajectories) const {
  // Build residual graph.
  Graph residual_graph = graph_;
  std::vector<std::vector<int>> node_index_to_incoming_edges =
      node_index_to_incoming_edges_;
  std::vector<std::vector<int>> node_index_to_outgoing_edges =
      node_index_to_outgoing_edges_;

  // Compute residual graph.
  static constexpr bool kStopOnBest = false;

  int final_flow = 0;
  switch (solver_type_) {
    case ShortestPathSolverType::kDijkstraLazyDeletion:
      final_flow =
          internal::RunSuccessiveShortestPathSearchDijkstraLazyDeletion(
              residual_graph, node_index_to_incoming_edges,
              node_index_to_outgoing_edges, flow, kStopOnBest);
      break;
#ifdef MCF_USE_Boost
    case ShortestPathSolverType::kDijkstraFibonacciHeap:
      final_flow =
          internal::RunSuccessiveShortestPathSearchDijkstraFibonacciHeap(
              residual_graph, node_index_to_incoming_edges,
              node_index_to_outgoing_edges, flow, kStopOnBest);
      break;
#endif
  }

  if (final_flow != flow) {
    // Infeasible, no solution for given flow.
    trajectories.clear();
    return std::numeric_limits<double>::infinity();
  }

  // Extract trajectories.
  return internal::ComputeTrajectoriesFromResidualGraph(
      residual_graph, node_index_to_outgoing_edges, trajectories);
}

double ShortestPathSolver::RunSearch(
    int min_flow, int max_flow,
    std::vector<std::vector<int>>& trajectories) const {
  // Build residual graph.
  Graph residual_graph = graph_;
  std::vector<std::vector<int>> node_index_to_incoming_edges =
      node_index_to_incoming_edges_;
  std::vector<std::vector<int>> node_index_to_outgoing_edges =
      node_index_to_outgoing_edges_;

  // Compute residual graph.
  static constexpr bool kStopOnBest = true;
  int final_flow = 0;
  switch (solver_type_) {
    case ShortestPathSolverType::kDijkstraLazyDeletion:
      final_flow =
          internal::RunSuccessiveShortestPathSearchDijkstraLazyDeletion(
              residual_graph, node_index_to_incoming_edges,
              node_index_to_outgoing_edges, max_flow, kStopOnBest);
      break;
#ifdef MCF_USE_Boost
    case ShortestPathSolverType::kDijkstraFibonacciHeap:
      final_flow =
          internal::RunSuccessiveShortestPathSearchDijkstraFibonacciHeap(
              residual_graph, node_index_to_incoming_edges,
              node_index_to_outgoing_edges, max_flow, kStopOnBest);
      break;
#endif
  }
  if (final_flow < min_flow) {
    // Infeasible, no solution in given bounds.
    trajectories.clear();
    return std::numeric_limits<double>::infinity();
  }

  // Extract trajectories.
  return internal::ComputeTrajectoriesFromResidualGraph(
      residual_graph, node_index_to_outgoing_edges, trajectories);
}

}  // namespace mcf
