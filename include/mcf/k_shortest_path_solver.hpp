// vim: expandtab:ts=2:sw=2
#ifndef MCF_K_SHORTEST_PATH_SOLVER_HPP
#define MCF_K_SHORTEST_PATH_SOLVER_HPP

#include <mcf/graph.hpp>
#include <mcf/solver.hpp>
#include <vector>

namespace mcf {

/**
 * An enumerator type for the shortest path solver used in the successive
 * shortest paths algorithm.
 *
 * * kDijkstraLazyDeletion is an implementation of Dijkstra's algorithm that
 *   uses an immutable priority queue. This is usually faster, but may have
 *   higher memory consumption.
 * * kDijkstraFibonacciHeap is an implementation of Dijkstra's algorithm that
 *   uses a Fibonacci heap structure with decrease-key functionality. This is
 *   usually slower, but may have lower memory consumption.
 */
enum class ShortestPathSolverType {
  kDijkstraLazyDeletion,
#ifdef MCF_USE_Boost
  kDijkstraFibonacciHeap,
#endif
};

/**
 * This class provides an implementation of the K shortest path algorithm [1],
 * a specialized algorithm for solving min-cost network flow problems with unit
 * capacities.
 *
 * The solver expects the graph to be ordered, such that no location has a
 * larger handle/index than any of its successors. This is typically the case
 * in multiple object tracking where the graph is ordered by time.
 *
 * [1] Berclaz, Fleuret, Tueretken, Fua: Multiple Object Tracking using
 *     K-shortes Paths, PAMI, 33(9), 2011.
 */
class ShortestPathSolver : public Solver {
 public:
  //! Empty constructor.
  ShortestPathSolver(ShortestPathSolverType solver_type =
                         ShortestPathSolverType::kDijkstraLazyDeletion);

  /**
   * Construct solver and initialize with given graph structure.
   *
   * @param graph The graph structure.
   * @param solver_type The type of shortest path algorithm to use during
   *        successive shortest path search.
   */
  ShortestPathSolver(const Graph& graph,
                     ShortestPathSolverType solver_type =
                         ShortestPathSolverType::kDijkstraLazyDeletion);

  void Build(const Graph& graph) override;

  double Run(int flow,
             std::vector<std::vector<int>>& trajectories) const override;

  double RunSearch(int min_flow, int max_flow,
                   std::vector<std::vector<int>>& trajectories) const override;

 private:
  ShortestPathSolverType solver_type_;
  Graph graph_;
  std::vector<std::vector<int>> node_index_to_incoming_edges_;
  std::vector<std::vector<int>> node_index_to_outgoing_edges_;
};

}  // namespace mcf

#endif
