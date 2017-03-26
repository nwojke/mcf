// vim: expandtab:ts=2:sw=2
#ifndef MCF_SOLVER_HPP
#define MCF_SOLVER_HPP

#include <mcf/graph.hpp>
#include <vector>

namespace mcf {

//! Abstract base class for all min-cost flow solvers.
class Solver {
 public:
  virtual ~Solver() {}

  /**
   * Build internal representation from given graph structure.
   *
   * Subsequent calls to Run() and RunSearch() operate on the provided
   * graph. Must be called once before any call to Run() or RunSearch().
   *
   * @param graph The graph structure.
   */
  virtual void Build(const Graph& graph) = 0;

  /**
   * Find min-cost solution for a given flow / number of trajectories.
   *
   * @param flow Flow value / number of trajectories.
   * @param trajectories Computed trajectories, where each trajectory is a
   *        sequence of location handles.
   * @return An internal cost representation that can be used to compare the
   *         quality of solutions obtained from the same graph (smaller is
   *         better). Returns std::numeric_limits<double>::infinity() if the
   *         solution is infeasible.
   *
   * @throws std::runtime_error If the solver is uninitialized (no graph was
   *         set using Build()).
   */
  virtual double Run(int flow,
                     std::vector<std::vector<int>>& trajectories) const = 0;

  /**
   * Find min-cost solution with flow / number of trajectories
   * in [min_flow, max_flow].
   *
   * @param min_flow Lower flow bound.
   * @param max_flow Upper flow bound.
   * @param trajectories Computed trajectories, where each trajectory is a
   *        sequence of location handles.
   * @return An internal cost representation that can be used to compare the
   *         quality of solutions obtained from the same graph (smaller is
   *         better). Returns std::numeric_limits<double>::infinity() if the
   *         solution is infeasible.
   *
   * @throws std::runtime_error If the solver is uninitialized (no graph was
   *         set using Build()).
   */
  virtual double RunSearch(
      int min_flow, int max_flow,
      std::vector<std::vector<int>>& trajectories) const = 0;
};

}  // namespace mcf

#endif
