// vim: expandtab:ts=2:sw=2
#ifndef MCF_CLP_SOLVER_HPP
#define MCF_CLP_SOLVER_HPP
#ifdef MCF_USE_Clp

#include <mcf/graph.hpp>
#include <mcf/solver.hpp>
#include <memory>
#include <vector>

class ClpSimplex;

namespace mcf {

/**
 * This class implements a min-cost flow solver that calls the COIN-OR linear
 * programming (CLP) library [1]. CLP uses a Simplex implementation to solve
 * the linear program.
 *
 * The RunSearch() method evaluates the linear program multiple times. It is
 * implemented using a binary search to find the min-cost solution in the given
 * range.
 *
 * The cost value returned by Run() and RunSearch() is (approximately) exact,
 * i.e., it is the sum of all edge costs in the multi-object trajectory.
 *
 * [1] https://projects.coin-or.org/Clp
 */
class ClpSolver : public Solver {
 public:
  //! Empty constructor.
  ClpSolver();

  //! Construct solver and initialize with given graph structure.
  ClpSolver(const Graph& graph);

  ~ClpSolver();

  void Build(const Graph& graph) override;

  double Run(int flow,
             std::vector<std::vector<int>>& trajectories) const override;

  double RunSearch(int min_flow, int max_flow,
                   std::vector<std::vector<int>>& trajectories) const override;

 private:
  std::unique_ptr<ClpSimplex> model_;
  std::vector<Edge> edges_;
  std::vector<std::vector<int>> node_index_to_incoming_edges_;
  std::vector<std::vector<int>> node_index_to_outgoing_edges_;
};

}  // namespace mcf

#endif  // MCF_USE_Clp
#endif  // MCF_CLP_SOLVER_HPP
