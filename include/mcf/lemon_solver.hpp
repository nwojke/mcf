// vim: expandtab:ts=2:sw=2
#ifndef MCF_LEMON_SOLVER_HPP
#define MCF_LEMON_SOLVER_HPP
#ifdef MCF_USE_Lemon

#include <lemon/cost_scaling.h>
#include <lemon/smart_graph.h>
#include <mcf/graph.hpp>
#include <mcf/solver.hpp>
#include <memory>

namespace mcf {

/**
 * This class implements a min-cost flow solver that calls a preflow
 * push-relabel algorithm from the COIN-OR LEMON library [1].
 *
 * The RunSearch() method evaluates the linear program multiple times. It is
 * implemented using a binary search to find the min-cost solution in the given
 * range.
 *
 * The cost value returned by Run() and RunSearch() is an internal cost
 * representation that can be used to compare solutions (smaller is better),
 * but differs from the exact cost of the multi-object trajectory (which is the
 * sum of all edge costs on the trajectory).
 *
 * [1] https://lemon.cs.elte.hu/trac/lemon
 */
class LemonSolver : public Solver {
 public:
  //! Empty constructor
  LemonSolver();

  //! Construct solver and initialize with given graph structure.
  LemonSolver(const Graph& graph);

  void Build(const Graph& graph) override;

  double Run(int flow,
             std::vector<std::vector<int>>& trajectories) const override;

  double RunSearch(int min_flow, int max_flow,
                   std::vector<std::vector<int>>& trajectories) const override;

 private:
  typedef lemon::SmartDigraph LemonGraph;
  typedef LemonGraph::Node Node;
  typedef LemonGraph::Arc Arc;

  typedef lemon::CostScaling<LemonGraph> Solver;
  typedef LemonGraph::ArcMap<Solver::Cost> CostMap;
  typedef LemonGraph::ArcMap<int> ArcMap;
  typedef LemonGraph::NodeMap<int> NodeMap;

  void Reset(std::size_t num_nodes, std::size_t num_arcs);

  std::vector<LemonGraph::Node> nodes_;
  std::vector<LemonGraph::Arc> arcs_;

  std::unique_ptr<LemonGraph> graph_;
  std::unique_ptr<CostMap> cost_map_;
  std::unique_ptr<ArcMap> lower_map_;
  std::unique_ptr<ArcMap> upper_map_;
  std::unique_ptr<Solver> solver_;
};

}  // namespace mcf

#endif  // MCF_USE_Lemon
#endif  // MCF_LEMON_SOLVER_HPP
