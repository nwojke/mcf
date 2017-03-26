// vim: expandtab:ts=2:sw=2
#ifdef MCF_USE_Lemon
#include <mcf/lemon_solver.hpp>

#include <cassert>
#include <cmath>

#include <mcf/internal/util.hpp>
#include <mcf/logging.hpp>

namespace {

constexpr int MaxIntCost = 50000;

//! Convert a floating-point cost value to an integer cost for use with LEMON.
int to_integer_cost(const double x, const double max_abs_cost) {
  const int cost = static_cast<int>(std::round(MaxIntCost * x / max_abs_cost));
  return std::min(MaxIntCost, std::max(-MaxIntCost, cost));
}

}  // unnamed namespace

namespace mcf {

LemonSolver::LemonSolver() = default;

LemonSolver::LemonSolver(const Graph& graph) { Build(graph); }

void LemonSolver::Build(const Graph& graph) {
  double max_abs_cost = std::numeric_limits<double>::lowest();
  for (const Edge& edge : graph.edges()) {
    max_abs_cost = std::max(max_abs_cost, std::abs(edge.cost));
  }

  const std::vector<Edge>& edges = graph.edges();
  Reset(graph.num_nodes(), edges.size());

  for (int node_index = 0; node_index < graph.num_nodes(); ++node_index) {
    nodes_[node_index] = graph_->addNode();
  }

  for (std::size_t edge_index = 0; edge_index < edges.size(); ++edge_index) {
    arcs_[edge_index] = graph_->addArc(nodes_[edges[edge_index].source_index],
                                       nodes_[edges[edge_index].target_index]);
  }

  upper_map_.reset(new ArcMap(*graph_, 1));
  lower_map_.reset(new ArcMap(*graph_, 0));
  cost_map_.reset(new ArcMap(*graph_));
  for (std::size_t i = 0; i < edges.size(); ++i) {
    (*cost_map_)[arcs_[i]] = to_integer_cost(edges[i].cost, max_abs_cost);
  }

  solver_.reset(new Solver(*graph_));
  solver_->costMap(*cost_map_);
  solver_->lowerMap(*lower_map_);
  solver_->upperMap(*upper_map_);
}

double LemonSolver::Run(const int flow,
                        std::vector<std::vector<int>>& trajectories) const {
  ASSERT_NOTNULL(solver_);  // Uninitialized

  solver_->stSupply(nodes_[Graph::InternalSourceNode],
                    nodes_[Graph::InternalSinkNode], flow);

  Solver::ProblemType solver_state = solver_->run();
  if (solver_state != Solver::OPTIMAL) {
    // Infeasible or unbounded.
    trajectories.clear();
    return std::numeric_limits<double>::infinity();
  }

  std::vector<int> trajectory_heads;
  trajectory_heads.reserve(flow);

  std::vector<int> node_index_to_successor(graph_->nodeNum(), -1);
  for (int edge_index = 0; edge_index < graph_->arcNum(); ++edge_index) {
    const Arc& arc = arcs_[edge_index];
    if (solver_->flow(arc) < 1) {
      continue;
    }

    const int source_index = graph_->id(graph_->source(arc));
    const int target_index = graph_->id(graph_->target(arc));
    if (source_index == Graph::InternalSourceNode) {
      trajectory_heads.push_back(target_index);
      continue;
    }

    assert(node_index_to_successor[source_index] < 0 &&
           "duplicated association");
    node_index_to_successor[source_index] = target_index;
  }

  trajectories =
      internal::ComputeTrajectories(trajectory_heads, node_index_to_successor);
  return solver_->totalCost<double>();
}

double LemonSolver::RunSearch(
    const int min_flow, const int max_flow,
    std::vector<std::vector<int>>& trajectories) const {
  ASSERT_NOTNULL(solver_);  // Uninitialized

  auto EvaluateFlow = [this](const int flow) {
    std::vector<std::vector<int>> trajectories_at_flow;
    double cost_at_flow;
    try {
      cost_at_flow = this->Run(flow, trajectories_at_flow);
    } catch (std::exception&) {
      // Infeasible solution
      cost_at_flow = std::numeric_limits<double>::infinity();
    }

    print("Network evaluation, flow: ", flow, " cost: ", cost_at_flow);
    return std::make_tuple(cost_at_flow, trajectories_at_flow);
  };

  double cost;
  std::tie(cost, trajectories) =
      internal::BinarySearch(EvaluateFlow, min_flow, max_flow);
  return cost;
}

void LemonSolver::Reset(const std::size_t num_nodes,
                        const std::size_t num_arcs) {
  graph_.reset(new LemonGraph());
  graph_->reserveNode(num_nodes);
  graph_->reserveArc(num_arcs);

  nodes_.clear();
  arcs_.clear();
  nodes_.reserve(num_nodes);
  arcs_.reserve(num_arcs);

  solver_.reset();
  cost_map_.reset();
  lower_map_.reset();
  upper_map_.reset();
}

}  // namespace mcf

#endif  // MCF_USE_Lemon
