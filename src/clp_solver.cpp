// vim: expandtab:ts=2:sw=2
#ifdef MCF_USE_Clp
#include <mcf/clp_solver.hpp>

#include <ClpSimplex.hpp>
#include <CoinBuild.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <sstream>
#include <vector>

#include <mcf/internal/util.hpp>
#include <mcf/logging.hpp>

namespace mcf {

namespace {

/**
 * Build linear program.
 *
 * @param edges Edges in the graph.
 * @param node_index_to_incoming_edges Contains the computed incoming edge
 *        indices for every node in the graph.
 * @param node_index_to_outgoing_edges Contains the computed outgoing edege
 *        indices for every node in the graph.
 * @param model The CLP simplex model to write to.
 */
void BuildLinearProgram(
    const std::vector<Edge>& edges,
    const std::vector<std::vector<int>>& node_index_to_incoming_edges,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    const int flow, ClpSimplex& model) {
  assert(node_index_to_incoming_edges.size() ==
             node_index_to_outgoing_edges.size() &&
         "node_index_to_incoming_edges.size() != "
         "node_index_to_outgoing_edges.size()");
  const int num_nodes = node_index_to_incoming_edges.size();

  // Set the objective coefficients and bounds.
  model.resize(0, edges.size());
  for (std::size_t edge_index = 0; edge_index < edges.size(); ++edge_index) {
    const Edge& edge = edges[edge_index];
    model.setObjectiveCoefficient(edge_index, edge.cost);
    model.setColumnBounds(edge_index, 0.0, 1.0);
  }

  // Build constraints.
  static constexpr int kBuildRows = 0;
  CoinBuild build_object(kBuildRows);
  std::vector<int> column_index_buffer(edges.size());
  std::vector<double> column_value_buffer(edges.size());

  // Constraint on the outgoing flow of the source node.
  const std::size_t num_outgoing_source_edges =
      node_index_to_outgoing_edges[Graph::InternalSourceNode].size();
  std::fill_n(column_value_buffer.begin(), num_outgoing_source_edges, 1.0);
  build_object.addRow(
      num_outgoing_source_edges,
      node_index_to_outgoing_edges[Graph::InternalSourceNode].data(),
      column_value_buffer.data(), 0.0, 0.0);

  // Constraint on the incoming flow of the sink node.
  const std::size_t num_incoming_sink_edges =
      node_index_to_incoming_edges[Graph::InternalSinkNode].size();
  std::fill_n(column_value_buffer.begin(), num_incoming_sink_edges, 1.0);
  build_object.addRow(
      num_incoming_sink_edges,
      node_index_to_incoming_edges[Graph::InternalSinkNode].data(),
      column_value_buffer.data(), 0.0, 0.0);

  // Flow conservation on remaining nodes.
  for (int node_index = Graph::FirstNonSourceSinkNode; node_index < num_nodes;
       ++node_index) {
    const std::size_t num_incoming_edges =
        node_index_to_incoming_edges[node_index].size();
    const std::size_t num_outgoing_edges =
        node_index_to_outgoing_edges[node_index].size();

    std::copy(node_index_to_incoming_edges[node_index].begin(),
              node_index_to_incoming_edges[node_index].end(),
              column_index_buffer.begin());
    std::fill_n(column_value_buffer.begin(), num_incoming_edges, 1.0);

    std::copy(node_index_to_outgoing_edges[node_index].begin(),
              node_index_to_outgoing_edges[node_index].end(),
              column_index_buffer.begin() + num_incoming_edges);
    std::fill_n(column_value_buffer.begin() + num_incoming_edges,
                num_outgoing_edges, -1.0);

    static constexpr int kLowerBound = 0.0;
    static constexpr int kUpperBound = 0.0;
    build_object.addRow(num_incoming_edges + num_outgoing_edges,
                        column_index_buffer.data(), column_value_buffer.data(),
                        kLowerBound, kUpperBound);
  }

  // Add rows to model.
  model.addRows(build_object);
}

/**
 * Solve CLP simplex model for a given flow / number of trajectories.
 *
 * @param flow The flow / number of trajectories.
 * @param model A CLP simplex model that has previously been constructed with
 *        BuildLinearProgram().
 * @return True if the optimal solution has been found, false otherwise.
 */
bool SolveLinearProgram(const int flow, ClpSimplex& model) {
  // Set flow and solve problem.
  model.setRowBounds(Graph::InternalSourceNode, flow, flow);
  model.setRowBounds(Graph::InternalSinkNode, flow, flow);
  model.primal();

  // Check that the result is integer.
  const double* solution = model.primalColumnSolution();
  for (int i = 0; i < model.numberColumns(); ++i) {
    static constexpr double kEpsilon = 1e-5;
    if (std::abs(solution[i]) < kEpsilon ||
        std::abs(solution[i] - 1.0) < kEpsilon) {
      continue;
    }
    return false;
  }

  return model.isProvenOptimal();
}

/**
 * Build a vector that maps each node index to its successor on the
 * multi-object trajectory.
 *
 * @param edges Edges in the graph.
 * @param node_index_to_incoming_edges Contains the computed incoming edge
 *        indices for every node in the graph.
 * @param node_index_to_outgoing_edges Contains the computed outgoing edege
 *        indices for every node in the graph.
 * @param model A CLP simplex model that has been solved by
 *        SolveLinearProgram().
 * @param flow The flow value / number of trajectories with which the model has
 *        been solved.
 * @param trajectory_heads Computed trajectory entry points (node indices).
 * @param node_index_to_successor Computed vector that maps each node index to
 *        the index of the subsequent node on the trajectory. Each trajectory
 *        ends with Graph::InternalSinkNode and nodes that are not part of any
 *        trajectory are set to -1.
 */
void ComputeSuccessorMap(
    const std::vector<Edge>& edges,
    const std::vector<std::vector<int>>& node_index_to_incoming_edges,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    const ClpSimplex& model, const int flow, std::vector<int>& trajectory_heads,
    std::vector<int>& node_index_to_successor) {
  node_index_to_successor.resize(model.numberColumns());
  std::fill(node_index_to_successor.begin(), node_index_to_successor.end(), -1);

  trajectory_heads.clear();
  trajectory_heads.reserve(flow);

  const double* solution = model.primalColumnSolution();
  for (int edge_index = 0; edge_index < model.numberColumns(); ++edge_index) {
    if (solution[edge_index] < 0.5) {
      continue;
    }

    const int& source_index = edges[edge_index].source_index;
    const int& target_index = edges[edge_index].target_index;

    if (source_index == Graph::InternalSourceNode) {
      trajectory_heads.push_back(edges[edge_index].target_index);
      continue;
    }

    assert(node_index_to_successor[source_index] < 0 &&
           "duplicate association");
    node_index_to_successor[source_index] = target_index;
  }
}

/**
 * Compute objective value of a solved model.
 *
 * @param model A CLP simplex model that has previously been solved by
 *        SolveLinearProgram().
 * @parm edges The graphs edge list corresponding to the simplex model.
 * @return Sum of all edge costs that are part of the multi-object trajectory.
 */
double ComputeSolutionCost(const ClpSimplex& model,
                           const std::vector<Edge>& edges) {
  const double* solution = model.primalColumnSolution();
  double cost = 0.0;
  for (std::size_t edge_index = 0; edge_index < edges.size(); ++edge_index) {
    cost += solution[edge_index] * edges[edge_index].cost;
  }
  return cost;
}

}  // unnamed namespace

ClpSolver::ClpSolver() = default;

ClpSolver::~ClpSolver() = default;

ClpSolver::ClpSolver(const Graph& graph) { Build(graph); }

void ClpSolver::Build(const Graph& graph) {
  model_.reset(new ClpSimplex());
  edges_ = graph.edges();
  internal::BuildEdgeMap(graph, node_index_to_incoming_edges_,
                         node_index_to_outgoing_edges_);

  static constexpr int kDefaultFlow = 0;
  BuildLinearProgram(edges_, node_index_to_incoming_edges_,
                     node_index_to_outgoing_edges_, kDefaultFlow, *model_);
}

double ClpSolver::Run(const int flow,
                      std::vector<std::vector<int>>& trajectories) const {
  ASSERT_NOTNULL(model_);  // Uninitialized

  model_->setLogLevel(mcf::printer().is_verbose() ? 1 : 0);
  if (!SolveLinearProgram(flow, *model_)) {
    // Infeasible or unbounded solution.
    trajectories.clear();
    return std::numeric_limits<double>::infinity();
  }

  std::vector<int> trajectory_heads;
  std::vector<int> node_index_to_successor;
  ComputeSuccessorMap(edges_, node_index_to_incoming_edges_,
                      node_index_to_outgoing_edges_, *model_, flow,
                      trajectory_heads, node_index_to_successor);
  trajectories =
      internal::ComputeTrajectories(trajectory_heads, node_index_to_successor);
  return ComputeSolutionCost(*model_, edges_);
}

double ClpSolver::RunSearch(const int min_flow, const int max_flow,
                            std::vector<std::vector<int>>& trajectories) const {
  ASSERT_NOTNULL(model_);  // Uninitialized

  auto EvaluateFlow = [this](const int flow) {
    std::vector<std::vector<int>> trajectories_at_flow;
    double cost_at_flow;
    try {
      cost_at_flow = this->Run(flow, trajectories_at_flow);
    } catch (std::exception&) {
      // Infeasible solution.
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

}  // namespace mcf

#endif  // MCF_USE_Clp
