// vim: expandtab:ts=2:sw=2
#include <mcf/internal/k_shortest_path.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <queue>

#ifdef MCF_USE_Boost
#include <boost/heap/fibonacci_heap.hpp>
#endif

#include <mcf/internal/util.hpp>
#include <mcf/logging.hpp>

namespace mcf {

namespace internal {

namespace {
//! A small value used for floating point comparison against 0.
constexpr double kEpsilon = 1e-8;

template <typename ShortestPathSolverFunction>
int RunSuccessiveShortestPathSearch(
    ShortestPathSolverFunction shortest_path_finder, Graph& residual_graph,
    std::vector<std::vector<int>>& node_index_to_incoming_edges,
    std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    const int max_flow, const bool stop_on_best) {
  print("Running K shortest path search, stop on best: ", stop_on_best,
        " max flow: ", max_flow);
  // If max_flow is 0, we have nothing to do.
  if (max_flow == 0) {
    print("K shortest path search has been invoked with maximum flow 0, ",
          "nothing to do.");
    return 0;
  }

  // Solve first iteration using dynamic programming and eliminate negative
  // cost, such that edges on the shortest path have cost 0.
  std::vector<double> node_index_to_distance;
  std::vector<int> node_index_to_predecessor_edge;
  std::vector<int> traversed_edges;
  double max_cost_sum = 0.0;

  ASSERT_TRUE(CheckGraphIsOrdered(residual_graph),
              "Expected an ordered graph, but graph is unordered.");
  ShortestPathInOrderedDAG(residual_graph, node_index_to_incoming_edges,
                           node_index_to_distance,
                           node_index_to_predecessor_edge);

  if (!FindPathFromSinkToSource(residual_graph, node_index_to_predecessor_edge,
                                traversed_edges, max_cost_sum)) {
    print("K shortest path search found infeasible solution for flow 1, ",
          "exiting.");
    return 0;
  }
  if (stop_on_best && max_cost_sum > 0) {
    print("Flow 1 has positive cost. Returning trivial solution (flow=0).");
    return 0;
  }
  max_cost_sum = -max_cost_sum;  // After conversion, all cost are positive.
  print("Flow 1, maximum cost sum: ", max_cost_sum);

  EliminateNegativeCost(node_index_to_distance, residual_graph);
  ReverseEdges(traversed_edges, residual_graph, node_index_to_incoming_edges,
               node_index_to_outgoing_edges);

  // Remaining iterations using Dijkstra.
  double cost_sum = 0.0;
  int flow = 1;
  while (flow < max_flow) {
    shortest_path_finder(residual_graph, node_index_to_outgoing_edges,
                         node_index_to_distance,
                         node_index_to_predecessor_edge);
    double path_cost = 0.0;
    if (!FindPathFromSinkToSource(residual_graph,
                                  node_index_to_predecessor_edge,
                                  traversed_edges, path_cost)) {
      print("K shortest path search stopped with infeasible solution at flow.");
      break;  // No feasible solution.
    }

    cost_sum += path_cost;
    if (stop_on_best && cost_sum > max_cost_sum) {
      print("K shortest path search stopped at optimal solution.");
      break;  // Previous iteration was the global minimum solution.
    }

    print("Flow ", flow + 1, " cost: ", path_cost, " sum: ", cost_sum);
    EliminateNegativeCost(node_index_to_distance, residual_graph);
    ReverseEdges(traversed_edges, residual_graph, node_index_to_incoming_edges,
                 node_index_to_outgoing_edges);
    ++flow;
  }

  print("K shortest path search final flow: ", flow);
  return flow;
}

}  // unnamed namespace

bool CheckGraphIsOrdered(const Graph& graph) {
  const std::vector<Edge>& edges = graph.edges();
  for (std::size_t edge_index = Graph::FirstNonSourceSinkNode;
       edge_index < edges.size(); ++edge_index) {
    const Edge& edge = edges[edge_index];
    if (edge.target_index != Graph::InternalSinkNode &&
        edge.target_index <= edge.source_index) {
      return false;
    }
  }
  return true;
}

void ShortestPathInOrderedDAG(
    const Graph& graph,
    const std::vector<std::vector<int>>& node_index_to_incoming_edges,
    std::vector<double>& node_index_to_distance,
    std::vector<int>& node_index_to_predecessor_edge) {
  // Initialize to "no predecessor" and infinity cost.
  node_index_to_predecessor_edge.resize(graph.num_nodes());
  std::fill_n(node_index_to_predecessor_edge.begin(), graph.num_nodes(), -1);
  node_index_to_distance.resize(graph.num_nodes());
  std::fill_n(node_index_to_distance.begin(), graph.num_nodes(),
              std::numeric_limits<double>::infinity());

  // Set source node to zero cost.
  node_index_to_distance[Graph::InternalSourceNode] = 0.0;

  // Sweep through graph and compute minimum distance from source to each node.
  const std::vector<Edge>& edges = graph.edges();
  auto ProcessNode = [&edges, &node_index_to_incoming_edges,
                      &node_index_to_distance,
                      &node_index_to_predecessor_edge](const int node_index) {
    double min_incoming_cost = std::numeric_limits<double>::infinity();
    int min_incoming_index = -1;
    for (const int edge_index : node_index_to_incoming_edges[node_index]) {
      const Edge& edge = edges[edge_index];
      const double path_cost =
          node_index_to_distance[edge.source_index] + edge.cost;

      if (path_cost < min_incoming_cost) {
        min_incoming_cost = path_cost;
        min_incoming_index = edge_index;
      }
    }
    node_index_to_distance[node_index] = min_incoming_cost;
    node_index_to_predecessor_edge[node_index] = min_incoming_index;
  };

  for (int node_index = Graph::FirstNonSourceSinkNode;
       node_index < graph.num_nodes(); ++node_index) {
    ProcessNode(node_index);
  }
  ProcessNode(Graph::InternalSinkNode);
}

void ShortestPathDijkstraLazyDeletion(
    const Graph& graph,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    std::vector<double>& node_index_to_distance,
    std::vector<int>& node_index_to_predecessor_edge) {
  // Set up heap structure and reserve space for all locations (each entry node
  // has an incoming edge originating from the source).
  std::vector<std::pair<double, int>> container;
  container.reserve(graph.num_nodes() / 2);

  std::priority_queue<std::pair<double, int>,
                      std::vector<std::pair<double, int>>,
                      std::greater<std::pair<double, int>>>
  queue(std::greater<std::pair<double, int>>(), std::move(container));

  // Initialize maps to "no predecessor" and infinity cost.
  node_index_to_predecessor_edge.resize(graph.num_nodes());
  std::fill_n(node_index_to_predecessor_edge.begin(), graph.num_nodes(), -1);
  node_index_to_distance.resize(graph.num_nodes());
  std::fill_n(node_index_to_distance.begin(), graph.num_nodes(),
              std::numeric_limits<double>::infinity());

  // Add source to queue.
  queue.emplace(0.0, Graph::InternalSourceNode);
  node_index_to_distance[Graph::InternalSourceNode] = 0.0;
  node_index_to_predecessor_edge[Graph::InternalSourceNode] = -1;

  // Iterate until all elements have been processed.
  while (!queue.empty()) {
    double distance;
    int node_index;
    std::tie(distance, node_index) = queue.top();
    queue.pop();

    if (distance > node_index_to_distance[node_index]) {
      // Lazy element deletion. This predecessor has been re-visited on a
      // shorter path, ignore it.
      continue;
    }

    for (const int edge_index : node_index_to_outgoing_edges[node_index]) {
      const Edge& edge = graph.edges()[edge_index];
      assert(edge.cost >= 0.0 && "Graph has negative edge cost");
      const int successor_index = edge.target_index;
      const double path_cost = distance + edge.cost;

      if (path_cost < node_index_to_distance[successor_index]) {
        node_index_to_distance[successor_index] = path_cost;
        node_index_to_predecessor_edge[successor_index] = edge_index;
        queue.emplace(path_cost, successor_index);
      }
    }
  }
}

#ifdef MCF_USE_Boost
void ShortestPathDijkstraFibonacciHeap(
    const Graph& graph,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    std::vector<double>& node_index_to_distance,
    std::vector<int>& node_index_to_predecessor_edge) {
  // Initialize to "no predecessor" and infinity cost.
  node_index_to_predecessor_edge.resize(graph.num_nodes());
  std::fill_n(node_index_to_predecessor_edge.begin(), graph.num_nodes(), -1);
  node_index_to_distance.resize(graph.num_nodes());
  std::fill_n(node_index_to_distance.begin(), graph.num_nodes(),
              std::numeric_limits<double>::infinity());

  // Set up heap structure.
  using FibonacciHeap = boost::heap::fibonacci_heap<
      std::tuple<double, int>,
      boost::heap::compare<std::greater<std::tuple<double, int>>>>;
  using ElementHandle = FibonacciHeap::handle_type;

  FibonacciHeap queue;
  std::map<int, ElementHandle> node_index_to_handle;

  // Add source to queue.
  auto handle = queue.emplace(0.0, Graph::InternalSourceNode);
  node_index_to_handle[Graph::InternalSourceNode] = handle;

  // Iterate until all elements have been processed.
  while (!queue.empty()) {
    double distance;
    int node_index;
    std::tie(distance, node_index) = queue.top();

    node_index_to_handle.erase(node_index);
    queue.pop();

    for (const int edge_index : node_index_to_outgoing_edges[node_index]) {
      const Edge& edge = graph.edges()[edge_index];
      const int successor_index = edge.target_index;
      const double path_cost = distance + edge.cost;
      assert(path_cost >= 0.0 && "Graph has negative path cost");

      if (path_cost < node_index_to_distance[successor_index]) {
        auto handle_it = node_index_to_handle.find(successor_index);
        if (handle_it != node_index_to_handle.end()) {
          const auto& handle = handle_it->second;
          std::get<0>(*handle) = path_cost;
          queue.decrease(handle);
        } else {
          auto handle = queue.emplace(path_cost, successor_index);
          node_index_to_handle[successor_index] = handle;
        }

        node_index_to_distance[successor_index] = path_cost;
        node_index_to_predecessor_edge[successor_index] = edge_index;
      }
    }
  }
}
#endif  // MCF_USE_Boost

void EliminateNegativeCost(const std::vector<double>& node_index_to_distance,
                           Graph& graph) {
  for (Edge& edge : graph.mutable_edges()) {
    edge.cost += node_index_to_distance[edge.source_index] -
                 node_index_to_distance[edge.target_index];
    edge.cost = std::max(0.0, edge.cost);

    if (edge.cost < kEpsilon) {
      // Due to floating point rounding errors, the cost on the shortest path
      // may not evaluate to exactly 0.0. If this is not handled, in consecutive
      // iterations the residual graph contains negative edge weights (and
      // possibly negative cycles) that prohibit application of Dijstra.
      edge.cost = 0.0;
    }
  }
}

void ReverseEdges(const std::vector<int>& edge_indices, Graph& graph,
                  std::vector<std::vector<int>>& node_index_to_incoming_edges,
                  std::vector<std::vector<int>>& node_index_to_outgoing_edges) {
  std::vector<Edge>& edges = graph.mutable_edges();

  std::vector<bool> has_been_swapped(edges.size(), false);
  for (const int edge_index : edge_indices) {
    Edge& edge = edges[edge_index];
    std::swap(edge.source_index, edge.target_index);
    edge.cost = -edge.cost;

    assert(edge.cost >= 0.0);
    has_been_swapped[edge_index] = true;
  }

  const int num_nodes = node_index_to_incoming_edges.size();
  assert(node_index_to_incoming_edges.size() ==
             node_index_to_outgoing_edges.size() &&
         "node_index_to_incoming_edges.size() != "
         "node_index_to_outgoing_edges.size()");
  for (int node_index = 0; node_index < num_nodes; ++node_index) {
    std::vector<int>& node_index_to_incoming_edges_n =
        node_index_to_incoming_edges[node_index];
    std::vector<int>& node_index_to_outgoing_edges_n =
        node_index_to_outgoing_edges[node_index];
    const std::size_t num_notswapped_node_index_to_outgoing_edges =
        node_index_to_outgoing_edges_n.size();

    // First, remove edges that have changed direction from
    // node_index_to_incoming_edges_n
    // and append them to node_index_to_outgoing_edges_n.
    for (std::size_t i = 0; i < node_index_to_incoming_edges_n.size();) {
      const int edge_index = node_index_to_incoming_edges_n[i];
      if (!has_been_swapped[edge_index]) {
        ++i;
        continue;
      }

      std::swap(node_index_to_incoming_edges_n[i],
                node_index_to_incoming_edges_n.back());
      node_index_to_incoming_edges_n.pop_back();
      node_index_to_outgoing_edges_n.push_back(edge_index);
    }

    // Now, check elements in node_index_to_outgoing_edges_n. Edges with
    // index >= num_notswapped_node_index_to_outgoing_edges are incoming edges
    // that have been swapped in the previous loop, i.e., leave them untouched.
    for (std::size_t i = 0; i < num_notswapped_node_index_to_outgoing_edges &&
                            i < node_index_to_outgoing_edges_n.size();
         ++i) {
      const int edge_index = node_index_to_outgoing_edges_n[i];
      if (!has_been_swapped[edge_index]) {
        continue;
      }

      node_index_to_incoming_edges_n.push_back(edge_index);
      std::swap(node_index_to_outgoing_edges_n[i],
                node_index_to_outgoing_edges_n.back());
      node_index_to_outgoing_edges_n.pop_back();
    }
  }
}

bool FindPathFromSinkToSource(
    const Graph& graph, const std::vector<int>& node_index_to_predecessor_edge,
    std::vector<int>& traversed_edges, double& path_cost) {
  traversed_edges.clear();
  path_cost = 0.0;

  int node_index = Graph::InternalSinkNode;
  while (node_index != Graph::InternalSourceNode) {
    const int edge_index = node_index_to_predecessor_edge[node_index];
    if (edge_index < 0) {
      return false;  // There is no path from T to S.
    }

    const Edge& edge = graph.edges()[edge_index];
    traversed_edges.push_back(edge_index);
    path_cost += edge.cost;
    node_index = edge.source_index;
  }

  return true;
}

void ComputeSuccessorMapFromResidualGraph(
    const Graph& residual_graph,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    std::vector<int>& trajectory_heads,
    std::vector<int>& node_index_to_successor) {
  auto IsReverseEdge = [](const Edge& edge) -> bool {
    if (edge.target_index == Graph::InternalSinkNode) {
      return false;
    }
    if (edge.source_index != Graph::InternalSinkNode &&
        edge.source_index <= edge.target_index) {
      return false;
    }
    return true;
  };

  node_index_to_successor.resize(residual_graph.num_nodes());
  std::fill_n(node_index_to_successor.begin(), residual_graph.num_nodes(), -1);
  trajectory_heads.clear();

  for (const int edge_index :
       node_index_to_outgoing_edges[Graph::InternalSinkNode]) {
    if (!IsReverseEdge(residual_graph.edges()[edge_index])) {
      continue;
    }

    int node_index = residual_graph.edges()[edge_index].target_index;
    node_index_to_successor[node_index] = Graph::InternalSinkNode;
    while (node_index != Graph::InternalSourceNode) {
      int next_node_index = -1;
      for (const int edge_index : node_index_to_outgoing_edges[node_index]) {
        const Edge& edge = residual_graph.edges()[edge_index];
        if (!IsReverseEdge(edge)) {
          continue;
        }
        next_node_index = edge.target_index;
        break;  // There is always exactly one outgoing reverse edge.
      }
      assert(next_node_index >= 0 && "Disconnected T-S path.");
      if (next_node_index == Graph::InternalSourceNode) {
        trajectory_heads.push_back(node_index);
      } else {
        node_index_to_successor[next_node_index] = node_index;
      }
    }
  }
}

double ComputeTrajectoriesFromResidualGraph(
    const Graph& residual_graph,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    std::vector<std::vector<int>>& trajectories,
    std::vector<int>* node_index_to_shortest_path_incoming_edge) {
  auto IsReverseEdge = [](const Edge& edge) -> bool {
    if (edge.target_index == Graph::InternalSinkNode) {
      return false;
    }
    if (edge.source_index != Graph::InternalSinkNode &&
        edge.source_index <= edge.target_index) {
      return false;
    }
    return true;
  };

  if (node_index_to_shortest_path_incoming_edge != nullptr) {
    node_index_to_shortest_path_incoming_edge->resize(
        residual_graph.num_nodes());
    std::fill_n(node_index_to_shortest_path_incoming_edge->begin(),
                residual_graph.num_nodes(), -1);
  }

  double total_cost = 0.0;
  trajectories.clear();
  for (const int edge_index :
       node_index_to_outgoing_edges[Graph::InternalSinkNode]) {
    if (!IsReverseEdge(residual_graph.edges()[edge_index])) {
      continue;
    }

    std::vector<int> trajectory;
    double cost_sum = 0.0;
    int node_index = residual_graph.edges()[edge_index].target_index;
    cost_sum += residual_graph.edges()[edge_index].cost;

    while (node_index != Graph::InternalSourceNode) {
      if (node_index % 2 == 0) {  // Every second edge is an observation edge.
        trajectory.push_back(node_index / 2);  // Convert to location index.
      }

      int next_node_index = -1;
      int next_incoming_edge_index = -1;
      for (const int edge_index : node_index_to_outgoing_edges[node_index]) {
        const Edge& edge = residual_graph.edges()[edge_index];
        if (!IsReverseEdge(edge)) {
          continue;
        }
        next_node_index = edge.target_index;
        next_incoming_edge_index = edge_index;
        cost_sum += edge.cost;
        break;  // There is always exactly one outgoing reverse edge.
      }
      assert(next_node_index >= 0 && "Disconnected T-S path.");
      if (node_index_to_shortest_path_incoming_edge != nullptr) {
        (*node_index_to_shortest_path_incoming_edge)[node_index] =
            next_incoming_edge_index;
      }
      node_index = next_node_index;
    }

    std::reverse(trajectory.begin(), trajectory.end());
    trajectories.push_back(std::move(trajectory));
    total_cost += cost_sum;
  }

  return total_cost;
}

int RunSuccessiveShortestPathSearchDijkstraLazyDeletion(
    Graph& residual_graph,
    std::vector<std::vector<int>>& node_index_to_incoming_edges,
    std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    const int max_flow, const bool stop_on_best) {
  return RunSuccessiveShortestPathSearch(
      ShortestPathDijkstraLazyDeletion, residual_graph,
      node_index_to_incoming_edges, node_index_to_outgoing_edges, max_flow,
      stop_on_best);
}

#ifdef MCF_USE_Boost
int RunSuccessiveShortestPathSearchDijkstraFibonacciHeap(
    Graph& residual_graph,
    std::vector<std::vector<int>>& node_index_to_incoming_edges,
    std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    const int max_flow, const bool stop_on_best) {
  return RunSuccessiveShortestPathSearch(
      ShortestPathDijkstraFibonacciHeap, residual_graph,
      node_index_to_incoming_edges, node_index_to_outgoing_edges, max_flow,
      stop_on_best);
}
#endif

}  // namespace internal

}  // namespace mcf
