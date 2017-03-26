// vim: expandtab:ts=2:sw=2
#ifndef MCF_K_SHORTEST_PATH_HPP
#define MCF_K_SHORTEST_PATH_HPP

#include <limits>
#include <mcf/graph.hpp>
#include <vector>

namespace mcf {

namespace internal {

/**
 * Check that a given graph is ordered such that no node has a larger index
 * than any of its successors, with Graph::InternalSinkNode being the only
 * exception.
 *
 * @param graph The graph to be checked.
 * @return True if the graph is ordered, false otherwise.
 */
bool CheckGraphIsOrdered(const Graph& graph);

/**
 * Compute the shortest path in an ordered directed acyclic graph.
 *
 * @param graph A directed acyclic graph that is ordered such that no node has
 *        a larger index than any of its successors (Graph::InternalSinkNode is
 *        the only allowed exception).
 * @param node_index_to_incoming_edges A vector that contains for every node
 *        the list of incoming edge indices.
 * @param node_index_to_predecessor_edge Contains the computed index of the
 *        incoming edge along shortest from the source to every node.
 * @param node_index_to_distance Contains the computed shortest distance
 *        between the source and every node.
 */
void ShortestPathInOrderedDAG(
    const Graph& graph,
    const std::vector<std::vector<int>>& node_index_to_incoming_edges,
    std::vector<double>& node_index_to_distance,
    std::vector<int>& node_index_to_predecessor_edge);

/**
 * Compute the shortest path in a directed graph with non-negative edge costs.
 *
 * This implementation uses an immutable std::priority_queue for the heap
 * structure.
 *
 * @param graph A directed graph with non-negative edge costs.
 * @param node_index_to_outgoing_edges A vector that contains for every node
 *        the list of outgoing edge indices.
 * @param node_index_to_predecessor_edge Contains the computed index of the
 *        incoming edge along shortest from the source to every node.
 * @param node_index_to_distance Contains the computed shortest distance
 *        between the source and every node.
 */
void ShortestPathDijkstraLazyDeletion(
    const Graph& graph,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    std::vector<double>& node_index_to_distance,
    std::vector<int>& node_index_to_predecessor_edge);

#ifdef MCF_USE_Boost
/**
 * Compute the shortest path in a directed graph with non-negative edge costs.
 *
 * This implementation uses a Fibonacci heap with decrease-key functionality.
 *
 * @param graph A directed graph with non-negative edge costs.
 * @param node_index_to_outgoing_edges A vector that contains for every node
 *        the list of outgoing edge indices.
 * @param node_index_to_predecessor_edge Contains the computed index of the
 *        incoming edge along shortest from the source to every node.
 * @param node_index_to_distance Contains the computed shortest distance
 *        between the source and every node.
 */
void ShortestPathDijkstraFibonacciHeap(
    const Graph& graph,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    std::vector<double>& node_index_to_distance,
    std::vector<int>& node_index_to_predecessor_edge);
#endif

/**
 * Transform a directed acyclic graph with arbitrary (possibly negative) edge
 * costs to a graph with positive costs only.
 *
 * For every edge (i, j) from node i to node j, the transformed cost is
 * computed as follows:
 *
 *     c_{i,j}' = c_{i,j} + d[i] - d[v],
 *
 * where d[i] contains the shortest path distance from the source to node i.
 * Consequently, this transformation creates edges with 0 cost along the
 * shortest paths.
 *
 * @param node_index_to_distance A vector containing the shortest distance from
 *        the source Graph::InternalSourceNode to every node.
 * @param graph A graph with arbitrary edge cost that, after transformation,
 *        contains positive cost terms only.
 */
void EliminateNegativeCost(const std::vector<double>& node_index_to_distance,
                           Graph& graph);

/**
 * Reverse direction of edges and change the sign of their associated costs.
 *
 * @param edge_indices Indices to edges to revert.
 * @param graph The graph structure.
 * @param node_index_to_incoming_edges A vector that contains for every node
 *        the list of incoming edge indices.
 * @param node_index_to_outgoing_edges A vector that contains for every node
 *        the list of outgoing edge indices.
 */
void ReverseEdges(const std::vector<int>& edge_indices, Graph& graph,
                  std::vector<std::vector<int>>& node_index_to_incoming_edges,
                  std::vector<std::vector<int>>& node_index_to_outgoing_edges);

/**
 * Traverse along predecessor map to extract a path from sink to source.
 *
 * @param graph The graph structure.
 * @param node_index_to_predecessor_edge The predecessor map contains the index
 *        of the incoming edge along the shortest path from the source to
 *        every node.
 * @param traversed_edges Contains the edges on the path in order from sink
 *        to source.
 * @param path_cost The sum of edge costs along the path.
 * @return True on success, false if there is no path from
 *         Graph::InternalSinkNode to Graph::InternalSourceNode.
 */
bool FindPathFromSinkToSource(
    const Graph& graph, const std::vector<int>& node_index_to_predecessor_edge,
    std::vector<int>& traversed_edges, double& path_cost);

/**
 * Run k shortest path search using ShortestPathDijkstraLazyDeletion as
 * shortest path solver.
 *
 * @param shortest_path_finder A shortest path solver function (must be any of
 *        ShortestPathDijkstra*).
 * @param residual_graph The residual graph to operate on. At input, the graph
 *        must be ordered such that no node has a larger index than any of its
 *        successors (Graph::InternalSinkNode is an exception). This graph is
 *        manipulated as part of the search.
 *        After completion, ComputeTrajectories can be called to extract
 *        trajectories from the final residual graph.
 * @param node_index_to_incoming_edges A vector that contains for every node in
 *        the residual graph the list of incoming edge indices. This vector is
 *        manipulated as the direction of edges changes as part of the search.
 * @param node_index_to_outgoing_edges A vector that contains for every node in
 *        the residual graph the list of outgoing edge indices. This vector is
 *        manipulated as the direction of edges changes as part of the search.
 * @param max_flow An upper bound for the total flow (number of trajectories).
 *        The shortest path search is stopped when this value is reached.
 * @param stop_on_best If true, search is stopped when the global best solution
 *        has been found. If false, search continues until max_flow is reached
 *        or until there exists no more feasible solution.
 * @return The final flow (number of trajectories). Should be max_flow if
 *         stop_on_best is false, otherwise the solution is infeasible.
 *
 * @throws std::runtime_error If the input graph is not ordered.
 */
int RunSuccessiveShortestPathSearchDijkstraLazyDeletion(
    Graph& residual_graph,
    std::vector<std::vector<int>>& node_index_to_incoming_edges,
    std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    int max_flow = std::numeric_limits<int>::max(), bool stop_on_best = true);

#ifdef MCF_USE_Boost
/**
 * Run k shortest path search using ShortestPathDijkstraFibonacciHeap as
 * shortest path solver.
 *
 * @param shortest_path_finder A shortest path solver function (must be any of
 *        ShortestPathDijkstra*).
 * @param residual_graph The residual graph to operate on. At input, the graph
 *        must be ordered such that no node has a larger index than any of its
 *        successors (Graph::InternalSinkNode is an exception). This graph is
 *        manipulated as part of the search.
 *        After completion, ComputeTrajectories can be called to extract
 *        trajectories from the final residual graph.
 * @param node_index_to_incoming_edges A vector that contains for every node in
 *        the residual graph the list of incoming edge indices. This vector is
 *        manipulated as the direction of edges changes as part of the search.
 * @param node_index_to_outgoing_edges A vector that contains for every node in
 *        the residual graph the list of outgoing edge indices. This vector is
 *        manipulated as the direction of edges changes as part of the search.
 * @param max_flow An upper bound for the total flow (number of trajectories).
 *        The shortest path search is stopped when this value is reached.
 * @param stop_on_best If true, search is stopped when the global best solution
 *        has been found. If false, search continues until max_flow is reached
 *        or until there exists no more feasible solution.
 * @return The final flow (number of trajectories). Should be max_flow if
 *         stop_on_best is false, otherwise the solution is infeasible.
 *
 * @throws std::runtime_error If the input graph is not ordered.
 */
int RunSuccessiveShortestPathSearchDijkstraFibonacciHeap(
    Graph& residual_graph,
    std::vector<std::vector<int>>& node_index_to_incoming_edges,
    std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    int max_flow = std::numeric_limits<int>::max(), bool stop_on_best = true);
#endif

/**
 * Compute trajectories from residual graph. This function traverses reversed
 * edges from sink to source to compute the set of trajectories.
 *
 * @param residual_graph The residual graph structure.
 * @param node_index_to_outgoing_edges Outgoing edge indices for every node
 *        in the residual graph.
 * @param trajectories On completion, contains the computed trajectories, where
 *        each trajectory is a sequence of location handles (Graph::ST is not
 *        part of the trajectory).
 * @param node_index_to_shortest_path_incoming_edge If not nullptr, contains
 *        for every the index of the incoming edge on the shortest path in the
 *        original graph (not residual graph). Entries for nodes not on a
 *        shortest path are set to -1.
 * @return Returns the sum of trajectory costs.
 */
double ComputeTrajectoriesFromResidualGraph(
    const Graph& residual_graph,
    const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
    std::vector<std::vector<int>>& trajectories,
    std::vector<int>* node_index_to_shortest_path_incoming_edge = nullptr);

}  // namespace internal

}  // namespace mcf

#endif
