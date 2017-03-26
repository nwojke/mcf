// vim: expandtab:ts=2:sw=2
#ifndef MCF_UTIL_HPP
#define MCF_UTIL_HPP

#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <mcf/graph.hpp>

#define ASSERT_NOTNULL(pointer)                                        \
  if ((pointer) == nullptr) {                                          \
    std::stringstream msg;                                             \
    msg << "Null pointer in " << __func__ << " at " << __FILE__ << ":" \
        << __LINE__;                                                   \
    throw std::runtime_error(msg.str().c_str());                       \
  }

#define ASSERT_CONTAINER_SIZE(expected_size, container)                    \
  if (static_cast<decltype(expected_size)>((container).size()) !=          \
      (expected_size)) {                                                   \
    std::stringstream msg;                                                 \
    msg << "Container size mismatch in " << __func__ << " at " << __FILE__ \
        << ":" << __LINE__;                                                \
    throw std::runtime_error(msg.str().c_str());                           \
  }

#define ASSERT_TRUE(expression, failure_msg) \
  if (!(expression)) {                       \
    throw std::runtime_error(failure_msg);   \
  }

namespace mcf {

namespace internal {

/**
 * Compute incoming and outgoing edges to and from every node in a graph.
 *
 * @param graph The graph structure.
 * @param node_index_to_incoming_edges Contains the computed incoming edge
 *        indices for every node in the graph.
 * @param node_index_to_outgoing_edges Contains the computed outgoing edege
 *        indices for every node in the graph.
 */
void BuildEdgeMap(const Graph& graph,
                  std::vector<std::vector<int>>& node_index_to_incoming_edges,
                  std::vector<std::vector<int>>& node_index_to_outgoing_edges);

/**
 * Compute trajectories from successor map.
 *
 * @param trajectory_heads Trajectory entry points (node indices).
 * @param node_index_to_successor Maps each node index to the index of the
 *        subsequent node on the trajectory. Each trajectory must end with
 *        Graph::InternalSinkNode and nodes that are not part of any trajectory
 *        should be set to -1.
 * @return trajectories A list of trajectories, where each trajectory is a
 *         sequence of location handles  (Graph::ST is not part of the
 *         trajectory).
 */
std::vector<std::vector<int>> ComputeTrajectories(
    const std::vector<int>& trajectory_heads,
    const std::vector<int>& node_index_to_successor);

/**
 * Find the global minimum of a convex function in a given
 * range [min_bound, max_bound].
 *
 * @param f A callable that takes a single argument of type int to evaluate a
 *        convex function and returns the function value at the given input.
 * @param min_bound Lower bound of the range.
 * @param max_bound Upper bound of the range.
 * @param comparator A callable that takes two function values (of return type
 *        of f) and returns true if the first argument is smaller than the
 *        second.
 * @return The minimum in range [min_bound, max_bound].
 */
template <typename F, typename Compare = std::less<std::result_of_t<F(int)>>>
std::result_of_t<F(int)> BinarySearch(F f, int min_bound, int max_bound,
                                      Compare comparator = Compare());

}  // namespace internal

}  // namespace mcf

#include <mcf/impl/util_inl.hpp>

#endif
