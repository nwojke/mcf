// vim: expandtab:ts=2:sw=2
#include <mcf/internal/util.hpp>

#include <cassert>

namespace mcf {

namespace internal {

void BuildEdgeMap(const Graph& graph,
                  std::vector<std::vector<int>>& node_index_to_incoming_edges,
                  std::vector<std::vector<int>>& node_index_to_outgoing_edges) {
  const int num_nodes = graph.num_nodes();
  const int num_edges = graph.edges().size();

  node_index_to_incoming_edges.resize(num_nodes);
  node_index_to_outgoing_edges.resize(num_nodes);
  for (int edge_index = 0; edge_index < num_edges; ++edge_index) {
    const Edge& edge = graph.edges()[edge_index];
    node_index_to_outgoing_edges[edge.source_index].push_back(edge_index);
    node_index_to_incoming_edges[edge.target_index].push_back(edge_index);
  }
}

std::vector<std::vector<int>> ComputeTrajectories(
    const std::vector<int>& trajectory_heads,
    const std::vector<int>& node_index_to_successor) {
  std::vector<std::vector<int>> trajectories(trajectory_heads.size());
  for (std::size_t t = 0; t < trajectory_heads.size(); ++t) {
    int v = trajectory_heads[t];
    while (v != Graph::InternalSinkNode) {
      assert(node_index_to_successor.at(v) >= 0 &&
             "inconsistent flow detected");
      if (v % 2 == 0) {  // Convert to location index
        trajectories[t].push_back(v / 2);
      }
      v = node_index_to_successor.at(v);
    }
  }
  return trajectories;
}

}  // namespace internal

}  // namespace mcf
