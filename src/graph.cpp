// vim: expandtab:ts=2:sw=2
#include <mcf/graph.hpp>

#include <limits>

namespace mcf {

const int Graph::InternalSourceNode = 1;

const int Graph::InternalSinkNode = 0;

const int Graph::FirstNonSourceSinkNode = 2;

const int Graph::ST = 0;

Graph::Graph() : next_id_(2) {}

void Graph::Reserve(int num_edges) { edges_.reserve(num_edges); }

int Graph::Add(double cost) {
  int node_id = next_id_;
  next_id_ += 2;
  edges_.push_back({node_id, node_id + 1, cost});
  return node_id / 2;
}

void Graph::Link(int src, int dst, double cost) {
  src = 2 * src + 1;  // This maps ST to 1, i.e., InternalSourceNode.
  dst = 2 * dst + 0;  // This maps ST to 0, i.e., InternalSinkNode.
  edges_.push_back({src, dst, cost});
}

}  // namespace mcf
