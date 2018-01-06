// vim: expandtab:ts=2:sw=2
#include <mcf/batch_processing.hpp>
#include <mcf/internal/k_shortest_path.hpp>
#include <mcf/internal/util.hpp>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <stdexcept>

namespace mcf {

namespace {

/**
 * Find k shortest paths in graph.
 *
 * @param solver_type The shortest path algorithm to use during the successive
 *        shortest path search.
 * @param graph The graph structure.
 * @param node_index_to_incoming_edges A vector that contains for every node
 *        the list of incoming edge indices.
 * @param node_index_to_outgoing_edges A vector that contanis for every node
 *        the list of outgoing edge indices.
 * @param locations_in_last_timestep The list of locations in the last time
 *        step (only used if ignore_last_exit_cost is true).
 * @param trajectories Contains the computed trajectories on success.
 * @param ignore_last_exit_cost If true, sets the exit cost for locations in
 *        the last time step to 0 prior to calling the solver.
 */
void Solve(ShortestPathSolverType solver_type, const Graph& graph,
           const std::vector<std::vector<int>>& node_index_to_incoming_edges,
           const std::vector<std::vector<int>>& node_index_to_outgoing_edges,
           const std::vector<int> locations_in_last_timestep,
           std::vector<std::vector<int>>& trajectories,
           std::vector<int>& node_index_to_shortest_path_incoming_edge,
           const bool ignore_last_exit_cost) {
  // Copy graph to residual graph structure to be modified by k shortest path
  // algorithm.
  Graph residual_graph = graph;
  std::vector<std::vector<int>> residual_node_index_to_incoming_edges =
      node_index_to_incoming_edges;
  std::vector<std::vector<int>> residual_node_index_to_outgoing_edges =
      node_index_to_outgoing_edges;

  if (ignore_last_exit_cost) {
    for (const int location : locations_in_last_timestep) {
      const int exit_node_index = 2 * location + 1;
      assert(exit_node_index >= Graph::FirstNonSourceSinkNode &&
             exit_node_index < residual_graph.num_nodes());
      for (const int edge_index :
           residual_node_index_to_outgoing_edges[exit_node_index]) {
        Edge& edge = residual_graph.mutable_edges()[edge_index];
        assert(edge.target_index == Graph::InternalSinkNode &&
               "Location in last time step has non-sink outgoing edge.");
        edge.cost = 0.0;
      }
    }
  }

  // Call k shortest path algorithm.
  switch (solver_type) {
    case ShortestPathSolverType::kDijkstraLazyDeletion:
      internal::RunSuccessiveShortestPathSearchDijkstraLazyDeletion(
          residual_graph, residual_node_index_to_incoming_edges,
          residual_node_index_to_outgoing_edges);
      break;
#ifdef MCF_USE_Boost
    case ShortestPathSolverType::kDijkstraFibonacciHeap:
      internal::RunSuccessiveShortestPathSearchDijkstraFibonacciHeap(
          residual_graph, residual_node_index_to_incoming_edges,
          residual_node_index_to_outgoing_edges);
      break;
#endif
  }

  // Extract succesor-map and trajectories.
  internal::ComputeTrajectoriesFromResidualGraph(
      residual_graph, residual_node_index_to_outgoing_edges, trajectories,
      &node_index_to_shortest_path_incoming_edge);
}

/**
 * Check if new_trajectory extends cached_trajectory.
 *
 * @param location_to_timestep A functor that maps from location to time step.
 * @param first_optimized_timestep The first time step in the current
 *        optimization window.
 * @param cached_trajectory A cached trajectory from the previous solver run.
 * @param new_trajectory A new trajectory found at the current solver run.
 * @param matched_cached_index Set to the index of the matched location on
 *        cached_trajectory.
 * @param new_trajectory Set to the index of the matched location on
 *        new_trajectory.
 * @return True if new_trajectory starts on cached_trajectory.
 */
template <typename LocationToTimestep>
bool FindMatchingLocation(
    LocationToTimestep location_to_timestep,
    const BatchProcessing::Index& first_optimized_timestep,
    const std::vector<BatchProcessing::Index>& cached_trajectory,
    const std::vector<BatchProcessing::Index>& new_trajectory,
    int& matched_cached_index, int& matched_new_index) {
  if (cached_trajectory.empty() || new_trajectory.empty()) {
    // Empty trajectory has nothing to match against.
    return false;
  }

  // A valid match is only found if the new_trajectory starts on
  // cached_trajectory.
  const BatchProcessing::Index matched_location = new_trajectory.front();

  // Find the first location in cached_trajectory within the current
  // optimization window.
  int matched_index = -1;
  for (int i = static_cast<int>(cached_trajectory.size()) - 1; i >= 0; --i) {
    if (location_to_timestep(cached_trajectory[i]) < first_optimized_timestep) {
      // We have moved beyond the optimization window.
      break;
    }
    matched_index = static_cast<int>(i);
  }
  if (matched_index < 0) {
    // Trajectory does not extend into the current optimization window.
    return false;
  }

  // Establish a match if the first location inside the current optimization
  // window equals new_trajectory.front().
  const bool locations_match =
      cached_trajectory[matched_index] == matched_location;
  if (locations_match) {
    matched_cached_index = matched_index;
    matched_new_index = 0;
  }
  return locations_match;
}

/**
 * Extend cached trajectory with a matched new_trajectory.
 *
 *param matched_new_index Index of the matching location in new_trajectory.
 * @param new_trajectory A new trajectory found at the current solver run.
 * @param matched_cached_index Index of the matching location in
 *        cached_trajectory.
 * @param cached_trajectory A cached trajectory from previous solver runs; will
 *        be extended with new_trajectory[matched_new_index:], i.e.,
 *        cached_trajectory' = cached_trajectory[:matched_cached_index] +
 *        new_trajectory[matched_new_index:].
 */
void ExtendTrajectory(int matched_new_index,
                      const std::vector<BatchProcessing::Index>& new_trajectory,
                      int matched_cached_index,
                      BatchProcessing::Trajectory& cached_trajectory) {
  assert(matched_cached_index >= 0 && matched_new_index >= 0);
  cached_trajectory.reserve(matched_cached_index + new_trajectory.size() -
                            matched_new_index);
  cached_trajectory.resize(matched_cached_index);
  cached_trajectory.insert(cached_trajectory.end(),
                           new_trajectory.begin() + matched_new_index,
                           new_trajectory.end());
}

/**
 * Merge cached trajectories and newly computed trajectories.
 *
 * @param location_to_timestep A functor that maps from location to time step.
 * @param first_optimized_timestep The first time step in the current
 *        optimization window.
 * @param new_trajectories The new trajectories, computed at the most recent
 *        time step.
 * @param cached_trajectories The cached trajectories, computed in previous
 *        time steps. The cached_trajectories are extended by matches in
 *        new_trajectories.
 * @param active Set to true if the trajectory is part of the current
 *        solution and false otherwise.
 * @param next_trajectory_index A counter used to create unique trajectory
 *        indices.
 */
template <typename LocationToTimestep>
void MergeTrajectories(
    LocationToTimestep location_to_timestep,
    const BatchProcessing::Index& first_optimized_timestep,
    const BatchProcessing::Index& clipping_timestep,
    const std::vector<std::vector<BatchProcessing::Index>>& new_trajectories,
    BatchProcessing::TrajectoryMap& cached_trajectories,
    std::unordered_map<BatchProcessing::Index, bool>& active,
    BatchProcessing::Index& next_trajectory_index) {
  // First, try to match existing trajectories.
  std::vector<bool> new_trajectory_matched(new_trajectories.size(), false);

  active.clear();

  for (auto it = cached_trajectories.begin(); it != cached_trajectories.end();
       ++it) {
    int matched_cached_index = -1;
    for (std::size_t k = 0; k < new_trajectories.size(); ++k) {
      int matched_new_index = -1;
      if (new_trajectory_matched[k] ||
          !FindMatchingLocation(location_to_timestep, first_optimized_timestep,
                                it->second, new_trajectories[k],
                                matched_cached_index, matched_new_index)) {
        // If a trajectory has been matched in an earlier iteration, assert
        // that we are not skipping on a potential assignment.
        assert(!FindMatchingLocation(
            location_to_timestep, first_optimized_timestep, it->second,
            new_trajectories[k], matched_cached_index, matched_new_index));
        continue;
      }

      ExtendTrajectory(matched_new_index, new_trajectories[k],
                       matched_cached_index, it->second);
      new_trajectory_matched[k] = true;
    }

    // Mark the trajectory active if (1) it has been matched against a
    // new_trajectory or (2) it is entirely outside of the current optimization
    // window.
    active[it->first] =
        matched_cached_index >= 0 ||
        location_to_timestep(it->second.back()) < first_optimized_timestep;
  }

  // Remove any trajectory from the cache which is only partially inside the
  // optimization window but is not currently detected. Simply setting
  // active=false can lead to ambiguities in MergeTrajectories() in future time
  // steps.
  for (auto it = cached_trajectories.begin();
       it != cached_trajectories.end();) {
    if (!active[it->first] &&
        location_to_timestep(it->second.front()) < clipping_timestep) {
      active.erase(it->first);
      it = cached_trajectories.erase(it);
    } else {
      ++it;
    }
  }

  // Now, add all unmatched trajectories to the cache.
  for (std::size_t i = 0; i < new_trajectories.size(); ++i) {
    if (new_trajectory_matched[i]) {
      continue;
    }
    const BatchProcessing::Index index = next_trajectory_index++;
    cached_trajectories[index] = new_trajectories[i];
    active[index] = true;
  }
}

/**
 * Collapse a trajectory into the trajectory head.
 *
 * This function traverses the shortest path from Graph::InternalSourceNode to
 * the entry node of trajectories[i][clipping_index[i]] to sum up the
 * trajectory cost. Then, it sets the entry cost of this node to the computed
 * cost sum and sets traversed edge costs to infinity.
 *
 * @param trajectories Trajectories of locations in the graph.
 * @param clipping_indices Indices of new trajectory heads.
 * @param node_index_to_incoming_edges A vector that contains for every node
 *        the list of incoming edge indices.
 * @param node_index_to_shorest_path_incoming_edge A vector that contains for
 *        every node the index of the incoming edge on the shortest path
 *        starting at Graph::InternalSourceNode.
 * @param The graph structure from which trajectories have been obtained; will
 *        be modified.
 */
void CollapseTrajectories(
    const std::vector<std::vector<int>>& trajectories,
    const std::vector<int>& clipping_indices,
    const std::vector<std::vector<int>>& node_index_to_incoming_edges,
    const std::vector<int>& node_index_to_shortest_path_incoming_edge,
    Graph& graph) {
  assert(static_cast<int>(node_index_to_incoming_edges.size()) ==
         graph.num_nodes());
  assert(static_cast<int>(node_index_to_shortest_path_incoming_edge.size()) ==
         graph.num_nodes());

  for (std::size_t trajectory_index = 0; trajectory_index < trajectories.size();
       ++trajectory_index) {
    const std::vector<int>& trajectory = trajectories[trajectory_index];
    const int head_index = clipping_indices[trajectory_index];
    if (head_index < 0) {
      // This trajectory starts inside the new optimization window, no need
      // to prune.
      continue;
    }
    assert(head_index < static_cast<int>(trajectory.size()));

    const int head_entry_node_index = 2 * trajectory[head_index];
    int node_index = head_entry_node_index;
    double cost_sum = 0.0;

    // Sum trajectory cost and add to new head entry node.
    while (node_index != Graph::InternalSourceNode) {
      assert(node_index >= Graph::FirstNonSourceSinkNode &&
             node_index < graph.num_nodes());
      const int edge_index =
          node_index_to_shortest_path_incoming_edge[node_index];
      assert(edge_index >= 0 &&
             edge_index < static_cast<int>(graph.edges().size()));

      Edge& edge = graph.mutable_edges()[edge_index];
      cost_sum += edge.cost;
      if (edge.source_index >= Graph::InternalSourceNode) {
        // Setting all edges to infinity cost, any predecessor of the new
        // head cannot be part of any trajectory in following solver steps.
        edge.cost = std::numeric_limits<double>::infinity();
      }
      node_index = edge.source_index;
    }

    // Set new head entry cost to trajectory cost. Set cost of remaining
    // incoming edges to infinity.
    for (const int edge_index :
         node_index_to_incoming_edges[head_entry_node_index]) {
      assert(edge_index >= 0 &&
             edge_index < static_cast<int>(graph.edges().size()));
      Edge& edge = graph.mutable_edges()[edge_index];
      edge.cost = edge.source_index == Graph::InternalSourceNode
                      ? cost_sum
                      : std::numeric_limits<double>::infinity();
    }
  }
}

/**
 * Remove inactive trajectories from the cache when they leave the optimization
 * window.
 *
 * @param location_to_timestep A functor that maps from location to time step.
 * @param clipping_timestep The first time step in the new optimization window.
 *        Inactive trajectories with starting time step smaller than this value
 *        are removed from the cache.
 * @param cached_trajectories A sparse trajectory cache that maps from index
 *        to trajectory.
 * @param active A map that indicates whether the corresponding entry in
 *        cached_trajectories is part of the current solution or if it
 *        corresponds to a cached, but inactive previous solution.
 */
template <typename LocationToTimestep>
void CleanupTrajectoryCache(
    LocationToTimestep location_to_timestep,
    const BatchProcessing::Index& clipping_timestep,
    BatchProcessing::TrajectoryMap& cached_trajectories,
    std::unordered_map<BatchProcessing::Index, bool>& active) {
  for (auto it = cached_trajectories.begin();
       it != cached_trajectories.end();) {
    assert(active.count(it->first) > 0);
    assert(!it->second.empty());

    if (active[it->first] ||
        location_to_timestep(it->second.front()) < clipping_timestep) {
      ++it;
      continue;
    }

    it = cached_trajectories.erase(it);
  }
}

/**
 * Pop first n non-source/sink nodes from the graph.
 *
 * @param n Number of nodes to remove.
 * @param graph The graph.
 */
void PopFirstN(int n, Graph& graph) {
  std::vector<Edge>& edges = graph.mutable_edges();
  const int num_nodes = graph.num_nodes();
  assert(num_nodes >= Graph::FirstNonSourceSinkNode + n &&
         "Graph contains less than n nodes");

  // Erase all edges pointing to or originating from the first n nodes.
  auto IsInFirstN = [n](const int x) {
    // Assume we have nodes [0, 1, 2, 3, 4, 5] and we want to prune n=2 nodes.
    // Then, the smallest to prune is Graph::FirstNonSourceSinkNode=2 and the
    // largest to prune is Graph::FirstNonSinkNode + n = 4, such that the
    // surviving nodes are [0, 1, 4, 5].
    if (x < Graph::FirstNonSourceSinkNode) {
      return false;
    }
    if (x >= Graph::FirstNonSourceSinkNode + n) {
      return false;
    }
    return true;
  };

  auto EdgePointsToFromFirstN = [&IsInFirstN](const Edge& edge) {
    return IsInFirstN(edge.source_index) || IsInFirstN(edge.target_index);
  };

  edges.erase(
      std::remove_if(edges.begin(), edges.end(), EdgePointsToFromFirstN),
      edges.end());

  // Decrement indices of remaining edges and overwrite number of nodes.
  for (Edge& edge : edges) {
    if (edge.source_index >= Graph::FirstNonSourceSinkNode) {
      edge.source_index -= n;
    }
    if (edge.target_index >= Graph::FirstNonSourceSinkNode) {
      edge.target_index -= n;
    }
  }

  graph.overwrite_num_nodes(num_nodes - n);
}

}  // namespace

const BatchProcessing::Index BatchProcessing::ST =
    static_cast<BatchProcessing::Index>(Graph::ST);

BatchProcessing::BatchProcessing(int window_len,
                                 ShortestPathSolverType solver_type)
    : solver_type_(solver_type),
      window_len_(window_len),
      current_timestep_(0),
      previous_clipping_timestep_(0),
      num_pruned_locations_(0),
      next_trajectory_index_(0) {
  location_to_timestep_[BatchProcessing::ST] =
      std::numeric_limits<Index>::max();
  timestep_to_locations_.emplace(0, std::vector<Index>());
}

void BatchProcessing::Reserve(int num_edges) { graph_.Reserve(num_edges); }

BatchProcessing::Index BatchProcessing::Add(double cost) {
  const int graph_index = graph_.Add(cost);
  const Index location = to_sequence_index(graph_index);

  location_to_timestep_[location] = current_timestep_;
  timestep_to_locations_[current_timestep_].push_back(location);
  return location;
}

void BatchProcessing::Link(Index src, Index dst, double cost) {
  if ((src != static_cast<Index>(ST) &&
       static_cast<Index>(src) < num_pruned_locations_) ||
      (dst != static_cast<Index>(ST) &&
       static_cast<Index>(dst) < num_pruned_locations_)) {
    std::stringstream msg;
    msg << "Cannot link a pruned location. Source: " << src
        << " target: " << dst
        << " first non-pruned location: " << num_pruned_locations_;
    throw std::invalid_argument(msg.str().c_str());
  }

  const int graph_src = to_graph_index(src);
  const int graph_dst = to_graph_index(dst);
  assert(graph_src >= 0 && graph_src <= graph_.num_nodes());
  assert(graph_dst >= 0 && graph_dst <= graph_.num_nodes());
  graph_.Link(graph_src, graph_dst, cost);
}

void BatchProcessing::FinalizeTimeStep() { ++current_timestep_; }

void BatchProcessing::RunSearch(std::vector<Trajectory>& trajectories,
                                bool ignore_last_exit_cost) {
  Update(ignore_last_exit_cost);

  trajectories.resize(next_trajectory_index_);
  for (const auto& index_and_trajectory : trajectories_) {
    assert(active_.count(index_and_trajectory.first) > 0);
    if (!active_.find(index_and_trajectory.first)->second) {
      continue;
    }

    Trajectory& trajectory = trajectories.at(index_and_trajectory.first);
    trajectory.reserve(index_and_trajectory.second.size());
    trajectory.insert(trajectory.end(), index_and_trajectory.second.begin(),
                      index_and_trajectory.second.end());
  }
}

BatchProcessing::TrajectoryMap BatchProcessing::ComputeTrajectories(
    bool ignore_last_exit_cost) {
  Update(ignore_last_exit_cost);

  TrajectoryMap trajectories;
  for (auto it = trajectories_.begin(); it != trajectories_.end(); ++it) {
    assert(active_.count(it->first) > 0);
    if (!active_.find(it->first)->second) {
      continue;
    }
    trajectories.emplace(it->first, it->second);
  }

  return trajectories;
}

void BatchProcessing::RemoveInactiveTracks() {
  for (auto it = trajectories_.begin(); it != trajectories_.end();) {
    assert(!it->second.empty());
    if (location_to_timestep_[it->second.back()] <
        previous_clipping_timestep_) {
      it = trajectories_.erase(it);
      continue;
    }
    ++it;
  }
}

void BatchProcessing::Update(bool ignore_last_exit_cost) {
  const Index num_timesteps_in_graph =
      current_timestep_ - previous_clipping_timestep_;
  print("Batch processor called at time step ", current_timestep_);
  print("Number of time steps in graph: ", num_timesteps_in_graph);
  if (num_timesteps_in_graph <= 0) {
    // No new data in the graph, no need to solve.
    return;
  }

  // Run solver.
  std::vector<std::vector<int>> node_index_to_incoming_edges;
  std::vector<std::vector<int>> node_index_to_outgoing_edges;
  internal::BuildEdgeMap(graph_, node_index_to_incoming_edges,
                         node_index_to_outgoing_edges);

  const auto& locations_in_last_timestep =
      timestep_to_locations_[current_timestep_ - 1];  // the last finalized
  std::vector<int> graph_locations_in_last_timestep(
      locations_in_last_timestep.size());
  std::transform(locations_in_last_timestep.begin(),
                 locations_in_last_timestep.end(),
                 graph_locations_in_last_timestep.begin(),
                 [this](const Index& index) { return to_graph_index(index); });

  std::vector<std::vector<int>> new_trajectories_in_graph;
  std::vector<int> node_index_to_shortest_path_incoming_edge;
  Solve(solver_type_, graph_, node_index_to_incoming_edges,
        node_index_to_outgoing_edges, graph_locations_in_last_timestep,
        new_trajectories_in_graph, node_index_to_shortest_path_incoming_edge,
        ignore_last_exit_cost);
  print("Found ", new_trajectories_in_graph.size(), " trajectories");

  // Convert locations in new trajectories to sequence indices and merge with
  // previously found trajectories. Also, copy into result.
  std::vector<std::vector<Index>> new_trajectories(
      new_trajectories_in_graph.size());
  for (std::size_t i = 0; i < new_trajectories.size(); ++i) {
    new_trajectories[i].resize(new_trajectories_in_graph[i].size());
    std::transform(new_trajectories_in_graph[i].begin(),
                   new_trajectories_in_graph[i].end(),
                   new_trajectories[i].begin(),
                   [this](int index) { return to_sequence_index(index); });
  }

  // The clipping time step is the first time step in the next optimization
  // window.
  // NOTE(nwojke): +1, because we optimize to the last finalized time step,
  // which is current_timestep_ - 1.
  const Index clipping_timestep =
      static_cast<Index>(std::max(0l, static_cast<long>(current_timestep_) -
                                          static_cast<long>(window_len_) + 1));

  auto location_to_timestep = [this](const Index& index) -> Index {
    // NOTE(nwojke): Mapping everything unknown to 0 is a bit unsafe, but
    // this way we can remove locations outside the optimization window from
    // the map.
    auto it = location_to_timestep_.find(index);
    return it != location_to_timestep_.end() ? it->second : 0;
  };

  MergeTrajectories(location_to_timestep, previous_clipping_timestep_,
                    clipping_timestep, new_trajectories, trajectories_, active_,
                    next_trajectory_index_);
  print("Number of trajectories in cache: ", trajectories_.size());

  // Collapse trajectories into the new trajectory head (index to new head is
  // stored in clipping_indices).
  std::vector<int> clipping_indices(new_trajectories.size());
  for (std::size_t i = 0; i < new_trajectories.size(); ++i) {
    for (int j = static_cast<int>(new_trajectories[i].size()) - 1; j >= 0;
         --j) {
      if (location_to_timestep(new_trajectories[i][j]) < clipping_timestep) {
        break;
      }
      clipping_indices[i] = j;
    }
  }
  CollapseTrajectories(new_trajectories_in_graph, clipping_indices,
                       node_index_to_incoming_edges,
                       node_index_to_shortest_path_incoming_edge, graph_);
  CleanupTrajectoryCache(location_to_timestep, clipping_timestep, trajectories_,
                         active_);

  // Pop time steps outside the next optimization window from graph.
  int num_locations_to_prune = 0;
  for (Index timestep = previous_clipping_timestep_;
       timestep < clipping_timestep; ++timestep) {
    num_locations_to_prune += timestep_to_locations_[timestep].size();
    for (const Index& location : timestep_to_locations_[timestep]) {
      location_to_timestep_.erase(location);
    }
    timestep_to_locations_.erase(timestep);
  }
  print("Pruning ", num_locations_to_prune, " locations.");
  PopFirstN(2 * num_locations_to_prune, graph_);  // Each location has 2 nodes.

  previous_clipping_timestep_ = clipping_timestep;
  num_pruned_locations_ += num_locations_to_prune;
  print("Total number of pruned locations ", num_pruned_locations_);
  print("Number of locations in graph: ", (graph_.num_nodes() - 2) / 2);
}

inline int BatchProcessing::to_graph_index(Index sequence_index) const {
  return sequence_index == ST ? sequence_index
                              : sequence_index - num_pruned_locations_;
}

inline BatchProcessing::Index BatchProcessing::to_sequence_index(
    int graph_index) const {
  return graph_index == static_cast<int>(ST)
             ? static_cast<Index>(graph_index)
             : static_cast<Index>(graph_index) + num_pruned_locations_;
}

}  // namespace mcf
