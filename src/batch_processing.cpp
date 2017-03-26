// vim: expandtab:ts=2:sw=2
#include <mcf/batch_processing.hpp>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <stdexcept>

#include <mcf/internal/k_shortest_path.hpp>
#include <mcf/internal/util.hpp>

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
 * Check if two trajectories overlap and find the matching location.
 *
 * Locations in new_trajectory with index in [new_begin, new_end) hould be
 * added to the cached trajectory, such that the trajectory extends into the
 * new optimization window. This function computes new_begin and new_end for
 * this purpose.
 *
 * @param cached_trajectory A cached trajectory from previous solver runs.
 * @param new_trajectory A new trajectory found at the current solver run.
 * @param location_to_timestep A mapping from location to time step.
 * @param clipping_timestep The first time step in the new optimization window.
 *        The new trajectory head is the first node in new_trajectory with
 *        time step larger or equal to this value.
 * @param new_begin On success, points one element behind the matched location
 *        in new_trajectory.
 *        Takes on new_trajectory.size() if the trajectory is fully contained
 *        in the new optimization window.
 * @param new_end On success, points one element behind the new head in
 *        new_trajectory, i.e., the first element with time step larger or
 *        equal to clipping_timestep.
 *        Takes on new_trajectory.size() if the trajectory does not extend into
 *        the new optimization window.
 * @return True if cached_trajectory and new_trajectory share a common
 *        location.
 */
bool ExtendTrajectory(const std::vector<int>& cached_trajectory,
                      const std::vector<int>& new_trajectory,
                      const std::vector<int>& location_to_timestep,
                      const int clipping_timestep, std::size_t& new_begin,
                      std::size_t& new_end) {
  if (cached_trajectory.empty()) {
    // Empty trajectory has nothing to match against.
    return false;
  }

  const int previous_head = cached_trajectory.back();
  const int previous_head_timestep = location_to_timestep[previous_head];

  // Search for a matching location on the trajectory.
  std::size_t matched_location_index = 0;
  for (; matched_location_index < new_trajectory.size();
       ++matched_location_index) {
    const int location = new_trajectory[matched_location_index];
    if (location_to_timestep[location] > previous_head_timestep) {
      // No need to iterate further, remaining locations are beyond the
      // previous trajectory head in time.
      return false;
    }

    if (location == previous_head) {
      // Match found.
      break;
    }
  }

  // Search for new trajectory head in the next optimization window.
  std::size_t new_head_index = matched_location_index;
  for (; new_head_index < new_trajectory.size(); ++new_head_index) {
    const int location = new_trajectory[new_head_index];
    if (location_to_timestep[location] >= clipping_timestep) {
      // We found a location in the new optimization window.
      break;
    }
  }

  // new_begin should point one behind the matched location, new_end should
  // point one behind the new trajectory head.
  new_begin = matched_location_index + 1;  // <= new_trajectory.size()
  new_end = std::min(new_head_index + 1, new_trajectory.size());

  // If matching_location points to the last element in new_trajectory, then
  // new_begin points to new_trajectory.end() and we return true with
  // new_begin == new_end == new_trajectory.size().
  return new_begin <= new_trajectory.size();
}

/**
 * Merge cached trajectories and newly computed trajectories.
 *
 * @param location_to_timestep A mapping from location to time step.
 * @param clipping_timestep The first time step in the new optimization window.
 *        Cached trajectories should be extended, such that the new head has a
 *        time step larger or equal to this value.
 * @param new_trajectories The new trajectories, computed at the most recent
 *        time step.
 * @param cached_trajectories The cached trajectories, computed in previous
 *        time steps. Will be extended to clipping_timestep using nodes from
 *        matching trajectories in new_trajectories.
 * @param cached_labels Unique trajectory labels associated with cached
 *        trajectories.
 * @param full_trajectories Contains all trajectories that have been found up
 *        until the most recent time step in their full length.
 * @param full_labels Unique trajectory labels associated with
 *        full_trajectories.
 * @param label_to_uncached_trajectory_head Maps from trajectory index to
 *        trajectory head; used to keep uncached trajectories, i.e.,
 *        trajectories that are fully contained in the optimization window,
 *        consistently labeled throughout subsequent solver runs.
 * @param clipping_indices Contains the index of the new trajectory head in
 *        new_trajectories or -1 if the trajectory is not part of the cache.
 */
void MergeTrajectories(const std::vector<int> location_to_timestep,
                       const int clipping_timestep,
                       const std::vector<std::vector<int>>& new_trajectories,
                       std::vector<std::vector<int>>& cached_trajectories,
                       std::vector<int>& cached_labels,
                       std::vector<std::vector<int>>& full_trajectories,
                       std::vector<int>& full_labels,
                       std::vector<int>& clipping_indices,
                       std::map<int, int>& label_to_uncached_trajectory_head) {
  assert(cached_trajectories.size() == cached_labels.size());
  // TODO(nwojke): Refactor to shorter functions.

  // Reserve some space and initialize clipping locations to default value.
  cached_trajectories.reserve(cached_trajectories.size() +
                              new_trajectories.size());
  full_trajectories.reserve(cached_trajectories.size() +
                            new_trajectories.size());

  clipping_indices.resize(new_trajectories.size());
  std::fill_n(clipping_indices.begin(), new_trajectories.size(), -1);

  // First, try to extend cached trajectories.
  for (std::size_t cache_index = 0; cache_index < cached_trajectories.size();
       ++cache_index) {
    std::vector<int>& cached_trajectory = cached_trajectories[cache_index];
    const int cached_label = cached_labels[cache_index];

    std::size_t new_begin = std::numeric_limits<std::size_t>::max();
    std::size_t new_end = std::numeric_limits<std::size_t>::max();
    std::size_t trajectory_index = 0;
    for (; trajectory_index < new_trajectories.size(); ++trajectory_index) {
      if (clipping_indices[trajectory_index] >= 0) {
        // This trajectory has been previously matched, ignore.
        continue;
      }
      if (ExtendTrajectory(
              cached_trajectory, new_trajectories[trajectory_index],
              location_to_timestep, clipping_timestep, new_begin, new_end)) {
        // Found a matching trajectory.
        break;
      }
    }

    // Copy trajectory over to full_trajectories. If a match was found, extend
    // cached trajectory to new end and extend the full trajectory.
    full_trajectories.push_back(cached_trajectory);
    full_labels.push_back(cached_label);
    if (trajectory_index >= new_trajectories.size()) {
      // No match.
      continue;
    }

    // Here,
    // * new_begin points one element behind the matched location (it is
    //   the beginning of the newly seen part of the trajectory).
    // * new_end points one element behind the new trajectory head (it is the
    //   end of the newly seen part of the trajectory). If new_end points to
    //   the end of the matched trajectory, it does not extend into the new
    //   optimization window.
    const std::vector<int>& matched_trajectory =
        new_trajectories[trajectory_index];
    const int num_new_in_cache = new_end - new_begin;
    assert(num_new_in_cache >= 0 && "Invalid new trajectory segment.");

    cached_trajectory.reserve(cached_trajectory.size() + num_new_in_cache);
    cached_trajectory.insert(cached_trajectory.end(),
                             matched_trajectory.begin() + new_begin,
                             matched_trajectory.begin() + new_end);

    const int num_new_in_full = matched_trajectory.size() - new_begin;
    assert(num_new_in_full >= 0 && "Invalid new trajectory segment.");

    full_trajectories.back().reserve(matched_trajectory.size() +
                                     num_new_in_full);
    full_trajectories.back().insert(full_trajectories.back().end(),
                                    matched_trajectory.begin() + new_begin,
                                    matched_trajectory.end());

    // Store index of new trajectory head in clipping_indices.
    clipping_indices[trajectory_index] = new_end - 1;
  }

  auto AddTrajectory =
      [&new_trajectories, &location_to_timestep, &cached_trajectories,
       &cached_labels, &full_trajectories, &full_labels, &clipping_indices,
       &clipping_timestep](const int trajectory_index, const int label) {
        // First, append the entire trajectory to full_trajectories.
        const std::vector<int>& trajectory = new_trajectories[trajectory_index];
        full_trajectories.push_back(trajectory);
        full_labels.push_back(label);

        // Now, find the clipping location and add partial trajectory to the
        // cache.
        std::size_t clipping_index = 0;
        for (; clipping_index < trajectory.size(); ++clipping_index) {
          if (location_to_timestep[trajectory[clipping_index]] >=
              clipping_timestep) {
            break;
          }
        }

        if (clipping_index == 0) {
          // The full trajectory is within the next optimization window, do not
          // add it to cache.
          return false;
        }

        if (clipping_index == trajectory.size()) {
          // The trajectory does not extent into the next optimization window,
          // do
          // not add it to cache.
          return false;
        }

        // Store clipping location.
        clipping_indices[trajectory_index] = clipping_index;

        const int end_index = std::min(clipping_index + 1, trajectory.size());
        cached_trajectories.resize(cached_trajectories.size() + 1);
        cached_trajectories.back().reserve(end_index);
        cached_trajectories.back().insert(cached_trajectories.back().end(),
                                          trajectory.begin(),
                                          trajectory.begin() + end_index);
        cached_labels.emplace_back(label);
        return true;
      };

  // Now try adding previously seen, but uncached trajectories with consistent
  // label.
  // TODO(nwojke): Clean this up, the array below should not be necessary.
  std::vector<bool> trajectory_is_new(clipping_indices.size());
  for (std::size_t index = 0; index < clipping_indices.size(); ++index) {
    trajectory_is_new[index] = clipping_indices[index] < 0;
  }

  for (auto index_and_head = label_to_uncached_trajectory_head.begin();
       index_and_head != label_to_uncached_trajectory_head.end();) {
    // Search for a matching trajectory that is not yet processed.
    std::size_t trajectory_index = 0;
    for (; trajectory_index < new_trajectories.size(); ++trajectory_index) {
      assert(!new_trajectories[trajectory_index].empty());
      if (index_and_head->second ==
              new_trajectories[trajectory_index].front() &&
          clipping_indices[trajectory_index] < 0) {
        break;
      }
    }

    bool trajectory_is_cached = false;
    if (trajectory_index < new_trajectories.size()) {
      trajectory_is_cached =
          AddTrajectory(trajectory_index, index_and_head->first);
      trajectory_is_new[trajectory_index] = false;
    }

    if (trajectory_is_cached ||
        location_to_timestep[index_and_head->second] < clipping_timestep) {
      index_and_head = label_to_uncached_trajectory_head.erase(index_and_head);
    } else {
      ++index_and_head;
    }
  }

  // Now, carry over remaining trajectories to cache and full_trajectories.
  for (std::size_t trajectory_index = 0;
       trajectory_index < new_trajectories.size(); ++trajectory_index) {
    if (!trajectory_is_new[trajectory_index]) {
      continue;
    }

    const int label = full_trajectories.size();
    if (!AddTrajectory(trajectory_index, label)) {
      const int index = full_trajectories.size() - 1;
      const int location = new_trajectories[trajectory_index].front();
      label_to_uncached_trajectory_head.insert(std::make_pair(index, location));
    }
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
 * Pop first n non-source/sink nodes from the graph.
 *
 * @param n Number of nodes to remove.
 * @param graph The graph.
 */
void PopFirstN(const int n, Graph& graph) {
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

}  // unnamed namespace

BatchProcessing::BatchProcessing(const int window_len,
                                 ShortestPathSolverType solver_type)
    : solver_type_(solver_type),
      window_len_(window_len),
      current_timestep_(0),
      previous_clipping_timestep_(0),
      num_pruned_locations_(0),
      location_to_timestep_(1, -1),
      timestep_to_locations_(1) {}

void BatchProcessing::Reserve(int num_edges) { graph_.Reserve(num_edges); }

int BatchProcessing::Add(const double cost) {
  const int graph_index = graph_.Add(cost);
  const int location = to_sequence_index(graph_index);

  location_to_timestep_.resize(location_to_timestep_.size() + 1);
  assert(static_cast<int>(timestep_to_locations_.size()) ==
         current_timestep_ + 1);
  assert(static_cast<int>(location_to_timestep_.size()) == location + 1);

  timestep_to_locations_[current_timestep_].push_back(location);
  location_to_timestep_[location] = current_timestep_;
  return location;
}

void BatchProcessing::Link(int src, int dst, double cost) {
  if ((src != ST && src < num_pruned_locations_) ||
      (dst != ST && dst < num_pruned_locations_)) {
    std::stringstream msg;
    msg << "Cannot link a pruned location. Source: " << src
        << " target: " << dst
        << " first non-pruned location: " << num_pruned_locations_;
    throw std::invalid_argument(msg.str().c_str());
  }

  src = to_graph_index(src);
  dst = to_graph_index(dst);
  assert(src >= 0 && src <= graph_.num_nodes());
  assert(dst >= 0 && dst <= graph_.num_nodes());

  graph_.Link(src, dst, cost);
}

void BatchProcessing::FinalizeTimeStep() {
  timestep_to_locations_.resize(timestep_to_locations_.size() + 1);
  ++current_timestep_;
}

void BatchProcessing::RunSearch(std::vector<std::vector<int>>& trajectories,
                                const bool ignore_last_exit_cost) {
  const int num_timesteps_in_graph =
      current_timestep_ - previous_clipping_timestep_;
  print("Batch processor called at time step ", current_timestep_);
  print("Number of time steps in graph: ", num_timesteps_in_graph);
  if (num_timesteps_in_graph <= 0) {
    // No new data in the graph, no need to solve.
    trajectories = trajectories_;
    return;
  }

  // Run solver on full graph up to current_timestep_.
  std::vector<std::vector<int>> node_index_to_incoming_edges;
  std::vector<std::vector<int>> node_index_to_outgoing_edges;
  internal::BuildEdgeMap(graph_, node_index_to_incoming_edges,
                         node_index_to_outgoing_edges);

  std::vector<int> locations_in_last_timestep =
      timestep_to_locations_[timestep_to_locations_.size() - 2];
  for (int& location : locations_in_last_timestep) {
    location = to_graph_index(location);
  }

  std::vector<std::vector<int>> new_trajectories_in_graph;
  std::vector<int> node_index_to_shortest_path_incoming_edge;
  Solve(solver_type_, graph_, node_index_to_incoming_edges,
        node_index_to_outgoing_edges, locations_in_last_timestep,
        new_trajectories_in_graph, node_index_to_shortest_path_incoming_edge,
        ignore_last_exit_cost);
  print("Found ", new_trajectories_in_graph.size(), " trajectories");

  // Convert locations in new trajectories to sequence indices and merge with
  // previously found trajectories. Also, copy into result.
  std::vector<std::vector<int>> new_trajectories_in_sequence =
      new_trajectories_in_graph;
  for (std::vector<int>& trajectory : new_trajectories_in_sequence) {
    for (int& location : trajectory) {
      location = to_sequence_index(location);
    }
  }

  // The clipping time step is the first time step in the next optimization
  // window. It is also the time step at which we want to place the new
  // cached trajectory head. Everything until clipping_timestep - 1 will be
  // forgotten.
  // NOTE(nwojke): +1, because we optimize to the last finalized time step,
  // which is current_timestep_ - 1.
  const int clipping_timestep =
      std::max(0, current_timestep_ - window_len_ + 1);
  std::vector<int> clipping_indices;
  std::vector<std::vector<int>> full_trajectories;
  std::vector<int> full_labels;
  MergeTrajectories(location_to_timestep_, clipping_timestep,
                    new_trajectories_in_sequence, trajectories_,
                    trajectory_labels_, full_trajectories, full_labels,
                    clipping_indices, label_to_uncached_trajectory_head_);

  trajectories.resize(full_trajectories.size());
  for (std::size_t i = 0; i < full_trajectories.size(); ++i) {
    assert(trajectories[full_labels[i]].empty());
    trajectories[full_labels[i]] = full_trajectories[i];
  }

  print("Total number of trajectories over the entire sequence: ",
        trajectories.size());

  // Collapse trajectories into the new trajectory head (index to new head is
  // stored in clipping_indices).
  CollapseTrajectories(new_trajectories_in_graph, clipping_indices,
                       node_index_to_incoming_edges,
                       node_index_to_shortest_path_incoming_edge, graph_);

  // Pop time steps outside the next optimization window from graph.
  int num_locations_to_prune = 0;
  for (int timestep = previous_clipping_timestep_; timestep < clipping_timestep;
       ++timestep) {
    num_locations_to_prune += timestep_to_locations_[timestep].size();
  }
  print("Pruning ", num_locations_to_prune, " locations.");
  PopFirstN(2 * num_locations_to_prune, graph_);  // Each location has 2 nodes.

  previous_clipping_timestep_ = clipping_timestep;
  num_pruned_locations_ += num_locations_to_prune;
  print("Total number of pruned locations ", num_pruned_locations_);
  print("Number of locations in graph: ", (graph_.num_nodes() - 2) / 2);
}

inline int BatchProcessing::to_graph_index(const int sequence_index) const {
  return sequence_index == ST ? sequence_index
                              : sequence_index - num_pruned_locations_;
}

inline int BatchProcessing::to_sequence_index(const int graph_index) const {
  return graph_index == ST ? graph_index : graph_index + num_pruned_locations_;
}

}  // namespace mcf
