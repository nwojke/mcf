// vim: expandtab:ts=2:sw=2
#ifndef MCF_BATCH_PROCESSING_HPP
#define MCF_BATCH_PROCESSING_HPP

#include <mcf/graph.hpp>
#include <mcf/k_shortest_path_solver.hpp>

#include <map>

namespace mcf {

/***
 * This class implements a multi-object trajectory solver that processes the
 * sequence in batches. For this purpose, the class mimics the interface of
 * the Graph class to build a flow network incrementally. Similar to the solver
 * interface, trajectories can be optimized using RunSearch().
 *
 * Internally, the class caches the solution obtained from RunSearch() to build
 * the set of trajectories incrementally, without solving over the entire
 * sequence at each call. On construction, a minimum optimization window length
 * can be specified.
 *
 * Exemplary use-cases
 * -------------------
 *
 * * To estimate trajectories in an on-line fashion, call RunSearch() at each
 *   time step (after calling FinalizeTimeStep()). Set the constructor's
 *   window_len parameter to the history of time steps that should be optimized
 *   (similar to N-back-scan pruning in multiple hypothesis tracking).
 * * To divide the sequence into fixed batches of 100 frames each, set the
 *   constructor's window_len parameter to 100 and call RunSearch() every
 *   100 frames.
 */
class BatchProcessing {
 public:
  /**
   * Constructor.
   *
   * @param window_len The minimum optimization window length. After each
   *        call to RunSearch(), nodes are pruned from the internal graph
   *        structure when they are older than this parameter.
   * @param solver_type The type of shortest path algorithm to use during
   *        successive shortest path search.
   */
  BatchProcessing(int window_len,
                  ShortestPathSolverType solver_type =
                      ShortestPathSolverType::kDijkstraLazyDeletion);

  //! Convenience typedef. See Graph::ST for more information.
  static const int ST = Graph::ST;

  //! Reserve memory for a given number of edges.
  void Reserve(int num_edges);

  //! Add a location with given cost. See Graph::Add() for more information.
  int Add(double cost);

  /**
   * Link two locations with given cost. See Graph::Link() for more
   * information.
   *
   * @param src location handle of the source location.
   * @param dst location handle of the target location.
   * @param cost Transition edge cost (see class description).
   *
   * @throws std::invalid_argument If src or dst are outside the current
   *         optimization window, i.e., have been added before the last call to
   *         RunSearch()).
   */
  void Link(int src, int dst, double cost);

  /**
   * Notify the batch processor that the current time step has been finalized.
   * Must be called exactly once after every time step. No location that has
   * been added before this call may appear as a target location in subsequent
   * calls to Link(). The only exception is a link to/from the source/sink
   * node.
   */
  void FinalizeTimeStep();

  /**
   * Find multi-object trajectory.
   *
   * @param trajectories Computed trajectories, where each trajectory is a
   *        sequence of location handles returned by Add().
   * @param ignore_last_exit_cost If true, the exit cost for locations in the
   *        last time step is set to 0 prior to calling the solver. This only
   *        affects the current trajectory search, future calls to RunSearch()
   *        work on the original graph structure.
   */
  void RunSearch(std::vector<std::vector<int>>& trajectories,
                 bool ignore_last_exit_cost = true);

 private:
  int to_graph_index(int sequence_index) const;

  int to_sequence_index(int graph_index) const;

  ShortestPathSolverType solver_type_;
  int window_len_;
  int current_timestep_;
  int previous_clipping_timestep_;
  int num_pruned_locations_;

  std::vector<int> location_to_timestep_;
  std::vector<std::vector<int>> timestep_to_locations_;
  std::vector<std::vector<int>> trajectories_;
  std::vector<int> trajectory_labels_;

  Graph graph_;

  std::map<int, int> label_to_noncached_trajectory_head_;
  int next_label_;
};

}  // namespace mcf

#endif
