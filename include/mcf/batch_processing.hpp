// vim: expandtab:ts=2:sw=2
#ifndef MCF_BATCH_PROCESSING_HPP
#define MCF_BATCH_PROCESSING_HPP

#include <mcf/graph.hpp>
#include <mcf/k_shortest_path_solver.hpp>

#include <map>
#include <unordered_map>

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
 * * To divide the sequence into fixed batches of 100 frames each, set the
 *   constructor's window_len parameter to 100 and call RunSearch() every
 *   100 frames.
 * * To estimate trajectories in an on-line fashion, call RunSearch() at each
 *   time step (after calling FinalizeTimeStep()). Set the constructor's
 *   window_len parameter to the history of time steps that should be optimized
 *   (similar to N-back-scan pruning in multiple hypothesis tracking).
 * * As alternative, call ComputeTrajectories() and RemoveInactiveTracks() at
 *   each time step (after a call to FinalizeTimeStep()) to only keep track of
 *   objects that are within the optimization window_len. This prevents the
 *   number of trajectories in the cache from growing unbounded.
 */
class BatchProcessing {
 public:
  using Index = uint64_t;
  using Trajectory = std::vector<Index>;
  using TrajectoryMap = std::map<Index, Trajectory>;

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
  static const Index ST;

  //! Reserve memory for a given number of edges.
  void Reserve(int num_edges);

  //! Add a location with given cost. See Graph::Add() for more information.
  Index Add(double cost);

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
  void Link(Index src, Index dst, double cost);

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
   *        sequence of location handles returned by Add(). Trajaectories are
   *        returned in a globally consistent order, such that each trajectory
   *        keeps its index over subsequent calls.
   * @param ignore_last_exit_cost If true, the exit cost for locations in the
   *        last time step is set to 0 prior to calling the solver. This only
   *        affects the current trajectory search, future calls to RunSearch()
   *        work on the original graph structure.
   */
  void RunSearch(std::vector<Trajectory>& trajectories,
                 bool ignore_last_exit_cost = true);

  /**
   * Find multi-object trajectory.
   *
   * @param ignore_last_exit_cost If true, the exit cost for locations in the
   *        last time step is set to 0 prior to calling the solver. This only
   *        affects the current trajectory search, future calls to RunSearch()
   *        work on the original graph structure.
   * @param removed_indices If not nullptr, filled with the list of trajectory
   *        indices that have been removed from the cache because they are
   *        not part of the final solution.
   * @return Maps from trajectory index/identifier to list of location handles
   *         returned by Add().
   *
   */
  TrajectoryMap ComputeTrajectories(
      bool ignore_last_exit_cost = true,
      std::vector<Index>* removed_indices = nullptr);

  /**
   * Remove trajectories from the cache that are not contained in the
   * optimization window. This affects the list of trajectories returned
   * by ComputeTrajectories() as well as RunSearch().
   *
   * @return A list of trajectory indices that are part of the final solution,
   *         but have been removed from the cache because they are outside of
   *         the optimization window.
   */
  std::vector<Index> RemoveInactiveTracks();

  /**
   * Returns the minimum active location handle. A location is active if it is
   * inside the current optimization window or part of a cached trajectory.
   * A location smaller than this value is guaranteed to not be part of a
   * trajectory in future calls to ComputeTrajectories().
   */
  Index ComputeMinActiveLocation() const;

 private:
  void Update(bool ignore_last_exit_cost, std::vector<Index>& removed_indices);

  int to_graph_index(Index sequence_index) const;

  Index to_sequence_index(int graph_index) const;

  ShortestPathSolverType solver_type_;
  Index window_len_;

  Index current_timestep_;
  Index previous_clipping_timestep_;
  Index num_pruned_locations_;

  std::unordered_map<Index, Index> location_to_timestep_;
  std::unordered_map<Index, std::vector<Index>> timestep_to_locations_;

  Graph graph_;

  TrajectoryMap trajectories_;
  std::unordered_map<Index, bool> active_;
  Index next_trajectory_index_;
};

}  // namespace mcf

#endif
