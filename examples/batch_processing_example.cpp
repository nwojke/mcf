// vim: expandtab:ts=2:sw=2
#include <iostream>

#include <mcf/batch_processing.hpp>
#include <mcf/graph.hpp>
#include <mcf/logging.hpp>

int main(int argc, char** argv) {
  mcf::printer().set_verbose(true);

  // We build two trajectories with nodes fully connected between successive
  // time steps.
  static constexpr int kWindowLen = 3;
  mcf::BatchProcessing batch_processing(kWindowLen);
  std::vector<int> heads = {mcf::BatchProcessing::ST, mcf::BatchProcessing::ST};

  static constexpr int kTrajectoryLength = 10;
  for (int timestep = 0; timestep < kTrajectoryLength; ++timestep) {
    std::vector<int> new_heads = {batch_processing.Add(-1.0),
                                  batch_processing.Add(-1.0)};

    for (int i = 0; i < 2; ++i) {
      batch_processing.Link(heads[i], new_heads[i], 0.0);
      batch_processing.Link(new_heads[i], mcf::BatchProcessing::ST, 0.5);
    }

    if (timestep > 0) {
      for (int i = 0; i < 2; ++i) {
        batch_processing.Link(mcf::BatchProcessing::ST, new_heads[i], 0.5);
      }
      batch_processing.Link(heads[0], new_heads[1], 1.0);
      batch_processing.Link(heads[1], new_heads[0], 1.0);
    }

    batch_processing.FinalizeTimeStep();
    heads = std::move(new_heads);

    // Solve at (almost) every time step, to check that this works as long as
    // the time since the last solver run is less than kWindowLen.
    if (timestep % 4 != 0) {
      continue;
    }

    std::vector<std::vector<int>> trajectories;
    batch_processing.RunSearch(trajectories);

    std::cout << "Time step " << timestep << std::endl;
    std::cout << "Number of trajectories: " << trajectories.size() << std::endl;
    for (std::size_t i = 0; i < trajectories.size(); ++i) {
      std::cout << "  Trajectory " << i << ":";
      for (const int& location : trajectories[i]) {
        std::cout << " " << location;
      }
      std::cout << std::endl;
    }
  }

  return 0;
}
