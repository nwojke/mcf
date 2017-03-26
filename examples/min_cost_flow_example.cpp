// vim: expandtab:ts=2:sw=2
#include <iostream>

#include <mcf/graph.hpp>
#include <mcf/k_shortest_path_solver.hpp>

#ifdef MCF_USE_Lemon
#include <mcf/lemon_solver.hpp>
#endif

#ifdef MCF_USE_Clp
#include <mcf/clp_solver.hpp>
#endif

int main(int argc, char** argv) {
  mcf::Graph g;

  int n[4];
  n[0] = g.Add(-1.0);
  n[1] = g.Add(1.0);
  n[2] = g.Add(1.0);
  n[3] = g.Add(-1.0);

  for (int i = 0; i < 4; ++i) {
    std::cout << "node[" << i << "] has ID " << n[i] << std::endl;
  }

  g.Link(g.ST, n[0], 0.0);
  g.Link(g.ST, n[1], 1.0);
  g.Link(n[2], g.ST, 1.0);
  g.Link(n[3], g.ST, 1.0);
  g.Link(n[0], n[2], 1.0);
  g.Link(n[1], n[3], 1.0);

  std::vector<std::vector<int>> trajectories;

#ifdef MCF_USE_Lemon
  std::cout << "Solving with Lemon" << std::endl;
  mcf::LemonSolver(g).RunSearch(1, 2, trajectories);

  std::cout << "num trajectories: " << trajectories.size() << std::endl;
  for (const std::vector<int>& trajectory : trajectories) {
    std::cout << "trajectory -----" << std::endl;
    for (int loc : trajectory) {
      std::cout << "\t node: " << loc << std::endl;
    }
  }
#endif

#ifdef MCF_USE_Clp
  std::cout << "Solving with CLP" << std::endl;
  mcf::ClpSolver(g).RunSearch(1, 2, trajectories);

  std::cout << "num trajectories: " << trajectories.size() << std::endl;
  for (const std::vector<int>& trajectory : trajectories) {
    std::cout << "trajectory -----" << std::endl;
    for (int loc : trajectory) {
      std::cout << "\t node: " << loc << std::endl;
    }
  }
#endif

  std::cout << "Solving with Successive Shortest Paths algorithm" << std::endl;
  mcf::ShortestPathSolver(g).RunSearch(1, 2, trajectories);
  for (const std::vector<int>& trajectory : trajectories) {
    std::cout << "trajectory -----" << std::endl;
    for (int loc : trajectory) {
      std::cout << "\t node: " << loc << std::endl;
    }
  }

  return 0;
}
