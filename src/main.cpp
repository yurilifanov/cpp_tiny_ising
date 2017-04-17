#include <iostream>             // std::cout
#include <cstdint>              // uint64_t
#include "ising.hpp"            // Ising 
#include "kinetic_ising.hpp"    // KineticIsing 
#include "timer.hpp"            // Timer
int main(int, char **argv) {
  KineticIsing<uint64_t> ising({3, 3}, std::stoll(argv[1]));
  Timer<double> timer;
  int64_t n_mv = std::stoll(argv[2]);
  ising.setkbt(std::stod(argv[3]));
  for(int64_t i = 0LL; i < n_mv; i++) {
    ising.mcmove();
  }
  std::cout << timer.toc() << "\n";
  std::cout << ising.get_num_ones() << " "
            << ising.get_t() << " "
            << ising.get_num_aligned() << "\n";
}
