#include <iostream>     // std::cout
#include <cstdint>      // uint64_t
#include "ising.hpp"    // Ising 
#include "timer.hpp"    // Timer
int main(int, char **argv) {
  Ising<uint64_t> ising({8, 8}, std::stoll(argv[1]));
  Timer<double> timer;
  int64_t n_mv = std::stoll(argv[2]);
  ising.setkbt(std::stod(argv[3]));
  for(int64_t i = 0LL; i < n_mv; i++) {
    ising.mcmove();
  }
  std::cout << timer.toc() << "\n";
  std::cout << ising.get_num_ones() << " "
            << ising.get_num_aligned() << "\n";
}
