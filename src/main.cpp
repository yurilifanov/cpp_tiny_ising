#include <iostream>             // std::cout
#include <cstdint>              // uint64_t
#include <exception>            // std::exception
#include <vector>               // std::vector 
#include <fstream>              // std::ofstream
#include <boost/progress.hpp>   // boost::progress_display
#include "kinetic_ising.hpp"    // KineticIsing 
#include "timer.hpp"            // Timer
template <typename Cfg>
void test(Cfg cfg_0, int64_t num_mv, double kbt) {
  KineticIsing<Cfg> ising({8, 8}, cfg_0);
  Timer<double> timer;
  int64_t n_mv = num_mv;
  ising.setkbt(kbt);
  std::cout << "Testing run time...\n";
  for(int64_t i = 0LL; i < n_mv; i++) {
    ising.mcmove();
  }
  std::cout << timer.toc() << "\n"
            << ising.get_num_ones() << " "
            << ising.get_t() << " "
            << ising.get_num_aligned() << "\n";
}
void fpt(int64_t n, double kbt, double h) {
  using Cfg = uint64_t;
  KineticIsing<Cfg> ising({3, 2});
  ising.setkbt(kbt);
  ising.seth(h);
  std::vector<double> fpts(n, 0.);
  std::cout << "Sampling first passage times...\n";
  boost::progress_display pd(n);
  for(auto & val : fpts) {
    ising.setcfg(Cfg(0));
    while(ising.get_num_ones() < ising.get_N()) {
      ising.mcmove();
    }
    val = ising.get_t();
    ++pd;
  }
  std::ofstream fl("fpts.bin", std::ios::binary);
  fl.write((const char *) fpts.data(), fpts.size() * sizeof(double));
  fl.close();
}
void tm_data() {
  Ising<uint64_t> ising({3, 2});
  ising.save_transition_matrix_data();
}
int main(int, char **argv) {
  using namespace std;
  try {
    //test<uint64_t>(stoull(argv[1]), stoll(argv[2]), stod(argv[3]));
    fpt(stoll(argv[1]), stod(argv[2]), stod(argv[3]));
    tm_data();
  } catch(exception & e) {
    cout << e.what() << "\n";
  }
}
