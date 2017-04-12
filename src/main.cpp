#include <iostream>     // std::cout
#include <cstdint>      // int64_t
#include <limits>       // std::numeric_limits
#include <bitset>       // std::bitset
#include <vector>       // std::vector
#include <functional>   // std::multiplies
#include <numeric>      // std::partial_sum
#include <algorithm>    // std::rotate
#include <string>       // std::stoll

template<typename T>
inline int popcount(T val) noexcept {
  using namespace std;
  return bitset<numeric_limits<T>::digits>(val).count();
}

template<typename T, int n = std::numeric_limits<T>::digits>
inline T rotl(T x, int k) { 
  return (x << k) | (x >> (n - k)); 
}

class Ising {
    using v_int64 = std::vector<int64_t>;
    const v_int64 len;
    const v_int64 cumprod_len;
    const int64_t N;
    const int64_t d;
    const int64_t d_times_2;
    const v_int64 masks;
    double kbt;
    double beta;
    double h;
    int64_t cfg;
    int num_aligned_pairs;
    int num_ones;
    v_int64 get_cumprod_len() const {
      auto ans = len;
      ans.emplace(ans.begin(), 1LL);
      std::partial_sum(
                      ans.begin(), 
                      ans.end(), 
                      ans.begin(), 
                      std::multiplies<int64_t>());
      return ans;
    }
    v_int64 get_masks() const {
      v_int64 ans(N, 0LL);
      for(int64_t i = 0LL; i < N; i++) {
        //std::cout << i << ": ";
        int64_t mask = 0LL;
        for(int64_t j = 0LL; j < d_times_2; j++) {
          //std::cout << get_neighbour_index(i, j) << ", ";
          mask |= 1LL << get_neighbour_index(i, j);
        }
        //std::cout << std::bitset<9>(mask) << "\n";
        ans[i] = mask;
      }
      return ans;
    }
    int64_t get_neighbour_index(int64_t i, int64_t j) const {
      int64_t pm[] = {-1LL, 1LL};
      int64_t i_diff = pm[j % 2LL];
      int64_t crd_id = j / 2LL;
      int64_t crd_old = (i / cumprod_len[crd_id]) % len[crd_id];
      int64_t crd_new = (len[crd_id] + crd_old + i_diff) % len[crd_id];
      return i + (crd_new - crd_old) * cumprod_len[crd_id];
    }
    int64_t compute_num_aligned_pairs() const {
      int64_t num = 0LL;
      for(int64_t i = 0; i < N; i++) {
        num += get_num_aligned_neighbours_for_site(cfg, i); 
      }
      return num / 2;
    }
    int64_t get_num_aligned_neighbours_for_site(int64_t mcfg, int64_t i) const {
      int64_t n1s = popcount(mcfg & masks[i]);
      int64_t lkp[] = {d_times_2 - n1s, n1s};
      //std::cout << lkp[(mcfg >> i) & 1LL] << "\n";
      return lkp[(mcfg >> i) & 1LL];
    }
  public:
    Ising(v_int64 len_in)
      : len(len_in)
      , cumprod_len(get_cumprod_len())
      , N(cumprod_len.back())
      , d(len_in.size())
      , d_times_2(2LL * d)
      , masks(get_masks())
      {
        if(N > 64) {
          std::cout << "Number of sites is too great for int64 representation."
                    << "\nYour should really throw an exception here.\n";
        }
        setkbt(1.);
        seth(0.);
        setcfg(0LL);
      }
    void setkbt(double val) { kbt = val; beta = 1. / val; }  
    void seth(double val) { h = val; }
    void setcfg(int64_t val) { 
      if(val >= (1LL << N)) {
        cfg = (1LL << N) - 1LL;
      } else if(val < 0LL) {
        cfg = 0LL;
      } else {
        cfg = val;
      }
      num_aligned_pairs = compute_num_aligned_pairs();
      num_ones = popcount(cfg);
      //std::cout << "Num aligned: " << num_aligned_pairs << "\n";
      //std::cout << "Num ones: " << num_ones << "\n";
    }
    void mcmove();
};

inline void Ising::mcmove() {
  int64_t i = 0LL; //N * rng();
  int64_t cfg_new = cfg ^ (1LL << i);
  int64_t mask = masks[i];
  // ...
  //int64_t lkp[] = {cfg, cfg_new};
  //cfg = lkp[rng() < acc_prob_lkp[aligned_delta][ones_delta]];
}

int main(int argc, char **argv) {
  Ising ising({3, 3});
  ising.setcfg(std::stoll(argv[1]));
}
