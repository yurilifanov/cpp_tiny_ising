#include <iostream>     // std::cout
#include <cstdint>      // int64_t
#include <limits>       // std::numeric_limits
#include <bitset>       // std::bitset
#include <vector>       // std::vector
#include <functional>   // std::multiplies
#include <numeric>      // std::partial_sum
#include <algorithm>    // std::rotate
#include <string>       // std::stoll
#include "rng.hpp"      // RNG_ull, Rng_d

/*
    TO DO:
      - Make Ising a template <typename Cfg>
      - Overload appropriate operators for arbitrary cfg type
      - e.g. __int128, __m256i, __m128i etc 
*/

template <typename T>
inline constexpr int num_digits() {
  return std::numeric_limits<T>::digits;
}

template <typename T>
inline int popcount(T val) noexcept {
  return std::bitset<num_digits<T>()>(val).count();
}

template <typename Cfg>
class Ising {
    using v_cfg = std::vector<Cfg>;
    using v_int64 = std::vector<int64_t>;

    const v_int64 len;
    const v_int64 cumprod_len;
    const int64_t N;
    const int64_t d;
    const int64_t d_times_2;
    const v_cfg   masks;
    double        kbt;
    double        beta;
    double        h;
    Cfg           cfg;
    int           num_aligned_pairs;
    int           num_ones;
    RNG_d<false>  rng;

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
    v_cfg get_masks() const {
      v_cfg ans(N, Cfg(0));
      for(int64_t i = 0LL; i < N; i++) {
        //std::cout << i << ": ";
        Cfg mask = Cfg(0);
        for(int64_t j = 0LL; j < d_times_2; j++) {
          //std::cout << get_neighbour_index(i, j) << ", ";
          mask |= Cfg(1) << get_neighbour_index(i, j);
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
    int64_t get_num_aligned_neighbours_for_site(Cfg mcfg, int64_t i) const {
      int64_t n1s = popcount(mcfg & masks[i]);
      int64_t lkp[] = {d_times_2 - n1s, n1s};
      return lkp[(mcfg >> i) & Cfg(1)];
    }
  public:
    Ising(v_int64 len_in, Cfg i_cfg = Cfg(0))
      : len(len_in)
      , cumprod_len(get_cumprod_len())
      , N(cumprod_len.back())
      , d(len_in.size())
      , d_times_2(2LL * d)
      , masks(get_masks())
      {
        if(N > num_digits<Cfg>()) {
          std::cout << "Number of sites is too great for representation."
                    << "\nYour should really throw an exception here.\n";
        }
        setkbt(1.);
        seth(0.);
        setcfg(i_cfg);
      }
    void setkbt(double val) { kbt = val; beta = 1. / val; }  
    void seth(double val) { h = val; }
    void setcfg(Cfg val) { 
      if(val >= (Cfg(1) << N)) {
        cfg = (Cfg(1) << N) - Cfg(1);
      } else if(val < Cfg(0)) {
        cfg = Cfg(0);
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

template <typename Cfg>
inline void Ising<Cfg>::mcmove() {
  int64_t i = N * rng();
  Cfg one(1);
  Cfg mask = masks[i];
  Cfg cfg_new = cfg ^ (one << i);
  int d_ones_id_lkp[] = {0L, 1L};
  int ones_delta_id = d_ones_id_lkp[(cfg_new >> i) & one];
  Cfg lkp[] = {cfg, cfg_new};
  //cfg = lkp[rng() < acc_prob_lkp[aligned_delta][ones_delta_id]];
}

int main(int, char **argv) {
  Ising<int64_t> ising({3, 3}, std::stoll(argv[1]));
  //std::cout << (int)__int128(3) << "\n";
}
