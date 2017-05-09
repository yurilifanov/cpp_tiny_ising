#ifndef _KINETIC_ISING_
#define _KINETIC_ISING_
#include <cstdint>          // int64_t
#include <vector>           // std::vector
#include <numeric>          // std::partial_sum
#include <functional>       // std::plus
#include <cmath>            // log
#include "ising.hpp"        // Ising
#include "rng.hpp"          // RNG_exp
/*
    This is a slightly optimised implementation of the n-fold method of 
    A. B. Bortz and M. H. Kalos and J. L. Lebowitz: 
    A new algorithm for Monte Carlo simulation of Ising spin systems (1975).
    It is a rejection free kMC scheme which works by first computing the set of
    the possible events followed by a random selection.
    Needless to say, it is quite slow per call to mcmove().
*/
template <typename Cfg>
class KineticIsing_extra_slow : public Ising<Cfg> {
    using v_int64 = std::vector<int64_t>;
    using v_double = std::vector<double>;
    RNG_exp rng_exp;
    double t;
    v_double flip_probs;
    void get_flip_delta(
      int64_t i, 
      int64_t & d_aligned, 
      int64_t & d_ones,
      bool & is_one
    ) const {
      static const Cfg one = Cfg(1);
      Cfg cfg_new = Ising<Cfg>::cfg ^ (one << i);
      is_one = (cfg_new >> i) & one;
      d_ones = is_one ? 1LL:-1LL;
      int64_t n_nbr_ones = popcount(cfg_new & Ising<Cfg>::masks[i]);
      d_aligned = d_ones * (2LL * n_nbr_ones - Ising<Cfg>::d_times_2);
    }
    double get_flip_weight(int64_t i) const {
      int64_t d_ones = 0LL;
      int64_t d_aligned = 0LL;
      bool is_one = false;
      get_flip_delta(i, d_aligned, d_ones, is_one);
      int64_t d_aligned_id = Ising<Cfg>::d_times_2 + d_aligned;
      return Ising<Cfg>::acc_weight_lkp[d_aligned_id][is_one];
    }
    double get_flip_prob(int64_t i) const {
      double fw = get_flip_weight(i);
      return fw > 1. ? 1.:fw;
    }
  public:
    KineticIsing_extra_slow(v_int64 len_in, Cfg i_cfg = Cfg(0))
      : Ising<Cfg>(len_in, i_cfg)
      , rng_exp(Ising<Cfg>::get_N())
      , t(0.)
      , flip_probs(Ising<Cfg>::N, 0.)
    {}
    //here self-loops ARE NOT counted as jumps
    void mcmove() { 
      int64_t i = 0LL;
      for(auto & val : flip_probs) {
        val = get_flip_prob(i++);
      }
      std::partial_sum(
        flip_probs.begin(), 
        flip_probs.end(), 
        flip_probs.begin(), 
        std::plus<double>()
      );
      double rnd = flip_probs.back() * Ising<Cfg>::rng();
      i = 0LL;
      for(auto val : flip_probs) {
        if(rnd < val) {
          break;
        }
        ++i;
      }
      static const Cfg one = Cfg(1);
      int64_t d_ones = 0LL;
      int64_t d_aligned = 0LL;
      bool is_one = false;
      get_flip_delta(i, d_aligned, d_ones, is_one);
      Ising<Cfg>::cfg ^= one << i;
      Ising<Cfg>::num_aligned_pairs += d_aligned;
      Ising<Cfg>::num_ones += d_ones;
      t += rng_exp();
    }
    void setcfg(Cfg i_cfg) { Ising<Cfg>::setcfg(i_cfg); t = 0.; }
    double get_t() { return t; }
};
/*
  This implementation works too, but is not rejection free and therefore a bit
  faster than the above (per call to mcmove). It is, however, significantly 
  slowed down by the call to log() in rng_exp
*/
template <typename Cfg>
class KineticIsing_slow : public Ising<Cfg> {
    using v_int64 = std::vector<int64_t>;
    RNG_exp rng_exp;
    double t;
  public:
    KineticIsing_slow(v_int64 len_in, Cfg i_cfg = Cfg(0))
      : Ising<Cfg>(len_in, i_cfg)
      , rng_exp(Ising<Cfg>::get_N())
      , t(0.)
    {}
    //here self-loops ARE NOT counted as jumps
    void mcmove() { 
      if(Ising<Cfg>::mcmove()) {
        t += rng_exp(); 
      } 
    }
    void setcfg(Cfg i_cfg) { Ising<Cfg>::setcfg(i_cfg); t = 0.; }
    double get_t() const { return t; }
};
/*
  This implementation is like the above, but exploits the sum property of 
  independent exponential variates. Calling get_t is more expensive, but 
  mcmove is almost as cheap as in Ising.
*/
template <typename Cfg>
class KineticIsing : public Ising<Cfg> {
    using v_int64 = std::vector<int64_t>;
    RNG_gamma rng_gamma;
    uint64_t num_jumps;
  public:
    KineticIsing(v_int64 len_in, Cfg i_cfg = Cfg(0))
      : Ising<Cfg>(len_in, i_cfg)
      , rng_gamma(Ising<Cfg>::get_N())
      , num_jumps(0ULL)
    {}
    //here self-loops ARE NOT counted as jumps
    void mcmove() { 
      if(Ising<Cfg>::mcmove()) {
        num_jumps++;
      }  
    }
    //here self-loops ARE counted as jumps
    void mcmove_lazy() {
      Ising<Cfg>::mcmove();
      num_jumps++;
    }
    void setcfg(Cfg i_cfg) { Ising<Cfg>::setcfg(i_cfg); num_jumps = 0LL; }
    double get_t() { return num_jumps > 0ULL ? rng_gamma(num_jumps):0.; }
};
#endif