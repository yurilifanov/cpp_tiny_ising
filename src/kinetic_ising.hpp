#ifndef _KINETIC_ISING_
#define _KINETIC_ISING_
#include <cstdint>          // int64_t
#include <vector>           // std::vector
#include "ising.hpp"        // Ising
#include "rng.hpp"          // RNG_exp
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
    void mcmove() { Ising<Cfg>::mcmove(); num_jumps++; }
    void setcfg(Cfg i_cfg) { Ising<Cfg>::setcfg(i_cfg); num_jumps = 0LL; }
    double get_t() { return num_jumps > 0ULL ? rng_gamma(num_jumps):0.; }
};
/*
  This implementation works too, but is slow because of the log operation
  in RNG_exp. The above implementation is much faster, given that get_t() is 
  rarely called.
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
    void mcmove() { Ising<Cfg>::mcmove(); t += rng_exp(); }
    void setcfg(Cfg i_cfg) { Ising<Cfg>::setcfg(i_cfg); t = 0.; }
    double get_t() const { return t; }
};
#endif