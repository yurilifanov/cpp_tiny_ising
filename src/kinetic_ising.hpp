#ifndef _KINETIC_ISING_
#define _KINETIC_ISING_
#include <cstdint>          // int64_t
#include <vector>           // std::vector
#include "ising.hpp"        // Ising
#include "rng.hpp"          // RNG_exp
template <typename Cfg>
class KineticIsing : public Ising<Cfg> {
    using v_int64 = std::vector<int64_t>;
    RNG_exp rng_exp;
    double t;
  public:
    KineticIsing(v_int64 len_in, Cfg i_cfg = Cfg(0))
      : Ising<Cfg>(len_in, i_cfg)
      , rng_exp(Ising<Cfg>::get_N())
      , t(0.)
    {}
    void mcmove() { Ising<Cfg>::mcmove(); t += rng_exp(); }
    double get_t() const { return t; }
};
#endif