/*
  This is a class wrapper around the C code for xoroshiro128+ generator
  available at http://xoroshiro.di.unimi.it/ (Accessed 2017-04-15).
*/
#ifndef _RNG_
#define _RNG_
#include <cstdint>
#include <random>
#include <ostream>
#include <iostream>
#include <utility> 
#include <cmath>
#include <limits>
/*  
  Random integer generator
*/
class RNG_ull {
    using result_type = uint64_t;
    using limits = std::numeric_limits<result_type>;
    uint64_t s0; 
    uint64_t s1;
    uint64_t rotl(const uint64_t x, int k) const { 
      return (x << k) | (x >> (64 - k)); 
    }
  public:
    RNG_ull();
    result_type operator()();
    static constexpr result_type min() { return limits::min(); }
    static constexpr result_type max() { return limits::max(); }
};
inline RNG_ull::RNG_ull() {
  std::random_device rd;
  union TmpUnion {
    uint32_t ul[2];
    uint64_t ull;
  } vals;
  uint64_t s[2];
  for(uint8_t k = 0; k < 2U; k++) {
    for(uint8_t i = 0U; i < 2U; i++) {
      uint32_t val = rd();
      while(!val) {
        val = rd();
      }
      vals.ul[i] = val;
    }
    s[k] = vals.ull;
  }
  s0 = s[0];
  s1 = s[1];
}
inline uint64_t RNG_ull::operator()() {
  const uint64_t result = s0 + s1;
  s1 ^= s0;
  s0 = rotl(s0, 55) ^ s1 ^ (s1 << 14);
  s1 = rotl(s1, 36);
  return result;
}
/*  
  Random double generator, template argument specifies interval type:
    - is_open = false : generates double in [0, 1)
    - is_open = true : generates double in (0, 1)
*/
template <bool is_open = false>
class RNG_d {
    RNG_ull rng_ull;
    static constexpr uint64_t dblprfx = 0x3ffULL << 52U;
  public:
    RNG_d() {}
    double operator()() {
      union DUll {
        double d;
        uint64_t ull;
      } val;
      val.ull = dblprfx | rng_ull() >> 12;
      if(is_open) {
        val.ull |= 1ULL;
      }
      return val.d - 1.0;
    }
};
/*  
  Exponential random number generator
*/
class RNG_exp {
    RNG_d<true> rng_d;
    const double inv_lambda;
  public:
    RNG_exp(double lambda) : rng_d(), inv_lambda(1. / lambda) {}
    double operator()() {
      return -inv_lambda * log(rng_d());
    }
};
/*  
  Gamma random number generator for sampling sums of exponential random
  variates.
*/
class RNG_gamma {
    RNG_ull rng_ull;
    using Gamma = std::gamma_distribution<double>;
    const double beta; //beta = 1 / theta
    Gamma gamma;
  public:
    RNG_gamma(double theta) 
      : rng_ull()
      , beta(1. / theta)
      , gamma(Gamma(1, beta)) 
    {}
    double operator()(uint64_t n) {
      gamma.param(Gamma::param_type(n, beta));
      return gamma(rng_ull);
    }
};
#endif