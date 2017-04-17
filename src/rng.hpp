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
/*  
  Random integer generator
*/
class RNG_ull {
    uint64_t s0; 
    uint64_t s1;
    uint64_t rotl(const uint64_t x, int k) const { 
      return (x << k) | (x >> (64 - k)); 
    }
  public:
    RNG_ull();
    uint64_t operator()();
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
class RNG_d : public RNG_ull {
    static constexpr uint64_t dblprfx = 0x3ffULL << 52U;
  public:
    RNG_d() : RNG_ull() {}
    double operator()() {
      union DUll {
        double d;
        uint64_t ull;
      } val;
      val.ull = dblprfx | this->RNG_ull::operator()() >> 12;
      if(is_open) {
        val.ull |= 1ULL;
      }
      return val.d - 1.0;
    }
};
/*  
  Exponential random number generator
*/
class RNG_exp : public RNG_d<true> {
    const double inv_lambda;
  public:
    RNG_exp(double lambda) : RNG_d<true>(), inv_lambda(1. / lambda) {}
    double operator()() {
      return -inv_lambda * log(this->RNG_d<true>::operator()());
    }
};
#endif