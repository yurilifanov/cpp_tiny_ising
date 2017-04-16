#ifndef _CFG_TYPES_
#define _CFG_TYPES_
#include <cstdint>      // int64_t
#include <limits>       // std::numeric_limits
#include <bitset>       // std::bitset
/*   TO DO:
 *     - Overload appropriate operators for arbitrary cfg type
 *     - e.g. __int128, __m256i, __m128i etc 
**/
template <typename T>
inline constexpr int num_digits() {
  return std::numeric_limits<T>::digits;
}

template <typename T>
inline int popcount(T val) noexcept {
  return std::bitset<num_digits<T>()>(val).count();
}
#endif
