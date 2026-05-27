#ifndef RIALCPP_RNG_HPP
#define RIALCPP_RNG_HPP

#include <cstdint>

namespace rial
{

// Generate a float from raw random bits, using the standard FP32 configuration.
// `rnd` carries the random bit pattern (low bits used as each generator needs).

float gen_float12 (std::uint64_t rnd); // uniform in [1, 2)
float gen_float01 (std::uint64_t rnd); // uniform in [0, 1)  (32-bit fraction)
float gen_gaussian(std::uint64_t rnd); // N(0, 1)            (64-bit word)

} // rial
#endif//RIALCPP_RNG_HPP
