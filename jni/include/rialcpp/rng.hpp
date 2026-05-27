#ifndef RIALCPP_RNG_HPP
#define RIALCPP_RNG_HPP

#include <array>
#include <cstdint>

namespace rial
{

// Generate a float from raw random bits, using the standard FP32 configuration.
// `rnd` carries the random bit pattern (low bits used as each generator needs).

float gen_float12 (std::uint64_t rnd); // uniform in [1, 2)
float gen_float01 (std::uint64_t rnd); // uniform in [0, 1)  (32-bit fraction)
float gen_gaussian(std::uint64_t rnd); // N(0, 1)            (64-bit word)

// Squares32 counter-based uniform RNG (Widynski 2022): the 32-bit random word
// for (count, key).
std::uint32_t squares32(std::uint64_t count, std::uint64_t key);

// Approx-gaussian via Irwin-Hall: 12 uniform 32-bit words (low 32 bits of each
// used) summed and converted to a float.
float irwin_hall(const std::array<std::uint64_t, 12>& rnds);

} // rial
#endif//RIALCPP_RNG_HPP
