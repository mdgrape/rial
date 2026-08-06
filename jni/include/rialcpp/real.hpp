#ifndef RIALCPP_REAL_HPP
#define RIALCPP_REAL_HPP

#include <rialcpp/config.hpp>

#include <cstdint>

namespace rial
{

// -----------------------------------------------------------------------------
// Conversion between the raw bit pattern of a RealSpec and an ordinary number.
//
// The bit-pattern interfaces (MathSim::sqrt_bits, RngSim::gen_float12_bits,
// ...) exist for formats a `float` cannot carry.  These turn such a bit pattern
// back into something printable, and round a number into the format.
//
// The spec is not part of the encoding, so it has to be passed in; MathSim and
// RngSim forward their own spec through `to_double(bits)` / `to_bits(x)`.
//
// Formats wider than 64 bits are out of reach here, and a `double` cannot
// represent more than 52 mantissa bits -- but a table-based approximation is
// not what you would use at that precision anyway.

/** the value of the number `spec` encodes in `bits` (bits above `spec` width
 *  are ignored). */
double to_double(const RealSpec& spec, std::uint64_t bits);
float  to_float (const RealSpec& spec, std::uint64_t bits);

/** `x` rounded into `spec`, as a bit pattern. */
std::uint64_t to_bits(const RealSpec& spec, double x);
std::uint64_t to_bits(const RealSpec& spec, float  x);

/** total width of the encoding: sign + exponent + mantissa. */
constexpr int bit_width(const RealSpec& spec) noexcept
{
    return 1 + spec.exW + spec.manW;
}

} // rial
#endif//RIALCPP_REAL_HPP
