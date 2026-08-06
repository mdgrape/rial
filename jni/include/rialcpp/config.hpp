#ifndef RIALCPP_CONFIG_HPP
#define RIALCPP_CONFIG_HPP

namespace rial
{

// -----------------------------------------------------------------------------
// Floating point format the software emulation works in.  Mirrors
// `rial.arith.RealSpec`; the sign / NaN / subnormal switches are not exposed
// because every preset on the Scala side keeps them at their defaults
// (sign and NaN enabled, subnormals disabled).
//
// Note that the C++ interface is `float`, so the emulated format has to be
// representable in an IEEE-754 binary32: inputs are rounded into the emulated
// format on the way in, and results are widened back on the way out.

struct RealSpec
{
    int exW;    // exponent width
    int exBias; // exponent bias
    int manW;   // mantissa width

    static constexpr RealSpec float32 () noexcept {return RealSpec{8, 0x7F, 23};}
    static constexpr RealSpec float16 () noexcept {return RealSpec{5, 0x0F, 10};}
    static constexpr RealSpec bfloat16() noexcept {return RealSpec{8, 0x7F,  7};}
};

// -----------------------------------------------------------------------------
// Table / polynomial configuration of the approximated math functions.  Mirrors
// `rial.math.MathFuncConfig` over the `standard` function set, so a MathSim
// built from the same numbers as a hardware `MathFunctions` module reproduces
// its results bit-exactly.

struct MathFuncConfig
{
    int nOrder;  // order of the approximating polynomial
    int addrW;   // table address width (== log2 of the interval count)
    int inBitW;  // width of the in-interval polynomial input
    int outBitW; // width of the polynomial output ((1 << outBitW) means 1.0)

    // The parameters of `MathFuncConfig.standardFP32`: (nOrder, addrW,
    // extraBits) = (2, 8, 3).
    static constexpr MathFuncConfig fp32() noexcept
    {
        return MathFuncConfig{2, 8,
            RealSpec::float32().manW - 8, RealSpec::float32().manW + 3};
    }
    // The parameters of `MathFuncConfig.standardBF16`: a plain lookup table
    // (nOrder = 0) addressed by the whole mantissa.
    static constexpr MathFuncConfig bf16() noexcept
    {
        return MathFuncConfig{0, RealSpec::bfloat16().manW, 0,
            RealSpec::bfloat16().manW + 1};
    }
};

// -----------------------------------------------------------------------------
// Configuration of the float-producing random number generators.  Mirrors the
// parameters `rial.rng.RngSim` passes on to the individual generator models.
//
// The raw bit generators (squares32, threefry4x32) are not configured here:
// they produce integers and know nothing about the floating point format, so
// they stay free functions and take what they need (the Threefry round count)
// as a call argument.  See rialcpp/rng.hpp.

struct RngConfig
{
    int float01W;            // width of the random fraction gen_float01 consumes
    int gaussianW;           // width of the random word gen_gaussian consumes
    int gaussianOrder;       // polynomial order of the probit approximation
    int gaussianLumpOctaves; // deepest tail octaves merged into one interval
    int irwinHallW;          // width of each of the 12 words irwin_hall sums

    static constexpr RngConfig standard() noexcept
    {
        return RngConfig{32, 64, 3, 2, 32};
    }
};

} // rial
#endif//RIALCPP_CONFIG_HPP
