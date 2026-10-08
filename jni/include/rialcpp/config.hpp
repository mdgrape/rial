#ifndef RIALCPP_CONFIG_HPP
#define RIALCPP_CONFIG_HPP

#include <array>
#include <limits>

#include <cstddef>
#include <cstdint>

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
// Which functions a MathSim supports.
//
// This is not just an access-control switch: the tables of the selected
// functions share one merged set of coefficient/intermediate bit widths, the
// same way the hardware shares one polynomial evaluator.  A `simple` MathSim
// therefore returns slightly different values (~1 ULP) from an `all` one, and
// matching a given MathFunctions module means giving MathSim the same function
// list the module was built with.
//
// Functions pull in what they need: selecting `acos` also gives you `sqrt`
// (its first phase reuses the sqrt table) and `atan2` gives you `reciprocal`,
// exactly as in the hardware.  `MathSim::config().funcs` reports the effective
// set after those additions.

enum class Func : std::uint32_t
{
    none       = 0,
    sqrt       = 1u << 0,
    invsqrt    = 1u << 1,
    reciprocal = 1u << 2,
    sin        = 1u << 3,
    cos        = 1u << 4,
    acos       = 1u << 5,
    atan2      = 1u << 6,
    exp        = 1u << 7,
    log        = 1u << 8,
    sigmoid    = 1u << 9,
    softplus   = 1u << 10,
    scale_mixture_gaussian = 1u << 11,

    // The presets of the `MathFuncConfig` companion object on the Scala side.

    // the functions that need no pre/post-process multiplier
    simple   = sqrt | invsqrt | reciprocal,
    // the frequently used ones
    standard = simple | sin | cos | acos | atan2 | exp | log,
    // everything, including the task-specific ones
    all      = standard | sigmoid | softplus | scale_mixture_gaussian,
};

// The number of functions `Func` names, i.e. the most per-function table
// address widths a configuration can carry.
inline constexpr std::size_t func_count = 12;

constexpr Func operator|(Func lhs, Func rhs) noexcept
{
    return static_cast<Func>(static_cast<std::uint32_t>(lhs) |
                             static_cast<std::uint32_t>(rhs));
}
constexpr Func operator&(Func lhs, Func rhs) noexcept
{
    return static_cast<Func>(static_cast<std::uint32_t>(lhs) &
                             static_cast<std::uint32_t>(rhs));
}
constexpr Func& operator|=(Func& lhs, Func rhs) noexcept
{
    lhs = lhs | rhs;
    return lhs;
}
/** true if every function in `f` is contained in `set`. */
constexpr bool has(Func set, Func f) noexcept
{
    return (set & f) == f;
}

// -----------------------------------------------------------------------------
// The width of the address of one function's table.  Mirrors an entry of
// `rial.math.PerFuncAddrWidth`.
//
// A function that reads the table of another one (`cos` reads the table of
// `sin`) is read at that table's width, so the two entries have to agree.

struct FuncAddrWidth
{
    Func func;  // the function whose table is read at this width
    int  addrW; // the width of its address, in bits
};

// The per-function widths of the FP32 presets below: each table is read at the
// narrowest address that keeps its maximum error -- sqrt 5 bits, 1/x and acos
// 7, the other standard functions 6 -- and every table left out at
// `MathFuncConfig::addrW`, 8 bits.  Mirrors `MathFuncConfig.standardFP32`.
//
// The narrowest of them decides the shared `inBitW`: a position that reaches
// the end of the fraction of every table is `manW` less that width.

constexpr std::array<FuncAddrWidth, func_count> fp32_addr_widths() noexcept
{
    return {{
        {Func::sqrt,       5}, {Func::invsqrt, 6}, {Func::reciprocal, 7},
        {Func::sin,        6}, {Func::cos,     6}, {Func::acos,       7},
        {Func::atan2,      6}, {Func::exp,     6}, {Func::log,        6},
    }};
}
inline constexpr int fp32_min_addrW = 5;

// -----------------------------------------------------------------------------
// Table / polynomial configuration of the approximated math functions.  Mirrors
// `rial.math.MathFuncConfig`, so a MathSim built from the same numbers as a
// hardware `MathFunctions` module reproduces its results bit-exactly.

struct MathFuncConfig
{
    Func funcs;   // which functions to support
    int  nOrder;  // order of the approximating polynomial
    int  addrW;   // table address width (== log2 of the interval count)
    int  inBitW;  // width of the in-interval polynomial input
    int  outBitW; // width of the polynomial output ((1 << outBitW) means 1.0)

    // The tables that are read at a width of their own; every other one is
    // read at `addrW`.  At most one entry per function, in any order; an entry
    // of `Func::none` is a hole and is skipped.
    std::array<FuncAddrWidth, func_count> addrWs = {};

    // (sigmaA, sigmaB) for Func::scale_mixture_gaussian; NaN takes the Scala
    // side's default, `MathFuncConfig.defaultSMGSigma` == (exp(-1), exp(-6)).
    double smgSigmaA = std::numeric_limits<double>::quiet_NaN();
    double smgSigmaB = std::numeric_limits<double>::quiet_NaN();

    // The width of the address of the table `f` reads.
    constexpr int addrW_of(Func f) const noexcept
    {
        for(const auto& w : this->addrWs)
        {
            if(w.func == Func::none) {continue;} // a hole
            if(w.func == f) {return w.addrW;}
        }
        return this->addrW;
    }

    // The parameters of `MathFuncConfig.standardFP32`: a 2nd order polynomial
    // with 3 bits of output below the mantissa, over the tables of
    // `fp32_addr_widths`.
    static constexpr MathFuncConfig fp32(Func funcs = Func::standard) noexcept
    {
        return MathFuncConfig{funcs, 2, 8,
            RealSpec::float32().manW - fp32_min_addrW,
            RealSpec::float32().manW + 3, fp32_addr_widths()};
    }
    // The parameters of `MathFuncConfig.standardBF16`: a plain lookup table
    // (nOrder = 0) addressed by the whole mantissa.
    static constexpr MathFuncConfig bf16(Func funcs = Func::standard) noexcept
    {
        return MathFuncConfig{funcs, 0, RealSpec::bfloat16().manW, 0,
            RealSpec::bfloat16().manW + 1};
    }

    // The parameters of `MathFuncConfig.allFP32` / `allBF16`: Func::all needs a
    // wider polynomial output, driven by scale_mixture_gaussian (and softplus
    // in BF16).  Reaching for `fp32(Func::all)` instead collapses the SMG
    // result -- see the scaladoc of MathFuncConfig.allFP32.
    static constexpr MathFuncConfig fp32_all() noexcept
    {
        return MathFuncConfig{Func::all, 2, 8,
            RealSpec::float32().manW - fp32_min_addrW,
            RealSpec::float32().manW + 8, fp32_addr_widths()};
    }
    static constexpr MathFuncConfig bf16_all() noexcept
    {
        return MathFuncConfig{Func::all, 0, RealSpec::bfloat16().manW, 0,
            RealSpec::bfloat16().manW + 10};
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
