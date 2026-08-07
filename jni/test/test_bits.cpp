// rial::to_bits / to_double / to_float against the encoding they are supposed
// to implement.
//
// Where the RealSpec is one of the ordinary IEEE formats the answer is fixed:
// the bit pattern is the one std::bit_cast produces, and anything else is a
// binding or a bias bug. Where the spec is narrower there is a choice of
// rounding, so what is checked is the part that no rounding mode can move: a
// value the format represents exactly has to come back exactly, and a value it
// does not has to land on one of the two patterns bracketing it.

#include <rialcpp/config.hpp>
#include <rialcpp/real.hpp>

#include <doctest/doctest.h>

#include <limits>
#include <random>
#include <vector>

#include <bit>
#include <cmath>
#include <cstdint>

namespace
{

std::uint64_t mask(int w)
{
    return (w >= 64) ? ~std::uint64_t(0) : ((std::uint64_t(1) << w) - 1);
}

/** The IEEE binary64 spec, for which to_bits is plain std::bit_cast. */
rial::RealSpec float64()
{
    return rial::RealSpec{11, 0x3FF, 52};
}

/** The value `spec` encodes in `bits`, assuming a normal number.  This is the
 *  definition to_double is measured against, so it is spelled out rather than
 *  going through the library. */
double decode(const rial::RealSpec& spec, std::uint64_t bits)
{
    const int           ex  = static_cast<int>((bits >> spec.manW) & mask(spec.exW));
    const std::uint64_t man = bits & mask(spec.manW);

    const double v = std::ldexp(
        1.0 + std::ldexp(static_cast<double>(man), -spec.manW), ex - spec.exBias);

    return ((bits >> (spec.exW + spec.manW)) & 1) ? -v : v;
}

/** A full precision value whose truncation into `spec` is `bits`: the same sign
 *  and exponent, with the mantissa bits below the format's last one filled from
 *  `low`.  Building the value this way keeps it inside the narrow format's
 *  exponent range, which a double drawn at random almost never is. */
double fill_below(const rial::RealSpec& spec, std::uint64_t bits, std::uint64_t low)
{
    const std::uint64_t b = std::bit_cast<std::uint64_t>(decode(spec, bits));
    return std::bit_cast<double>(b | (low & mask(52 - spec.manW)));
}

/** Bit patterns of `spec`, all of them normal numbers. */
std::vector<std::uint64_t> normal_patterns(const rial::RealSpec& spec,
                                           std::size_t n, std::uint32_t seed)
{
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::uint64_t> bit(0, mask(rial::bit_width(spec)));

    std::vector<std::uint64_t> patterns;
    while (patterns.size() < n)
    {
        const std::uint64_t b  = bit(rng);
        const std::uint64_t ex = (b >> spec.manW) & mask(spec.exW);
        if (ex == 0 || ex == mask(spec.exW)) {continue;} // zero, subnormal, inf, nan
        patterns.push_back(b);
    }
    return patterns;
}

} // namespace

TEST_CASE("to_bits and to_float are bit_cast in FP32")
{
    const rial::RealSpec spec = rial::RealSpec::float32();

    SUBCASE("the values a float takes")
    {
        for (const std::uint64_t b : normal_patterns(spec, 1000, 20220809))
        {
            const float x = std::bit_cast<float>(static_cast<std::uint32_t>(b));
            CAPTURE(x);
            CHECK(rial::to_bits(spec, x) == b);
            CHECK(rial::to_float(spec, b) == x);
            // to_double widens the same number, so it has to agree
            CHECK(rial::to_double(spec, b) == static_cast<double>(x));
            // and the double overload of to_bits must not take a detour
            CHECK(rial::to_bits(spec, static_cast<double>(x)) == b);
        }
    }

    SUBCASE("infinity and NaN")
    {
        for (const float x : {std::numeric_limits<float>::infinity(),
                             -std::numeric_limits<float>::infinity(),
                              std::numeric_limits<float>::quiet_NaN()})
        {
            CAPTURE(x);
            const auto bits = std::bit_cast<std::uint32_t>(x);
            CHECK(rial::to_bits(spec, x) == bits);
            // NaN never compares equal, so the pattern is what is compared
            CHECK(std::bit_cast<std::uint32_t>(rial::to_float(spec, bits)) == bits);
        }
    }

    SUBCASE("signed zero")
    {
        // to_bits keeps the sign of a zero, as bit_cast does.
        CHECK(rial::to_bits(spec, 0.0f)  == 0u);
        CHECK(rial::to_bits(spec, -0.0f) == 0x80000000u);

        // Decoding only gets it back through the FP32 fast path of
        // RealGeneric.toFloat; to_double returns an unsigned zero in every
        // format, so only the value is checked here.
        CHECK(std::bit_cast<std::uint32_t>(rial::to_float(spec, 0x80000000u))
              == 0x80000000u);
        CHECK(rial::to_double(spec, 0x80000000u) == 0.0);
    }

    SUBCASE("subnormals flush to zero")
    {
        // The presets keep subnormals off, so this is the one place the
        // encoding departs from bit_cast -- deliberately.
        for (const float x : {std::numeric_limits<float>::denorm_min(),
                              std::numeric_limits<float>::min() / 2.0f})
        {
            CAPTURE(x);
            CHECK(rial::to_bits(spec, x) == 0u);
            CHECK(rial::to_bits(spec, -x) == (std::uint64_t(1) << 31));
        }
        // the smallest normal, on the other hand, is still exact
        const float lo = std::numeric_limits<float>::min();
        CHECK(rial::to_bits(spec, lo) == std::bit_cast<std::uint32_t>(lo));
    }
}

TEST_CASE("to_bits and to_double are bit_cast in FP64")
{
    const rial::RealSpec spec = float64();

    SUBCASE("the values a double takes")
    {
        for (const std::uint64_t b : normal_patterns(spec, 1000, 19660809))
        {
            const double x = std::bit_cast<double>(b);
            CAPTURE(x);
            CHECK(rial::to_bits(spec, x) == b);
            CHECK(rial::to_double(spec, b) == x);
        }
    }

    SUBCASE("infinity and NaN")
    {
        for (const double x : {std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()})
        {
            CAPTURE(x);
            const auto bits = std::bit_cast<std::uint64_t>(x);
            CHECK(rial::to_bits(spec, x) == bits);
            CHECK(std::bit_cast<std::uint64_t>(rial::to_double(spec, bits)) == bits);
        }
    }

    SUBCASE("signed zero")
    {
        const auto neg = std::uint64_t(1) << 63;
        CHECK(rial::to_bits(spec, 0.0)  == 0u);
        CHECK(rial::to_bits(spec, -0.0) == neg);
        // as in FP32: the value comes back, the sign of the zero does not
        CHECK(rial::to_double(spec, neg) == 0.0);
    }
}

TEST_CASE("a narrower format keeps what it can represent exactly")
{
    // Values taken from the format itself: whatever the rounding mode, a value
    // the format already holds must come back as the pattern it came from.
    for (const rial::RealSpec spec : {rial::RealSpec::bfloat16(),
                                      rial::RealSpec::float16()})
    {
        CAPTURE(spec.exW);
        CAPTURE(spec.manW);

        for (const std::uint64_t b : normal_patterns(spec, 500, 20260807))
        {
            const double x = decode(spec, b);
            CAPTURE(x);

            CHECK(rial::to_double(spec, b) == x);
            CHECK(rial::to_float (spec, b) == static_cast<float>(x));

            CHECK(rial::to_bits(spec, x) == b);
            CHECK(rial::to_bits(spec, static_cast<float>(x)) == b);
        }
    }
}

TEST_CASE("a narrower format rounds to one of the two patterns around the value")
{
    // Full-precision values, which the format generally does not hold. Every
    // rounding mode -- to nearest either way, toward zero, toward an infinity
    // -- picks one of the two neighbours, so the pattern is the truncation or
    // the one above it. That is 1 ULP of the narrow format wide, and it is as
    // tight as this can be without pinning the mode down.
    std::mt19937_64 rng(20260808);

    for (const rial::RealSpec spec : {rial::RealSpec::bfloat16(),
                                      rial::RealSpec::float16()})
    {
        CAPTURE(spec.exW);
        CAPTURE(spec.manW);

        for (const std::uint64_t b : normal_patterns(spec, 500, 20260808))
        {
            // b is the truncation of x by construction
            const double x = fill_below(spec, b, rng());
            CAPTURE(x);

            const std::uint64_t got = rial::to_bits(spec, x);
            CHECK(got >= b);
            CHECK(got <= b + 1);

            // the float overload lands in the same pair: narrowing x to a
            // float first moves it by well under half an ULP of either format
            const std::uint64_t via_float = rial::to_bits(spec, static_cast<float>(x));
            CHECK(via_float >= b);
            CHECK(via_float <= b + 1);
        }
    }
}

TEST_CASE("BF16 widens into FP32 by filling the low mantissa bits with zero")
{
    // BF16 is the top half of a float32 -- same exponent width, same bias --
    // so widening is a shift and nothing else.
    const rial::RealSpec spec = rial::RealSpec::bfloat16();

    for (const std::uint64_t b : normal_patterns(spec, 500, 20260809))
    {
        CAPTURE(b);
        const auto widened = static_cast<std::uint32_t>(b) << 16;
        CHECK(std::bit_cast<std::uint32_t>(rial::to_float(spec, b)) == widened);

        // and narrowing back drops exactly those bits again
        CHECK(rial::to_bits(spec, std::bit_cast<float>(widened)) == b);
    }
}
