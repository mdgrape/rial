// The float-producing generators of RngSim.
//
// strict distribution tests are done in scala side. here we roughly check if
// the function generates the distribution.

#include <rialcpp/config.hpp>
#include <rialcpp/rng.hpp>

#include <doctest/doctest.h>

#include <array>
#include <random>
#include <vector>

#include <cmath>
#include <cstdint>

std::uint64_t mask(int w)
{
    return (w >= 64) ? ~std::uint64_t(0) : ((std::uint64_t(1) << w) - 1);
}

/** Relative tolerance of one ULP or so of `spec`.
 *
 *  The generators splice the random bits into a mantissa, and whether the bits
 *  that do not fit are dropped or rounded is the hardware's business. Either
 *  way the answer is within an ULP of the exact fraction, so that is what is
 *  checked -- the value, not the encoding. */
double ulp_tol(const rial::RealSpec& spec)
{
    return 2.0 * std::ldexp(1.0, -spec.manW);
}

/** The standard normal CDF, for the gen_gaussian histogram. */
double normal_cdf(double x)
{
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

TEST_CASE("gen_float12 is 1 + the mantissa the random bits spell out")
{
    const rial::RngSim rng;
    const int manW = rng.spec().manW;

    std::mt19937_64 gen(19660809);
    for (std::size_t i = 0; i < 1000; ++i)
    {
        const std::uint64_t rnd = gen();
        CAPTURE(rnd);

        // rng.hpp: "uniform in [1, 2); the low manW bits of rnd become the
        // mantissa". Bits above manW must not reach the result.
        const double ref = 1.0 + std::ldexp(static_cast<double>(rnd & mask(manW)), -manW);
        const double got = rng.gen_float12(rnd);

        CHECK(std::abs(got - ref) <= ulp_tol(rng.spec()));
        CHECK(got >= 1.0);
        CHECK(got <  2.0);
    }
}

TEST_CASE("gen_float01 is the random bits read as a fraction")
{
    const rial::RngSim rng;
    const int rndW = rng.config().float01W;

    // rng.hpp: "uniform in [0, 1) from the low float01W bits of rnd". Inputs
    // wider than that are left out on purpose -- they are not masked off on the
    // way in, so the result leaves [0, 1) (gen_float01(1 << 32) == 1.0).
    std::mt19937_64 gen(123456789);
    for (std::size_t i = 0; i < 1000; ++i)
    {
        const std::uint64_t rnd = gen() & mask(rndW);
        CAPTURE(rnd);

        const double ref = std::ldexp(static_cast<double>(rnd), -rndW);
        const double got = rng.gen_float01(rnd);

        CHECK(std::abs(got - ref) <= ulp_tol(rng.spec()) * std::max(ref, 0.0));
        CHECK(got >= 0.0);
        CHECK(got <  1.0);
    }

    SUBCASE("the ends")
    {
        CHECK(rng.gen_float01(0) == 0.0f);
        CHECK(rng.gen_float01(mask(rndW)) <  1.0f);
        CHECK(rng.gen_float01(mask(rndW)) >= 1.0f - ulp_tol(rng.spec()));
    }
}

TEST_CASE("irwin_hall is the 12 uniforms summed and centred")
{
    const rial::RngSim rng;
    const int rndW = rng.config().irwinHallW;

    std::mt19937_64 gen(19660809);
    for (std::size_t i = 0; i < 500; ++i)
    {
        std::array<std::uint64_t, 12> rnds{};
        double ref = -6.0;
        for (auto& r : rnds)
        {
            r = gen(); // the low irwinHallW bits of each are what is used
            ref += std::ldexp(static_cast<double>(r & mask(rndW)), -rndW);
        }

        const double got = rng.irwin_hall(rnds);
        CAPTURE(ref);

        // the sum spans [-6, 6], so an ULP there is the scale of the tolerance
        CHECK(std::abs(got - ref) <= 6.0 * ulp_tol(rng.spec()));
    }
}

TEST_CASE("gen_gaussian is distributed as a standard normal")
{
    const rial::RngSim rng;

    constexpr std::size_t n = 100000;

    std::mt19937_64 gen(19660809);
    std::vector<double> zs;
    zs.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        zs.push_back(rng.gen_gaussian(gen()));
    }

    SUBCASE("mean and variance")
    {
        double sum = 0.0;
        for (const double z : zs) {sum += z;}
        const double mean = sum / n;

        double sq = 0.0;
        for (const double z : zs) {sq += (z - mean) * (z - mean);}
        const double sd = std::sqrt(sq / (n - 1));

        // the standard error of each is ~1/sqrt(n) == 0.0032; the bound is
        // wide enough that only a real shift trips it
        CAPTURE(mean);
        CAPTURE(sd);
        CHECK(std::abs(mean) < 0.02);
        CHECK(std::abs(sd - 1.0) < 0.02);
    }

    SUBCASE("histogram")
    {
        // 32 bins over [-4, 4], plus one for each tail.
        constexpr int    nbins = 32;
        constexpr double lo = -4.0;
        constexpr double hi =  4.0;
        const     double w  = (hi - lo) / nbins;

        std::vector<std::size_t> count(nbins + 2, 0);
        for (const double z : zs)
        {
            if      (z <  lo) {++count[0];}
            else if (z >= hi) {++count[nbins + 1];}
            else              {++count[1 + static_cast<int>((z - lo) / w)];}
        }

        for (int i = 0; i < nbins + 2; ++i)
        {
            const double from = (i == 0)         ? -40.0 : lo + (i - 1) * w;
            const double to   = (i == nbins + 1) ?  40.0 : lo + i * w;
            const double expected = n * (normal_cdf(to) - normal_cdf(from));

            // a bin count is ~Poisson, so its spread is sqrt(expected); six of
            // those, plus a floor for the nearly empty tail bins
            CAPTURE(from);
            CAPTURE(to);
            CAPTURE(expected);
            CHECK(std::abs(count[i] - expected) <= 6.0 * std::sqrt(expected) + 5.0);
        }
    }
}
