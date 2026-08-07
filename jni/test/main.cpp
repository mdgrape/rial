#include <rialcpp/math.hpp>
#include <rialcpp/rng.hpp>

#include <array>
#include <iomanip>
#include <iostream>

#include <cmath>
#include <cstdint>

template<typename Fref, typename Frial>
void check(std::string name, float from, float to,
        Fref f_ref, Frial f_rial)
{
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << name << std::endl;
    std::cout << "x libm rial diff" << std::endl;
    for(std::size_t i=0; i<100; ++i)
    {
        const float x = 0.01f * (from * i + to * (100 - i));
        const auto libm = f_ref(x);
        const auto rial = f_rial(x);
        const auto diff = std::abs(libm - rial);

        std::cout << std::setw(15) << std::setprecision(10) << x    << " "
                  << std::setw(15) << std::setprecision(10) << libm << " "
                  << std::setw(15) << std::setprecision(10) << rial << " "
                  << diff << std::endl;
    }
}

int main()
{
    constexpr float pi = 3.14159265;

    // The approximation depends on the float format and on the table
    // parameters, so pick a configuration and keep the instance around: the
    // tables are generated once, when it is constructed.
    const rial::MathSim math; // == MathSim(RealSpec::float32(), MathFuncConfig::fp32())

    check("sqrt", 0.0, 4.0,
        [](float x) {return std::sqrt(x); },
        [&](float x) {return math.sqrt(x);}
    );

    check("invsqrt", 0.5, 4.0,
        [](float x) {return 1.0f / std::sqrt(x); },
        [&](float x) {return math.invsqrt(x);}
    );

    check("reciprocal", 0.5, 4.0,
        [](float x) {return 1.0f / x; },
        [&](float x) {return math.reciprocal(x);}
    );

    check("sin", -pi, pi,
        [](float x) {return std::sin(x); },
        [&](float x) {return math.sin(x);}
    );

    check("cos", -pi, pi,
        [](float x) {return std::cos(x); },
        [&](float x) {return math.cos(x);}
    );

    check("acos", -1.0, 1.0,
        [](float x) {return std::acos(x); },
        [&](float x) {return math.acos(x);}
    );


    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "atan2" << std::endl;
    std::cout << "x y libm rial diff" << std::endl;
    for(std::size_t i=0; i<100; ++i)
    {
    for(std::size_t j=0; j<100; ++j)
    {
        const float x = 0.01f * (-1.0f * i + 1.0f * (100 - i));
        const float y = 0.01f * (-1.0f * j + 1.0f * (100 - j));
        const auto libm = std::atan2(y, x);
        const auto rial = math.atan2(y, x);
        const auto diff = std::abs(libm - rial);

        std::cout << std::setw(15) << std::setprecision(10) << x    << " "
                  << std::setw(15) << std::setprecision(10) << y    << " "
                  << std::setw(15) << std::setprecision(10) << libm << " "
                  << std::setw(15) << std::setprecision(10) << rial << " "
                  << diff << std::endl;
    }
    }

    check("exp", 0.1, 3.1,
        [](float x) {return std::exp(x); },
        [&](float x) {return math.exp(x);}
    );

    check("log", 0.1, 3.1,
        [](float x) {return std::log(x); },
        [&](float x) {return math.log(x);}
    );

    // -----------------------------------------------------------------------
    // Other configurations, living alongside the FP32 one. Each instance keeps
    // its own tables and none of them shares any state.

    // BFloat16 with a plain lookup table (nOrder = 0).
    const rial::MathSim bf16(rial::RealSpec::bfloat16(),
                             rial::MathFuncConfig::bf16());

    // FP32 again, but supporting sqrt/invsqrt/reciprocal only. The tables of
    // the supported functions share one merged set of bit widths, so a narrower
    // function set gives (slightly) different results -- which is the point:
    // this is what a MathFunctions module built with `simple` produces.
    const rial::MathSim simple(rial::RealSpec::float32(),
                               rial::MathFuncConfig::fp32(rial::Func::simple));

    // ... and one supporting everything, which needs a wider polynomial output
    // (see MathFuncConfig::fp32_all).
    const rial::MathSim all(rial::RealSpec::float32(),
                            rial::MathFuncConfig::fp32_all());

    // sqrt alone: the shared widths now come from the sqrt table only.
    const rial::MathSim sqrt_only(rial::RealSpec::float32(),
                                  rial::MathFuncConfig::fp32(rial::Func::sqrt));

    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "sqrt under different configurations" << std::endl;
    std::cout << "x libm standard sqrt-only simple all bf16" << std::endl;
    for(std::size_t i=1; i<=10; ++i)
    {
        const float x = 0.4f * i;
        std::cout << std::setw(15) << std::setprecision(10) << x                 << " "
                  << std::setw(15) << std::setprecision(10) << std::sqrt(x)      << " "
                  << std::setw(15) << std::setprecision(10) << math.sqrt(x)      << " "
                  << std::setw(15) << std::setprecision(10) << sqrt_only.sqrt(x) << " "
                  << std::setw(15) << std::setprecision(10) << simple.sqrt(x)    << " "
                  << std::setw(15) << std::setprecision(10) << all.sqrt(x)       << " "
                  << std::setw(15) << std::setprecision(10) << bf16.sqrt(x)      << std::endl;
    }

    // How often each configuration's sqrt disagrees with the `standard` one.
    // Not every function set changes the widths -- `simple` already contains
    // reciprocal, whose table dominates the merge in FP32, so it matches
    // `standard` exactly; sqrt on its own does not.
    {
        constexpr std::size_t n = 20000;
        const auto count = [&](const rial::MathSim& sim)
        {
            std::size_t differ = 0;
            for(std::size_t i=0; i<=n; ++i)
            {
                const float x = 4.0f * i / n;
                if(sim.sqrt(x) != math.sqrt(x)) {++differ;}
            }
            return differ;
        };
        std::cout << "sqrt differing from `standard` over " << n+1 << " inputs: "
                  << "sqrt-only=" << count(sqrt_only)
                  << ", simple="  << count(simple)
                  << ", all="     << count(all) << std::endl;
    }

    // functions outside the configured set are rejected, and a function pulled
    // in as a dependency (acos needs the sqrt table) becomes available too.
    const rial::MathSim acos_only(rial::RealSpec::float32(),
                                  rial::MathFuncConfig::fp32(rial::Func::acos));
    std::cout << "acos-only config: supports sqrt = "
              << acos_only.supports(rial::Func::sqrt) << ", supports exp = "
              << acos_only.supports(rial::Func::exp)  << std::endl;
    try
    {
        std::cout << acos_only.exp(1.0f) << std::endl;
        std::cout << "FAILED: exp should not be available" << std::endl;
        return 1;
    }
    catch(const std::exception& e)
    {
        std::cout << "  rejected as expected: " << e.what() << std::endl;
    }

    std::cout << "sigmoid(0.5)  = " << all.sigmoid(0.5f)  << std::endl;
    std::cout << "softplus(0.5) = " << all.softplus(0.5f) << std::endl;
    std::cout << "smg(0.5)      = " << all.scale_mixture_gaussian(0.5f) << std::endl;

    // -----------------------------------------------------------------------
    // The bit-pattern interface, for formats a float cannot carry. The spec is
    // the instance's own, so to_bits()/to_double() get you in and out.

    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "bit-pattern interface" << std::endl;

    // it agrees with the float interface where both apply
    if (math.to_float(math.sqrt_bits(math.to_bits(2.0f))) != math.sqrt(2.0f))
    {
        std::cout << "FAILED: sqrt_bits disagrees with sqrt in FP32" << std::endl;
        return 1;
    }
    std::cout << "fp32: sqrt_bits(0x" << std::hex << math.to_bits(2.0f) << ") = 0x"
              << math.sqrt_bits(math.to_bits(2.0f)) << std::dec
              << " = " << math.to_double(math.sqrt_bits(math.to_bits(2.0f)))
              << " (matches sqrt(2.0f))" << std::endl;

    // an FP64-shaped format: a float would throw most of this away
    const rial::MathSim fp64(rial::RealSpec{11, 0x3FF, 52},
        rial::MathFuncConfig{rial::Func::sqrt, /*nOrder=*/2, /*addrW=*/8,
                             /*inBitW=*/44, /*outBitW=*/55});

    std::cout << "fp64-shaped (exW=11, manW=52), sqrt:" << std::endl;
    std::cout << "x libm rial diff" << std::endl;
    for(const double x : {0.5, 2.0, 3.0, 10.0, 1234.5})
    {
        const double r = fp64.to_double(fp64.sqrt_bits(fp64.to_bits(x)));
        std::cout << std::setw(15) << std::setprecision(17) << x         << " "
                  << std::setw(24) << std::setprecision(17) << std::sqrt(x) << " "
                  << std::setw(24) << std::setprecision(17) << r         << " "
                  << std::abs(std::sqrt(x) - r) << std::endl;
    }

    // -----------------------------------------------------------------------
    // RNG generators: input is raw random bits, output is the generated float.
    const rial::RngSim rng; // == RngSim(RealSpec::float32(), RngConfig::standard())

    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "rng" << std::endl;
    std::cout << "gen_float12(0)          = " << rng.gen_float12(0)             << " (expect 1)"    << std::endl;
    std::cout << "gen_float01(0x80000000) = " << rng.gen_float01(0x80000000ull) << " (expect 0.5)"  << std::endl;
    std::cout << "gen_float01(0xC0000000) = " << rng.gen_float01(0xC0000000ull) << " (expect 0.75)" << std::endl;
    for (std::uint64_t s = 1; s <= 5; ++s)
    {
        const std::uint64_t bits = s * 0x9E3779B97F4A7C15ull; // arbitrary spread
        std::cout << "gen_gaussian(0x" << std::hex << bits << std::dec << ") = "
                  << rng.gen_gaussian(bits) << std::endl;
    }

    const std::uint64_t key = 0x2467cb532b5ce8d1ull;
    for (std::uint64_t cnt = 0; cnt < 5; ++cnt)
    {
        std::cout << "squares32(count=" << cnt << ", key) = 0x" << std::hex
                  << rial::squares32(cnt, key) << std::dec << std::endl;
    }

    // feed 12 squares32 words into Irwin-Hall to draw an approx-gaussian sample.
    std::array<std::uint64_t, 12> rnds{};
    for (std::uint64_t i = 0; i < 12; ++i) { rnds[i] = rial::squares32(i, key); }
    std::cout << "irwin_hall(squares32 x12) = " << rng.irwin_hall(rnds) << std::endl;

    // the generators have a bit-pattern interface too
    if (rng.to_float(rng.gen_gaussian_bits(0x9E3779B97F4A7C15ull)) !=
        rng.gen_gaussian(0x9E3779B97F4A7C15ull))
    {
        std::cout << "FAILED: gen_gaussian_bits disagrees with gen_gaussian" << std::endl;
        return 1;
    }
    std::cout << "gen_gaussian_bits(0x9e3779b97f4a7c15) = 0x" << std::hex
              << rng.gen_gaussian_bits(0x9E3779B97F4A7C15ull) << std::dec
              << " (matches gen_gaussian)" << std::endl;

    return 0;
}
