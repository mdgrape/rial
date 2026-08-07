// Loose comparison of MathSim against libm. (strict tests are done in Scala)

#include <rialcpp/config.hpp>
#include <rialcpp/math.hpp>

#include <doctest/doctest.h>

#include <stdexcept>
#include <string>
#include <vector>

#include <cmath>
#include <cstddef>

double smg_ref(double x)
{
    const double sgm_a = std::exp(-1.0);
    const double sgm_b = std::exp(-6.0);
    const double a2 = sgm_a * sgm_a;
    const double b2 = sgm_b * sgm_b;
    const double p2 = (a2 * b2) / (a2 - b2);
    const double g  = std::exp(-x * x / (2.0 * p2));
    return -x / a2 * ((a2 / p2) / (sgm_b / (sgm_a * g) + 1.0) + 1.0);
}

/** One unary function: how to call it, what it should be, and where to look. */
struct UnaryFunc
{
    rial::Func  func;
    std::string name;
    /** what the function pulls in as a dependency, itself included. */
    rial::Func  effective;

    float  (rial::MathSim::*sim)(float) const;
    double (*ref)(double);

    double from;
    double to;
};

std::vector<UnaryFunc> unary_funcs()
{
    constexpr double pi = 3.141592653589793;

    return {
        {rial::Func::sqrt, "sqrt", rial::Func::sqrt,
         &rial::MathSim::sqrt, [](double x){return std::sqrt(x);}, 0.0, 100.0},

        {rial::Func::invsqrt, "invsqrt", rial::Func::invsqrt,
         &rial::MathSim::invsqrt, [](double x){return 1.0 / std::sqrt(x);}, 0.001, 100.0},

        {rial::Func::reciprocal, "reciprocal", rial::Func::reciprocal,
         &rial::MathSim::reciprocal, [](double x){return 1.0 / x;}, 0.001, 100.0},

        // sin and cos share one polynomial table, so either pulls in the other.
        {rial::Func::sin, "sin", rial::Func::sin | rial::Func::cos,
         &rial::MathSim::sin, [](double x){return std::sin(x);}, -pi, pi},

        {rial::Func::cos, "cos", rial::Func::sin | rial::Func::cos,
         &rial::MathSim::cos, [](double x){return std::cos(x);}, -pi, pi},

        // the first phase of acos reuses the sqrt table.
        {rial::Func::acos, "acos", rial::Func::acos | rial::Func::sqrt,
         &rial::MathSim::acos, [](double x){return std::acos(x);}, -1.0, 1.0},

        {rial::Func::exp, "exp", rial::Func::exp,
         &rial::MathSim::exp, [](double x){return std::exp(x);}, -10.0, 10.0},

        {rial::Func::log, "log", rial::Func::log,
         &rial::MathSim::log, [](double x){return std::log(x);}, 0.001, 100.0},

        // past |x| ~ 6 sigmoid rounds to 1 and the Scala side says so on stdout.
        {rial::Func::sigmoid, "sigmoid", rial::Func::sigmoid,
         &rial::MathSim::sigmoid, [](double x){return 1.0 / (1.0 + std::exp(-x));}, -6.0, 6.0},

        {rial::Func::softplus, "softplus", rial::Func::softplus,
         &rial::MathSim::softplus, [](double x){return std::log1p(std::exp(x));}, -10.0, 10.0},

        {rial::Func::scale_mixture_gaussian, "scale_mixture_gaussian",
         rial::Func::scale_mixture_gaussian,
         &rial::MathSim::scale_mixture_gaussian, smg_ref, -1.0, 1.0},
    };
}

/** Error relative to the magnitude of the result, but never dividing by less
 *  than one: an absolute check where the result is small, a relative one where
 *  it is large. */
double error(double ref, double got)
{
    return std::abs(ref - got) / std::max(1.0, std::abs(ref));
}

/** The worst point of a sweep, kept so a failure can name it. */
struct Worst
{
    double err = -1.0;
    double x   = 0.0;
    double y   = 0.0; // second argument, for atan2
    double ref = 0.0;
    double got = 0.0;
};

/** Sweeps `n` points over [from, to] and returns the worst disagreement. */
template<typename Sim, typename Ref>
Worst sweep(Sim sim, Ref ref, double from, double to, std::size_t n)
{
    Worst worst;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double x = from + (to - from) * i / (n - 1);
        const double r = ref(x);
        if (!std::isfinite(r)) {continue;} // the endpoint of log, invsqrt, ...

        const double g = (sim)(static_cast<float>(x));
        const double e = error(r, g);
        if (worst.err < e) {worst = Worst{e, x, 0.0, r, g};}
    }
    return worst;
}

/** The single-function config to test `f` with.
 *
 *  scale_mixture_gaussian needs the wider polynomial output of `fp32_all`; with
 *  the standard one its result collapses (see MathFuncConfig::fp32_all). */
rial::MathFuncConfig single_func_config(rial::Func f)
{
    rial::MathFuncConfig config = (f == rial::Func::scale_mixture_gaussian)
        ? rial::MathFuncConfig::fp32_all()
        : rial::MathFuncConfig::fp32();
    config.funcs = f;
    return config;
}

TEST_CASE("MathSim agrees with libm to within 1e-3 in FP32")
{
    // Func::all so that every function is available from one instance; the
    // per-function configs are exercised by the config test cases below.
    const rial::MathSim sim(rial::RealSpec::float32(),
                            rial::MathFuncConfig::fp32_all());

    constexpr double tol = 1.0e-3;

    for (const auto& f : unary_funcs())
    {
        INFO("func = ", f.name);
        REQUIRE(sim.supports(f.func));

        const auto worst = sweep([&](float x){return (sim.*(f.sim))(x);},
                                 f.ref, f.from, f.to, 1000);

        INFO("worst at x = ", worst.x,
             ": libm = ", worst.ref, ", rial = ", worst.got);
        CHECK(worst.err <= tol);
    }

    SUBCASE("atan2")
    {
        Worst worst;
        for (std::size_t i = 0; i < 64; ++i)
        {
        for (std::size_t j = 0; j < 64; ++j)
        {
            const double x = -1.0 + 2.0 * i / 63;
            const double y = -1.0 + 2.0 * j / 63;
            if (x == 0.0 && y == 0.0) {continue;}

            const double r = std::atan2(y, x);
            const double g = sim.atan2(static_cast<float>(y), static_cast<float>(x));
            const double e = error(r, g);
            if (worst.err < e) {worst = Worst{e, x, y, r, g};}
        }
        }
        INFO("worst at (x, y) = (", worst.x, ", ", worst.y,
             "): libm = ", worst.ref, ", rial = ", worst.got);
        CHECK(worst.err <= tol);
    }
}

TEST_CASE("MathSim in BF16 answers in BF16, and still with the right function")
{
    // The same functions in a 7 bit mantissa and a plain lookup table. That the
    // RealSpec crossed JNI is checked directly -- every result has to be a BF16
    // value, so it round-trips through the format unchanged -- and the loose
    // libm comparison then says it is still the function it should be. The
    // tolerance is the format's and the table's, not the approximation's.
    const rial::MathSim sim(rial::RealSpec::bfloat16(),
                            rial::MathFuncConfig::bf16_all());

    constexpr double tol = 0.25;

    for (const auto& f : unary_funcs())
    {
        INFO("func = ", f.name);
        REQUIRE(sim.supports(f.func));

        for (std::size_t i = 0; i <= 8; ++i)
        {
            const float x = static_cast<float>(f.from + (f.to - f.from) * i / 8);
            const float y = (sim.*(f.sim))(x);
            INFO("at x = ", x, ", rial = ", y);
            CHECK(sim.to_float(sim.to_bits(y)) == y);
        }

        const auto worst = sweep([&](float x){return (sim.*(f.sim))(x);},
                                 f.ref, f.from, f.to, 1000);

        INFO("worst at x = ", worst.x,
             ": libm = ", worst.ref, ", rial = ", worst.got);
        CHECK(worst.err <= tol);
    }
}

TEST_CASE("a config gives exactly the requested functions and their dependencies")
{
    for (const auto& f : unary_funcs())
    {
        INFO("func = ", f.name);
        const rial::MathSim sim(rial::RealSpec::float32(), single_func_config(f.func));

        // the requested function, plus whatever it reuses the table of
        CHECK(sim.config().funcs == f.effective);
        CHECK(sim.supports(f.func));

        // and nothing else: every other function is rejected rather than
        // silently answered by whichever table happens to be there
        for (const auto& g : unary_funcs())
        {
            if (has(f.effective, g.func)) {continue;}
            INFO("unsupported = ", g.name);
            CHECK_FALSE(sim.supports(g.func));
            CHECK_THROWS_AS((sim.*(g.sim))(0.5f), std::runtime_error);
        }
    }

    SUBCASE("atan2 pulls in reciprocal")
    {
        const rial::MathSim sim(rial::RealSpec::float32(),
                                rial::MathFuncConfig::fp32(rial::Func::atan2));
        CHECK(sim.config().funcs == (rial::Func::atan2 | rial::Func::reciprocal));
        CHECK(sim.supports(rial::Func::reciprocal));
        CHECK_FALSE(sim.supports(rial::Func::sqrt));
        CHECK_THROWS_AS(sim.sqrt(0.5f), std::runtime_error);
    }
}

TEST_CASE("a single-function config still computes that function")
{
    // The point of the per-function configs: selecting one function must give
    // that function's table, not the first one in the list.
    constexpr double tol = 1.0e-3;

    for (const auto& f : unary_funcs())
    {
        INFO("func = ", f.name);
        const rial::MathSim sim(rial::RealSpec::float32(), single_func_config(f.func));

        const auto worst = sweep([&](float x){return (sim.*(f.sim))(x);},
                                 f.ref, f.from, f.to, 100);

        INFO("worst at x = ", worst.x,
             ": libm = ", worst.ref, ", rial = ", worst.got);
        CHECK(worst.err <= tol);
    }
}
