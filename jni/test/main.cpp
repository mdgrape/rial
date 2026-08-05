#include <rialcpp/math.hpp>
#include <rialcpp/rng.hpp>

#include <Random123/threefry.h>

// Random123 defines threefry4x32 as an object-like macro, which would rewrite
// rial::threefry4x32 too (macros do not respect namespaces). Drop it and call
// the underlying function directly. that also pins the round count here.
#undef threefry4x32

#include <array>
#include <iomanip>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include <cmath>
#include <cstdint>

// checks if threefry reference impl and scala emulation in rial matches.
std::size_t check_threefry(std::size_t n)
{
    const auto reference = [](const std::array<std::uint32_t, 4>& c,
                              const std::array<std::uint32_t, 4>& k)
    {
        threefry4x32_ctr_t ctr{{c[0], c[1], c[2], c[3]}};
        threefry4x32_key_t key{{k[0], k[1], k[2], k[3]}};
        const threefry4x32_ctr_t r = threefry4x32_R(20, ctr, key);
        return std::array<std::uint32_t, 4>{{r.v[0], r.v[1], r.v[2], r.v[3]}};
    };

    constexpr std::uint32_t lo = 0x00000000u;
    constexpr std::uint32_t hi = 0xffffffffu;
    constexpr std::uint32_t sb = 0x80000000u;

    // edge cases first, then a random sweep over the full parameter space.
    std::vector<std::pair<std::array<std::uint32_t, 4>, std::array<std::uint32_t, 4>>> cases{
        {{{lo, lo, lo, lo}}, {{lo, lo, lo, lo}}},
        {{{hi, hi, hi, hi}}, {{hi, hi, hi, hi}}},
        {{{lo, lo, lo, lo}}, {{hi, hi, hi, hi}}},
        {{{hi, hi, hi, hi}}, {{lo, lo, lo, lo}}},
        {{{sb, lo, lo, sb}}, {{0x11111111u, 0x22222222u, lo, lo}}},
        {{{lo, lo, lo, lo}}, {{0x11111111u, 0x22222222u, lo, lo}}},
        {{{1u, lo, lo, lo}}, {{0x11111111u, 0x22222222u, lo, lo}}},
        {{{hi, lo, lo, lo}}, {{0x11111111u, 0x22222222u, lo, lo}}},
    };

    std::mt19937 rng(19660809);
    std::uniform_int_distribution<std::uint32_t> word(lo, hi);
    while (cases.size() < n)
    {
        cases.push_back({
                {{word(rng), word(rng), word(rng), word(rng)}},
                {{word(rng), word(rng), word(rng), word(rng)}}
            });
    }

    std::size_t mismatch = 0;
    for (const auto& [ctr, key] : cases)
    {
        const auto rial = rial::threefry4x32(ctr, key);
        const auto ref  = reference(ctr, key);
        if (rial != ref)
        {
            if (mismatch < 10)
            {
                std::cout << std::hex
                          << "  ctr=" << ctr [0] << "," << ctr [1] << "," << ctr [2] << "," << ctr [3]
                          << " key="  << key [0] << "," << key [1] << "," << key [2] << "," << key [3]
                          << " rial=" << rial[0] << "," << rial[1] << "," << rial[2] << "," << rial[3]
                          << " ref="  << ref [0] << "," << ref [1] << "," << ref [2] << "," << ref [3]
                          << std::dec << std::endl;
            }
            ++mismatch;
        }
    }
    std::cout << "threefry4x32: " << cases.size() - mismatch << "/" << cases.size()
              << " cases match Random123" << std::endl;
    return mismatch;
}

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

    check("sqrt", 0.0, 4.0,
        [](float x) {return std::sqrt(x); },
        [](float x) {return rial::sqrt(x);}
    );

    check("invsqrt", 0.5, 4.0,
        [](float x) {return 1.0f / std::sqrt(x); },
        [](float x) {return rial::invsqrt(x);}
    );

    check("reciprocal", 0.5, 4.0,
        [](float x) {return 1.0f / x; },
        [](float x) {return rial::reciprocal(x);}
    );

    check("sin", -pi, pi,
        [](float x) {return std::sin(x); },
        [](float x) {return rial::sin(x);}
    );

    check("cos", -pi, pi,
        [](float x) {return std::cos(x); },
        [](float x) {return rial::cos(x);}
    );

    check("acos", -1.0, 1.0,
        [](float x) {return std::acos(x); },
        [](float x) {return rial::acos(x);}
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
        const auto rial = rial::atan2(y, x);
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
        [](float x) {return rial::exp(x);}
    );

    check("log", 0.1, 3.1,
        [](float x) {return std::log(x); },
        [](float x) {return rial::log(x);}
    );

    // -----------------------------------------------------------------------
    // RNG generators: input is raw random bits, output is the generated float.
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "rng" << std::endl;
    std::cout << "gen_float12(0)          = " << rial::gen_float12(0)               << " (expect 1)"    << std::endl;
    std::cout << "gen_float01(0x80000000) = " << rial::gen_float01(0x80000000ull)   << " (expect 0.5)"  << std::endl;
    std::cout << "gen_float01(0xC0000000) = " << rial::gen_float01(0xC0000000ull)   << " (expect 0.75)" << std::endl;
    for (std::uint64_t s = 1; s <= 5; ++s)
    {
        const std::uint64_t bits = s * 0x9E3779B97F4A7C15ull; // arbitrary spread
        std::cout << "gen_gaussian(0x" << std::hex << bits << std::dec << ") = "
                  << rial::gen_gaussian(bits) << std::endl;
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
    std::cout << "irwin_hall(squares32 x12) = " << rial::irwin_hall(rnds) << std::endl;

    // -----------------------------------------------------------------------
    // Threefry emulation vs the reference implementation.

    std::cout << "---------------------------------------------------" << std::endl;
    const std::size_t mismatch = check_threefry(4096);
    if (mismatch != 0)
    {
        std::cout << "FAILED: " << mismatch
                  << " threefry4x32 mismatches against Random123" << std::endl;
        return 1;
    }

    return 0;
}
