#include <rialcpp/rng.hpp>

#include <Random123/threefry.h>

// Random123 defines threefry4x32 as an object-like macro, which would rewrite
// rial::threefry4x32 too (macros do not respect namespaces). Drop it and call
// the underlying function directly. that also pins the round count here.
#undef threefry4x32

#include <doctest/doctest.h>

#include <array>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <cstdint>

using word4 = std::array<std::uint32_t, 4>;

word4 reference_threefry4x32(const word4& c, const word4& k, int rounds)
{
    threefry4x32_ctr_t ctr{{c[0], c[1], c[2], c[3]}};
    threefry4x32_key_t key{{k[0], k[1], k[2], k[3]}};
    const threefry4x32_ctr_t r = threefry4x32_R(rounds, ctr, key);
    return word4{{r.v[0], r.v[1], r.v[2], r.v[3]}};
}

std::string to_hex(const word4& w)
{
    std::ostringstream oss;
    oss << std::hex << w[0] << "," << w[1] << "," << w[2] << "," << w[3];
    return oss.str();
}

// edge cases: all-zero / all-one counters and keys, and the sign bit alone.
std::vector<std::pair<word4, word4>> edge_cases()
{
    constexpr std::uint32_t lo = 0x00000000u;
    constexpr std::uint32_t hi = 0xffffffffu;
    constexpr std::uint32_t sb = 0x80000000u;

    return {
        {{{lo, lo, lo, lo}}, {{lo, lo, lo, lo}}},
        {{{hi, hi, hi, hi}}, {{hi, hi, hi, hi}}},
        {{{lo, lo, lo, lo}}, {{hi, hi, hi, hi}}},
        {{{hi, hi, hi, hi}}, {{lo, lo, lo, lo}}},
        {{{sb, lo, lo, sb}}, {{0x11111111u, 0x22222222u, lo, lo}}},
        {{{lo, lo, lo, lo}}, {{0x11111111u, 0x22222222u, lo, lo}}},
        {{{1u, lo, lo, lo}}, {{0x11111111u, 0x22222222u, lo, lo}}},
        {{{hi, lo, lo, lo}}, {{0x11111111u, 0x22222222u, lo, lo}}},
    };
}

// a deterministic sweep over the full parameter space, `n` cases of it.
std::vector<std::pair<word4, word4>> random_cases(std::size_t n, std::uint32_t seed)
{
    std::mt19937 rng(seed);
    std::uniform_int_distribution<std::uint32_t> word(0x00000000u, 0xffffffffu);

    std::vector<std::pair<word4, word4>> cases;
    cases.reserve(n);
    while (cases.size() < n)
    {
        cases.push_back({
                {{word(rng), word(rng), word(rng), word(rng)}},
                {{word(rng), word(rng), word(rng), word(rng)}}
            });
    }
    return cases;
}

// compares the rial emulation against Random123 and reports the first
// mismatching case; returns the number of mismatches.
std::size_t compare(const std::vector<std::pair<word4, word4>>& cases, int rounds)
{
    std::size_t mismatch = 0;
    for (const auto& [ctr, key] : cases)
    {
        const auto rial = rial::threefry4x32(ctr, key, rounds);
        const auto ref  = reference_threefry4x32(ctr, key, rounds);
        if (rial != ref)
        {
            if (mismatch == 0)
            {
                INFO("rounds = ", rounds);
                INFO("ctr    = ", to_hex(ctr));
                INFO("key    = ", to_hex(key));
                INFO("rial   = ", to_hex(rial));
                INFO("ref    = ", to_hex(ref));
                CHECK(rial == ref);
            }
            ++mismatch;
        }
    }
    return mismatch;
}

TEST_CASE("threefry4x32-20 matches the Random123 reference")
{
    SUBCASE("edge cases")
    {
        CHECK(compare(edge_cases(),20) == 0);
    }

    SUBCASE("random counters and keys")
    {
        CHECK(compare(random_cases(4096, 19660809), 20) == 0);
    }
}

TEST_CASE("threefry4x32 matches the reference for other round counts")
{
    for (const int rounds : {1, 2, 3, 4, 5, 13, 32, 72})
    {
        CAPTURE(rounds);
        CHECK(compare(edge_cases(),rounds) == 0);
        CHECK(compare(random_cases(64, 19660809 + rounds), rounds) == 0);
    }
}
