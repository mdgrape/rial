#ifndef RIALCPP_RNG_HPP
#define RIALCPP_RNG_HPP

#include <rialcpp/config.hpp>

#include <array>
#include <memory>

#include <cstdint>

namespace rial
{

/** Software emulation of the float-producing random number generators.
 *
 *  Each of them turns raw random bits into a float, so the result depends on
 *  the floating point format and on how many bits are consumed: as with
 *  rial::MathSim, one instance carries one configuration.
 *
 *  {{{
 *  const rial::RngSim rng; // the standard FP32 setup
 *  const float u = rng.gen_float01(rial::squares32(0, key));
 *  }}}
 *
 *  Not thread-safe.  Throws std::runtime_error on JVM / Scala-side errors.
 */
class RngSim
{
  public:

    /** The standard FP32 setup, i.e. what the fixed free functions used to
     *  provide. */
    RngSim() : RngSim(RealSpec::float32(), RngConfig::standard()) {}

    RngSim(const RealSpec& spec, const RngConfig& config);
    ~RngSim();

    RngSim(RngSim&&) noexcept;
    RngSim& operator=(RngSim&&) noexcept;
    RngSim(const RngSim&)            = delete;
    RngSim& operator=(const RngSim&) = delete;

    /** uniform in [1, 2); the low `manW` bits of `rnd` become the mantissa. */
    float gen_float12 (std::uint64_t rnd) const;
    /** uniform in [0, 1) from the low `float01W` bits of `rnd`. */
    float gen_float01 (std::uint64_t rnd) const;
    /** N(0, 1) from the low `gaussianW` bits of `rnd`. */
    float gen_gaussian(std::uint64_t rnd) const;

    /** approx-gaussian via Irwin-Hall: 12 uniform words (the low `irwinHallW`
     *  bits of each are used) summed and converted to a float. */
    float irwin_hall(const std::array<std::uint64_t, 12>& rnds) const;

    const RealSpec&  spec  () const noexcept {return this->spec_;}
    const RngConfig& config() const noexcept {return this->config_;}

  private:

    struct Impl;

    RealSpec              spec_;
    RngConfig             config_;
    std::unique_ptr<Impl> impl_;
};

// -----------------------------------------------------------------------------
// Counter-based raw bit generators.  They produce integers and never touch the
// floating point format, so they need no configuration object and stay free
// functions.

/** Squares32 counter-based uniform RNG (Widynski 2022): the 32-bit random word
 *  for (count, key).  The hardware is a fixed 4-round pipeline, so there is no
 *  round count to choose. */
std::uint32_t squares32(std::uint64_t count, std::uint64_t key);

/** Threefry-4x32 (Salmon et al. 2011): the four 32-bit random words for
 *  (ctr, key).
 *
 *  @param rounds number of rounds; 20 is the standard Threefry-4x32-20.
 */
std::array<std::uint32_t, 4> threefry4x32(
        const std::array<std::uint32_t, 4>& ctr,
        const std::array<std::uint32_t, 4>& key,
        int rounds = 20);

} // rial
#endif//RIALCPP_RNG_HPP
