#ifndef RIALCPP_MATH_HPP
#define RIALCPP_MATH_HPP

#include <rialcpp/config.hpp>
#include <rialcpp/real.hpp>

#include <memory>

#include <cstdint>

namespace rial
{

/** Software emulation of the approximated math functions.
 *
 *  The result of every function depends on the floating point format, on the
 *  table/polynomial parameters and on which functions are supported, so there
 *  is no single `rial::sqrt(x)`: build a MathSim for the configuration you want
 *  and call it through that.
 *
 *  {{{
 *  const rial::MathSim fp32; // the standard FP32 setup
 *  const rial::MathSim bf16(rial::RealSpec::bfloat16(),
 *                           rial::MathFuncConfig::bf16());
 *
 *  // a module that only does sqrt/invsqrt/reciprocal: its narrower shared
 *  // tables make it differ from `fp32` by about a ULP.
 *  const rial::MathSim simple(rial::RealSpec::float32(),
 *                             rial::MathFuncConfig::fp32(rial::Func::simple));
 *
 *  std::cout << fp32.sqrt(2.0f) << " " << bf16.sqrt(2.0f) << " "
 *            << simple.sqrt(2.0f) << std::endl;
 *  }}}
 *
 *  Constructing one generates the approximation tables on the Scala side, which
 *  is not cheap; keep the instance around instead of rebuilding it per call.
 *  Instances with different configurations may coexist, and none of them shares
 *  any state, so mixing precisions in one program is fine.
 *
 *  Not thread-safe: every call goes into the single JNI environment.
 *  Throws std::runtime_error if the JVM or the Scala side reports a problem
 *  (an invalid configuration is reported when the MathSim is constructed), or
 *  if a function outside `config().funcs` is called.
 */
class MathSim
{
  public:

    /** The standard FP32 setup (`MathFuncConfig.standardFP32` on the Scala
     *  side), i.e. what the fixed `rial::sqrt(x)` export used to provide. */
    MathSim() : MathSim(RealSpec::float32(), MathFuncConfig::fp32()) {}

    MathSim(const RealSpec& spec, const MathFuncConfig& config);
    ~MathSim();

    MathSim(MathSim&&) noexcept;
    MathSim& operator=(MathSim&&) noexcept;
    MathSim(const MathSim&)            = delete;
    MathSim& operator=(const MathSim&) = delete;

    float sqrt       (float x)          const;
    float invsqrt    (float x)          const;
    float reciprocal (float x)          const;
    float sin        (float x)          const;
    float cos        (float x)          const;
    float acos       (float x)          const;
    float atan2      (float y, float x) const;
    float exp        (float x)          const;
    float log        (float x)          const;
    float sigmoid    (float x)          const;
    float softplus   (float x)          const;
    float scale_mixture_gaussian(float x) const;

    // -------------------------------------------------------------------------
    // The same functions on the raw bit pattern of spec(), for formats a float
    // cannot carry (manW > 23, exW > 8, ...).  Bits above the spec width are
    // ignored.  Use to_double() / to_bits() to get in and out.
    //
    // {{{
    // const rial::MathSim fp64(rial::RealSpec{11, 0x3FF, 52},
    //                          rial::MathFuncConfig{rial::Func::sqrt, 2, 8, 44, 55});
    // std::cout << fp64.to_double(fp64.sqrt_bits(fp64.to_bits(2.0))) << std::endl;
    // }}}

    std::uint64_t sqrt_bits       (std::uint64_t x) const;
    std::uint64_t invsqrt_bits    (std::uint64_t x) const;
    std::uint64_t reciprocal_bits (std::uint64_t x) const;
    std::uint64_t sin_bits        (std::uint64_t x) const;
    std::uint64_t cos_bits        (std::uint64_t x) const;
    std::uint64_t acos_bits       (std::uint64_t x) const;
    std::uint64_t atan2_bits      (std::uint64_t y, std::uint64_t x) const;
    std::uint64_t exp_bits        (std::uint64_t x) const;
    std::uint64_t log_bits        (std::uint64_t x) const;
    std::uint64_t sigmoid_bits    (std::uint64_t x) const;
    std::uint64_t softplus_bits   (std::uint64_t x) const;
    std::uint64_t scale_mixture_gaussian_bits(std::uint64_t x) const;

    // conversion between spec() bit patterns and ordinary numbers
    double        to_double(std::uint64_t bits) const {return rial::to_double(spec_, bits);}
    float         to_float (std::uint64_t bits) const {return rial::to_float (spec_, bits);}
    std::uint64_t to_bits  (double x)           const {return rial::to_bits  (spec_, x);}
    std::uint64_t to_bits  (float  x)           const {return rial::to_bits  (spec_, x);}

    const RealSpec& spec() const noexcept {return this->spec_;}

    /** The configuration, with `funcs` updated to the effective function set:
     *  the dependencies a requested function pulled in are included. */
    const MathFuncConfig& config() const noexcept {return this->config_;}

    /** Shorthand for `has(config().funcs, f)`. */
    bool supports(Func f) const noexcept {return has(this->config_.funcs, f);}

  private:

    /** throws unless `f` is in the effective function set. */
    void require(Func f) const;

    struct Impl;

    RealSpec              spec_;
    MathFuncConfig        config_;
    std::unique_ptr<Impl> impl_;
};

} // rial
#endif//RIALCPP_MATH_HPP
