#ifndef RIALCPP_MATH_HPP
#define RIALCPP_MATH_HPP

#include <rialcpp/config.hpp>

#include <memory>

namespace rial
{

/** Software emulation of the approximated math functions.
 *
 *  The result of every function depends on the floating point format and on the
 *  table/polynomial parameters, so there is no single `rial::sqrt(x)`: build a
 *  MathSim for the configuration you want and call it through that.
 *
 *  {{{
 *  const rial::MathSim fp32; // the standard FP32 setup
 *  const rial::MathSim bf16(rial::RealSpec::bfloat16(),
 *                           rial::MathFuncConfig::bf16());
 *
 *  std::cout << fp32.sqrt(2.0f) << " " << bf16.sqrt(2.0f) << std::endl;
 *  }}}
 *
 *  Constructing one generates the approximation tables on the Scala side, which
 *  is not cheap; keep the instance around instead of rebuilding it per call.
 *  Instances with different configurations may coexist, and none of them shares
 *  any state, so mixing precisions in one program is fine.
 *
 *  Not thread-safe: every call goes into the single JNI environment.
 *  Throws std::runtime_error if the JVM or the Scala side reports a problem
 *  (an invalid configuration is reported when the MathSim is constructed).
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

    const RealSpec&       spec  () const noexcept {return this->spec_;}
    const MathFuncConfig& config() const noexcept {return this->config_;}

  private:

    struct Impl;

    RealSpec              spec_;
    MathFuncConfig        config_;
    std::unique_ptr<Impl> impl_;
};

} // rial
#endif//RIALCPP_MATH_HPP
