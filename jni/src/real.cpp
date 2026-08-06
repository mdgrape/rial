#include "rialcpp/JVMManager.hpp"
#include "rialcpp/real.hpp"

namespace rial
{

// The conversion lives in the stateless scala `object rial.arith.RealConvSim`
// -- the spec goes with each call, so one process-wide handle is enough.
namespace
{

struct RealConvHandle
{
    RealConvHandle()
    {
        const auto& jvm = detail::JVMManager::get_instance();
        this->obj_ = jvm.singleton_object("rial/arith/RealConvSim");
        this->to_double_   = jvm.method_of(this->obj_, "toDouble",   "(IIIJ)D");
        this->to_float_    = jvm.method_of(this->obj_, "toFloat",    "(IIIJ)F");
        this->from_double_ = jvm.method_of(this->obj_, "fromDouble", "(IIID)J");
        this->from_float_  = jvm.method_of(this->obj_, "fromFloat",  "(IIIF)J");
    }
    // the reference is intentionally never released: it is process-wide and
    // outlives every caller.

    jobject   obj_;
    jmethodID to_double_;
    jmethodID to_float_;
    jmethodID from_double_;
    jmethodID from_float_;
};

const RealConvHandle& conv()
{
    static const RealConvHandle handle;
    return handle;
}

} // namespace

double to_double(const RealSpec& spec, std::uint64_t bits)
{
    const auto& jvm = detail::JVMManager::get_instance();
    const auto& cnv = conv();

    const double result = jvm.env()->CallDoubleMethod(cnv.obj_, cnv.to_double_,
            spec.exW, spec.exBias, spec.manW, static_cast<jlong>(bits));
    jvm.check_exception("to_double");
    return result;
}

float to_float(const RealSpec& spec, std::uint64_t bits)
{
    const auto& jvm = detail::JVMManager::get_instance();
    const auto& cnv = conv();

    const float result = jvm.env()->CallFloatMethod(cnv.obj_, cnv.to_float_,
            spec.exW, spec.exBias, spec.manW, static_cast<jlong>(bits));
    jvm.check_exception("to_float");
    return result;
}

std::uint64_t to_bits(const RealSpec& spec, double x)
{
    const auto& jvm = detail::JVMManager::get_instance();
    const auto& cnv = conv();

    const jlong result = jvm.env()->CallLongMethod(cnv.obj_, cnv.from_double_,
            spec.exW, spec.exBias, spec.manW, x);
    jvm.check_exception("to_bits");
    return static_cast<std::uint64_t>(result);
}

std::uint64_t to_bits(const RealSpec& spec, float x)
{
    const auto& jvm = detail::JVMManager::get_instance();
    const auto& cnv = conv();

    const jlong result = jvm.env()->CallLongMethod(cnv.obj_, cnv.from_float_,
            spec.exW, spec.exBias, spec.manW, x);
    jvm.check_exception("to_bits");
    return static_cast<std::uint64_t>(result);
}

} // rial
