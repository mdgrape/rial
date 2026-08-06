#include "rialcpp/JVMManager.hpp"
#include "rialcpp/math.hpp"

namespace rial
{

// The Scala counterpart is `class rial.math.MathSim`, whose constructor takes
// the floating point spec and the table parameters as plain ints (see
// src/main/scala/math/MathSim.scala). One Impl == one Scala-side instance,
// holding the generated tables.
struct MathSim::Impl
{
    Impl(const RealSpec& spec, const MathFuncConfig& config)
        : obj_(nullptr)
    {
        const auto& jvm = detail::JVMManager::get_instance();

        this->obj_ = jvm.create_object("rial/math/MathSim", {
                spec.exW,      spec.exBias,   spec.manW,
                config.nOrder, config.addrW,  config.inBitW, config.outBitW
            });

        this->sqrt_       = jvm.method_of(this->obj_, "sqrt",       "(F)F");
        this->invsqrt_    = jvm.method_of(this->obj_, "invsqrt",    "(F)F");
        this->reciprocal_ = jvm.method_of(this->obj_, "reciprocal", "(F)F");
        this->sin_        = jvm.method_of(this->obj_, "sin",        "(F)F");
        this->cos_        = jvm.method_of(this->obj_, "cos",        "(F)F");
        this->acos_       = jvm.method_of(this->obj_, "acos",       "(F)F");
        this->atan2_      = jvm.method_of(this->obj_, "atan2",     "(FF)F");
        this->exp_        = jvm.method_of(this->obj_, "exp",        "(F)F");
        this->log_        = jvm.method_of(this->obj_, "log",        "(F)F");
    }
    ~Impl()
    {
        // during static destruction the JVM may already be gone; the global
        // reference dies with it, so releasing it is neither needed nor safe.
        if(this->obj_ != nullptr && detail::jvm_is_alive())
        {
            detail::JVMManager::get_instance().release_object(this->obj_);
        }
    }

    Impl(const Impl&)            = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&)                 = delete;
    Impl& operator=(Impl&&)      = delete;

    template<typename ... Args>
    float call(jmethodID method, const char* name, Args ... args) const
    {
        const auto& jvm = detail::JVMManager::get_instance();
        const float result = jvm.env()->CallFloatMethod(this->obj_, method, args...);
        jvm.check_exception(name);
        return result;
    }

    jobject   obj_;
    jmethodID sqrt_;
    jmethodID invsqrt_;
    jmethodID reciprocal_;
    jmethodID sin_;
    jmethodID cos_;
    jmethodID acos_;
    jmethodID atan2_;
    jmethodID exp_;
    jmethodID log_;
};

MathSim::MathSim(const RealSpec& spec, const MathFuncConfig& config)
    : spec_(spec), config_(config),
      impl_(std::make_unique<Impl>(spec, config))
{}

MathSim::~MathSim() = default;
MathSim::MathSim(MathSim&&) noexcept = default;
MathSim& MathSim::operator=(MathSim&&) noexcept = default;

float MathSim::sqrt      (float x) const {return impl_->call(impl_->sqrt_,       "MathSim::sqrt",       x);}
float MathSim::invsqrt   (float x) const {return impl_->call(impl_->invsqrt_,    "MathSim::invsqrt",    x);}
float MathSim::reciprocal(float x) const {return impl_->call(impl_->reciprocal_, "MathSim::reciprocal", x);}
float MathSim::sin       (float x) const {return impl_->call(impl_->sin_,        "MathSim::sin",        x);}
float MathSim::cos       (float x) const {return impl_->call(impl_->cos_,        "MathSim::cos",        x);}
float MathSim::acos      (float x) const {return impl_->call(impl_->acos_,       "MathSim::acos",       x);}
float MathSim::exp       (float x) const {return impl_->call(impl_->exp_,        "MathSim::exp",        x);}
float MathSim::log       (float x) const {return impl_->call(impl_->log_,        "MathSim::log",        x);}

float MathSim::atan2(float y, float x) const
{
    return impl_->call(impl_->atan2_, "MathSim::atan2", y, x);
}

} // rial
