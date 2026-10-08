#include "rialcpp/JVMManager.hpp"
#include "rialcpp/math.hpp"

#include <string>
#include <utility>
#include <vector>

namespace rial
{

namespace
{

// The names `rial.math.MathSim` selects functions by (MathSim.funcNames).
constexpr std::pair<Func, const char*> func_names[] = {
    {Func::sqrt,                   "sqrt"},
    {Func::invsqrt,                "invsqrt"},
    {Func::reciprocal,             "reciprocal"},
    {Func::sin,                    "sin"},
    {Func::cos,                    "cos"},
    {Func::acos,                   "acos"},
    {Func::atan2,                  "atan2"},
    {Func::exp,                    "exp"},
    {Func::log,                    "log"},
    {Func::sigmoid,                "sigmoid"},
    {Func::softplus,               "softplus"},
    {Func::scale_mixture_gaussian, "scaleMixtureGaussian"},
};

// `funcs` and the tables they read, the way `rial.math.FuncKind.normalize`
// adds them: acos reads the sqrt table and atan2 the reciprocal one, so naming
// either names that table too and hands it the width it is configured at.
constexpr Func with_dependencies(Func funcs) noexcept
{
    if(has(funcs, Func::acos))  {funcs |= Func::sqrt;}
    if(has(funcs, Func::atan2)) {funcs |= Func::reciprocal;}
    return funcs;
}

std::vector<std::string> to_names(Func funcs)
{
    std::vector<std::string> names;
    for(const auto& [flag, name] : func_names)
    {
        if(has(funcs, flag)) {names.emplace_back(name);}
    }
    if(names.empty())
    {
        throw std::runtime_error(
            "rial: MathSim: at least one function must be selected");
    }
    return names;
}

// The width each of those tables is read at (MathSim's `addrWs`), in the same
// order; a function the configuration gives no width of its own is read at
// `MathFuncConfig::addrW`.
std::vector<jint> to_addr_widths(const MathFuncConfig& config, Func funcs)
{
    std::vector<jint> widths;
    for(const auto& fn : func_names)
    {
        if(has(funcs, fn.first)) {widths.push_back(config.addrW_of(fn.first));}
    }
    return widths;
}

Func to_funcs(const std::vector<std::string>& names)
{
    Func funcs = Func::none;
    for(const auto& name : names)
    {
        for(const auto& [flag, known] : func_names)
        {
            if(name == known) {funcs |= flag; break;}
        }
    }
    return funcs;
}

const char* name_of(Func f) noexcept
{
    for(const auto& [flag, name] : func_names)
    {
        if(flag == f) {return name;}
    }
    return "(unknown)";
}

} // namespace

// The Scala counterpart is `class rial.math.MathSim`, whose constructor takes
// the floating point spec, the table parameters (the address width of each
// table among them) and the function names (see
// src/main/scala/math/MathSim.scala). One Impl == one Scala-side instance,
// holding the generated tables.
struct MathSim::Impl
{
    Impl(const RealSpec& spec, const MathFuncConfig& config)
        : obj_(nullptr)
    {
        const auto& jvm = detail::JVMManager::get_instance();
        JNIEnv* env = jvm.env();

        const Func funcs = with_dependencies(config.funcs);
        const jobjectArray names = jvm.new_string_array(to_names(funcs));
        jintArray addrWs = nullptr;

        try
        {
            addrWs = jvm.new_int_array(to_addr_widths(config, funcs));

            std::vector<jvalue> argv(11);
            argv.at(0).i  = spec.exW;
            argv.at(1).i  = spec.exBias;
            argv.at(2).i  = spec.manW;
            argv.at(3).i  = config.nOrder;
            argv.at(4).i  = config.addrW;
            argv.at(5).l  = addrWs;
            argv.at(6).i  = config.inBitW;
            argv.at(7).i  = config.outBitW;
            argv.at(8).l  = names;
            argv.at(9).d  = config.smgSigmaA;
            argv.at(10).d = config.smgSigmaB;

            this->obj_ = jvm.create_object("rial/math/MathSim",
                    "(IIIII[III[Ljava/lang/String;DD)V", argv);
        }
        catch(...)
        {
            if(addrWs != nullptr) {env->DeleteLocalRef(addrWs);}
            env->DeleteLocalRef(names);
            throw;
        }
        env->DeleteLocalRef(addrWs);
        env->DeleteLocalRef(names);

        this->sqrt_       = jvm.method_of(this->obj_, "sqrt",       "(F)F");
        this->invsqrt_    = jvm.method_of(this->obj_, "invsqrt",    "(F)F");
        this->reciprocal_ = jvm.method_of(this->obj_, "reciprocal", "(F)F");
        this->sin_        = jvm.method_of(this->obj_, "sin",        "(F)F");
        this->cos_        = jvm.method_of(this->obj_, "cos",        "(F)F");
        this->acos_       = jvm.method_of(this->obj_, "acos",       "(F)F");
        this->atan2_      = jvm.method_of(this->obj_, "atan2",     "(FF)F");
        this->exp_        = jvm.method_of(this->obj_, "exp",        "(F)F");
        this->log_        = jvm.method_of(this->obj_, "log",        "(F)F");
        this->sigmoid_    = jvm.method_of(this->obj_, "sigmoid",    "(F)F");
        this->softplus_   = jvm.method_of(this->obj_, "softplus",   "(F)F");
        this->smg_        = jvm.method_of(this->obj_, "scaleMixtureGaussian", "(F)F");

        this->sqrt_bits_       = jvm.method_of(this->obj_, "sqrtBits",       "(J)J");
        this->invsqrt_bits_    = jvm.method_of(this->obj_, "invsqrtBits",    "(J)J");
        this->reciprocal_bits_ = jvm.method_of(this->obj_, "reciprocalBits", "(J)J");
        this->sin_bits_        = jvm.method_of(this->obj_, "sinBits",        "(J)J");
        this->cos_bits_        = jvm.method_of(this->obj_, "cosBits",        "(J)J");
        this->acos_bits_       = jvm.method_of(this->obj_, "acosBits",       "(J)J");
        this->atan2_bits_      = jvm.method_of(this->obj_, "atan2Bits",     "(JJ)J");
        this->exp_bits_        = jvm.method_of(this->obj_, "expBits",        "(J)J");
        this->log_bits_        = jvm.method_of(this->obj_, "logBits",        "(J)J");
        this->sigmoid_bits_    = jvm.method_of(this->obj_, "sigmoidBits",    "(J)J");
        this->softplus_bits_   = jvm.method_of(this->obj_, "softplusBits",   "(J)J");
        this->smg_bits_        = jvm.method_of(this->obj_, "scaleMixtureGaussianBits", "(J)J");

        this->supported_  = jvm.method_of(this->obj_, "supportedFuncs",
                                          "()[Ljava/lang/String;");
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

    /** The function set the Scala side ended up with: what was requested plus
     *  the tables its dependencies pulled in. */
    Func supported_funcs() const
    {
        const auto& jvm = detail::JVMManager::get_instance();
        JNIEnv* env = jvm.env();

        const jobjectArray names = static_cast<jobjectArray>(
                env->CallObjectMethod(this->obj_, this->supported_));
        jvm.check_exception("MathSim::supportedFuncs");

        const Func funcs = to_funcs(jvm.to_string_vector(names));
        env->DeleteLocalRef(names);
        return funcs;
    }

    template<typename ... Args>
    float call(jmethodID method, const char* name, Args ... args) const
    {
        const auto& jvm = detail::JVMManager::get_instance();
        const float result = jvm.env()->CallFloatMethod(this->obj_, method, args...);
        jvm.check_exception(name);
        return result;
    }

    template<typename ... Args>
    std::uint64_t call_bits(jmethodID method, const char* name, Args ... args) const
    {
        const auto& jvm = detail::JVMManager::get_instance();
        const jlong result = jvm.env()->CallLongMethod(this->obj_, method,
                static_cast<jlong>(args)...);
        jvm.check_exception(name);
        return static_cast<std::uint64_t>(result);
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
    jmethodID sigmoid_;
    jmethodID softplus_;
    jmethodID smg_;

    jmethodID sqrt_bits_;
    jmethodID invsqrt_bits_;
    jmethodID reciprocal_bits_;
    jmethodID sin_bits_;
    jmethodID cos_bits_;
    jmethodID acos_bits_;
    jmethodID atan2_bits_;
    jmethodID exp_bits_;
    jmethodID log_bits_;
    jmethodID sigmoid_bits_;
    jmethodID softplus_bits_;
    jmethodID smg_bits_;

    jmethodID supported_;
};

MathSim::MathSim(const RealSpec& spec, const MathFuncConfig& config)
    : spec_(spec), config_(config),
      impl_(std::make_unique<Impl>(spec, config))
{
    // report what is really available, dependencies included
    this->config_.funcs = this->impl_->supported_funcs();
}

MathSim::~MathSim() = default;
MathSim::MathSim(MathSim&&) noexcept = default;
MathSim& MathSim::operator=(MathSim&&) noexcept = default;

void MathSim::require(Func f) const
{
    if( ! this->supports(f))
    {
        throw std::runtime_error(std::string("rial: MathSim: ") + name_of(f) +
                " is not supported by this configuration");
    }
}

float MathSim::sqrt      (float x) const {require(Func::sqrt);       return impl_->call(impl_->sqrt_,       "MathSim::sqrt",       x);}
float MathSim::invsqrt   (float x) const {require(Func::invsqrt);    return impl_->call(impl_->invsqrt_,    "MathSim::invsqrt",    x);}
float MathSim::reciprocal(float x) const {require(Func::reciprocal); return impl_->call(impl_->reciprocal_, "MathSim::reciprocal", x);}
float MathSim::sin       (float x) const {require(Func::sin);        return impl_->call(impl_->sin_,        "MathSim::sin",        x);}
float MathSim::cos       (float x) const {require(Func::cos);        return impl_->call(impl_->cos_,        "MathSim::cos",        x);}
float MathSim::acos      (float x) const {require(Func::acos);       return impl_->call(impl_->acos_,       "MathSim::acos",       x);}
float MathSim::exp       (float x) const {require(Func::exp);        return impl_->call(impl_->exp_,        "MathSim::exp",        x);}
float MathSim::log       (float x) const {require(Func::log);        return impl_->call(impl_->log_,        "MathSim::log",        x);}
float MathSim::sigmoid   (float x) const {require(Func::sigmoid);    return impl_->call(impl_->sigmoid_,    "MathSim::sigmoid",    x);}
float MathSim::softplus  (float x) const {require(Func::softplus);   return impl_->call(impl_->softplus_,   "MathSim::softplus",   x);}

float MathSim::scale_mixture_gaussian(float x) const
{
    require(Func::scale_mixture_gaussian);
    return impl_->call(impl_->smg_, "MathSim::scale_mixture_gaussian", x);
}

float MathSim::atan2(float y, float x) const
{
    require(Func::atan2);
    return impl_->call(impl_->atan2_, "MathSim::atan2", y, x);
}

// -----------------------------------------------------------------------------
// bit-pattern interface

std::uint64_t MathSim::sqrt_bits      (std::uint64_t x) const {require(Func::sqrt);       return impl_->call_bits(impl_->sqrt_bits_,       "MathSim::sqrt_bits",       x);}
std::uint64_t MathSim::invsqrt_bits   (std::uint64_t x) const {require(Func::invsqrt);    return impl_->call_bits(impl_->invsqrt_bits_,    "MathSim::invsqrt_bits",    x);}
std::uint64_t MathSim::reciprocal_bits(std::uint64_t x) const {require(Func::reciprocal); return impl_->call_bits(impl_->reciprocal_bits_, "MathSim::reciprocal_bits", x);}
std::uint64_t MathSim::sin_bits       (std::uint64_t x) const {require(Func::sin);        return impl_->call_bits(impl_->sin_bits_,        "MathSim::sin_bits",        x);}
std::uint64_t MathSim::cos_bits       (std::uint64_t x) const {require(Func::cos);        return impl_->call_bits(impl_->cos_bits_,        "MathSim::cos_bits",        x);}
std::uint64_t MathSim::acos_bits      (std::uint64_t x) const {require(Func::acos);       return impl_->call_bits(impl_->acos_bits_,       "MathSim::acos_bits",       x);}
std::uint64_t MathSim::exp_bits       (std::uint64_t x) const {require(Func::exp);        return impl_->call_bits(impl_->exp_bits_,        "MathSim::exp_bits",        x);}
std::uint64_t MathSim::log_bits       (std::uint64_t x) const {require(Func::log);        return impl_->call_bits(impl_->log_bits_,        "MathSim::log_bits",        x);}
std::uint64_t MathSim::sigmoid_bits   (std::uint64_t x) const {require(Func::sigmoid);    return impl_->call_bits(impl_->sigmoid_bits_,    "MathSim::sigmoid_bits",    x);}
std::uint64_t MathSim::softplus_bits  (std::uint64_t x) const {require(Func::softplus);   return impl_->call_bits(impl_->softplus_bits_,   "MathSim::softplus_bits",   x);}

std::uint64_t MathSim::scale_mixture_gaussian_bits(std::uint64_t x) const
{
    require(Func::scale_mixture_gaussian);
    return impl_->call_bits(impl_->smg_bits_, "MathSim::scale_mixture_gaussian_bits", x);
}

std::uint64_t MathSim::atan2_bits(std::uint64_t y, std::uint64_t x) const
{
    require(Func::atan2);
    return impl_->call_bits(impl_->atan2_bits_, "MathSim::atan2_bits", y, x);
}

} // rial
