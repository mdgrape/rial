#include "rialcpp/JVMManager.hpp"
#include "rialcpp/rng.hpp"

namespace rial
{

// The Scala counterpart is `class rial.rng.RngSim`, whose constructor takes the
// floating point spec and the generator parameters as plain ints (see
// src/main/scala/rng/RngSim.scala).
struct RngSim::Impl
{
    Impl(const RealSpec& spec, const RngConfig& config)
        : obj_(nullptr)
    {
        const auto& jvm = detail::JVMManager::get_instance();

        this->obj_ = jvm.create_object("rial/rng/RngSim", {
                spec.exW,            spec.exBias,          spec.manW,
                config.float01W,     config.gaussianW,     config.gaussianOrder,
                config.gaussianLumpOctaves, config.irwinHallW
            });

        this->gen_float12_  = jvm.method_of(this->obj_, "genFloat12",  "(J)F");
        this->gen_float01_  = jvm.method_of(this->obj_, "genFloat01",  "(J)F");
        this->gen_gaussian_ = jvm.method_of(this->obj_, "genGaussian", "(J)F");
        this->irwin_hall_   = jvm.method_of(this->obj_, "irwinHall",  "([J)F");

        this->gen_float12_bits_  = jvm.method_of(this->obj_, "genFloat12Bits",  "(J)J");
        this->gen_float01_bits_  = jvm.method_of(this->obj_, "genFloat01Bits",  "(J)J");
        this->gen_gaussian_bits_ = jvm.method_of(this->obj_, "genGaussianBits", "(J)J");
        this->irwin_hall_bits_   = jvm.method_of(this->obj_, "irwinHallBits",  "([J)J");
    }
    ~Impl()
    {
        if(this->obj_ != nullptr && detail::jvm_is_alive())
        {
            detail::JVMManager::get_instance().release_object(this->obj_);
        }
    }

    Impl(const Impl&)            = delete;
    Impl& operator=(const Impl&) = delete;
    Impl(Impl&&)                 = delete;
    Impl& operator=(Impl&&)      = delete;

    float call(jmethodID method, const char* name, std::uint64_t rnd) const
    {
        const auto& jvm = detail::JVMManager::get_instance();
        const float result = jvm.env()->CallFloatMethod(
                this->obj_, method, static_cast<jlong>(rnd));
        jvm.check_exception(name);
        return result;
    }

    std::uint64_t call_bits(jmethodID method, const char* name, std::uint64_t rnd) const
    {
        const auto& jvm = detail::JVMManager::get_instance();
        const jlong result = jvm.env()->CallLongMethod(
                this->obj_, method, static_cast<jlong>(rnd));
        jvm.check_exception(name);
        return static_cast<std::uint64_t>(result);
    }

    /** builds the long[] Irwin-Hall takes; the local reference is the
     *  caller's. */
    jlongArray new_rnd_array(const std::array<std::uint64_t, 12>& rnds,
                             const char* name) const
    {
        const auto& jvm = detail::JVMManager::get_instance();
        JNIEnv* env = jvm.env();

        const jsize      n   = static_cast<jsize>(rnds.size());
        const jlongArray arr = env->NewLongArray(n);
        if(arr == nullptr)
        {
            jvm.check_exception(name);
            throw std::runtime_error(std::string("rial: ") + name +
                    ": NewLongArray failed");
        }
        env->SetLongArrayRegion(arr, 0, n,
                reinterpret_cast<const jlong*>(rnds.data()));
        return arr;
    }

    jobject   obj_;
    jmethodID gen_float12_;
    jmethodID gen_float01_;
    jmethodID gen_gaussian_;
    jmethodID irwin_hall_;
    jmethodID gen_float12_bits_;
    jmethodID gen_float01_bits_;
    jmethodID gen_gaussian_bits_;
    jmethodID irwin_hall_bits_;
};

RngSim::RngSim(const RealSpec& spec, const RngConfig& config)
    : spec_(spec), config_(config),
      impl_(std::make_unique<Impl>(spec, config))
{}

RngSim::~RngSim() = default;
RngSim::RngSim(RngSim&&) noexcept = default;
RngSim& RngSim::operator=(RngSim&&) noexcept = default;

float RngSim::gen_float12 (std::uint64_t rnd) const {return impl_->call(impl_->gen_float12_,  "RngSim::gen_float12",  rnd);}
float RngSim::gen_float01 (std::uint64_t rnd) const {return impl_->call(impl_->gen_float01_,  "RngSim::gen_float01",  rnd);}
float RngSim::gen_gaussian(std::uint64_t rnd) const {return impl_->call(impl_->gen_gaussian_, "RngSim::gen_gaussian", rnd);}

float RngSim::irwin_hall(const std::array<std::uint64_t, 12>& rnds) const
{
    const auto& jvm = detail::JVMManager::get_instance();
    JNIEnv* env = jvm.env();

    const jlongArray arr = impl_->new_rnd_array(rnds, "RngSim::irwin_hall");
    const float result = env->CallFloatMethod(this->impl_->obj_,
            this->impl_->irwin_hall_, arr);
    env->DeleteLocalRef(arr);
    jvm.check_exception("RngSim::irwin_hall");
    return result;
}

std::uint64_t RngSim::gen_float12_bits (std::uint64_t rnd) const {return impl_->call_bits(impl_->gen_float12_bits_,  "RngSim::gen_float12_bits",  rnd);}
std::uint64_t RngSim::gen_float01_bits (std::uint64_t rnd) const {return impl_->call_bits(impl_->gen_float01_bits_,  "RngSim::gen_float01_bits",  rnd);}
std::uint64_t RngSim::gen_gaussian_bits(std::uint64_t rnd) const {return impl_->call_bits(impl_->gen_gaussian_bits_, "RngSim::gen_gaussian_bits", rnd);}

std::uint64_t RngSim::irwin_hall_bits(const std::array<std::uint64_t, 12>& rnds) const
{
    const auto& jvm = detail::JVMManager::get_instance();
    JNIEnv* env = jvm.env();

    const jlongArray arr = impl_->new_rnd_array(rnds, "RngSim::irwin_hall_bits");
    const jlong result = env->CallLongMethod(this->impl_->obj_,
            this->impl_->irwin_hall_bits_, arr);
    env->DeleteLocalRef(arr);
    jvm.check_exception("RngSim::irwin_hall_bits");
    return static_cast<std::uint64_t>(result);
}

// -----------------------------------------------------------------------------
// The raw bit generators live in the stateless scala `object rial.rng.BitRngSim`
// -- no floating point format, no table, so one process-wide handle is enough.

namespace
{

struct BitRngHandle
{
    BitRngHandle()
    {
        const auto& jvm = detail::JVMManager::get_instance();
        this->obj_ = jvm.singleton_object("rial/rng/BitRngSim");
        this->squares32_    = jvm.method_of(this->obj_, "squares32",       "(JJ)I");
        this->threefry4x32_ = jvm.method_of(this->obj_, "threefry4x32", "([I[II)[I");
    }
    // the reference is intentionally never released: it is process-wide and
    // outlives every caller.

    jobject   obj_;
    jmethodID squares32_;
    jmethodID threefry4x32_;
};

const BitRngHandle& bitrng()
{
    static const BitRngHandle handle;
    return handle;
}

} // namespace

std::uint32_t squares32(std::uint64_t count, std::uint64_t key)
{
    const auto& jvm = detail::JVMManager::get_instance();
    const auto& brg = bitrng();

    const jint result = jvm.env()->CallIntMethod(brg.obj_, brg.squares32_,
            static_cast<jlong>(count), static_cast<jlong>(key));
    jvm.check_exception("squares32");
    return static_cast<std::uint32_t>(result);
}

std::array<std::uint32_t, 4> threefry4x32(
        const std::array<std::uint32_t, 4>& ctr,
        const std::array<std::uint32_t, 4>& key,
        int rounds)
{
    const auto& jvm = detail::JVMManager::get_instance();
    const auto& brg = bitrng();
    JNIEnv* env = jvm.env();

    const jintArray jctr = env->NewIntArray(4);
    const jintArray jkey = env->NewIntArray(4);
    if(jctr == nullptr || jkey == nullptr)
    {
        jvm.check_exception("threefry4x32");
        throw std::runtime_error("rial: threefry4x32: NewIntArray failed");
    }
    env->SetIntArrayRegion(jctr, 0, 4, reinterpret_cast<const jint*>(ctr.data()));
    env->SetIntArrayRegion(jkey, 0, 4, reinterpret_cast<const jint*>(key.data()));

    const jintArray result = static_cast<jintArray>(env->CallObjectMethod(
            brg.obj_, brg.threefry4x32_, jctr, jkey, static_cast<jint>(rounds)));

    env->DeleteLocalRef(jkey);
    env->DeleteLocalRef(jctr);
    jvm.check_exception("threefry4x32");

    std::array<std::uint32_t, 4> r{};
    env->GetIntArrayRegion(result, 0, 4, reinterpret_cast<jint*>(r.data()));
    env->DeleteLocalRef(result);
    return r;
}

} // rial
