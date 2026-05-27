#ifndef RIALCPP_JVM_MANAGER_HPP
#define RIALCPP_JVM_MANAGER_HPP

#include <jni.h>

#include <array>
#include <atomic>
#include <mutex>
#include <stdexcept>
#include <string>
#include <optional>

#include <cassert>

namespace rial
{
namespace detail
{

class JVMManager
{
  public:

    static JVMManager& get_instance()
    {
        static JVMManager instance;
        return instance;
    }

  private:

    JVMManager();
    ~JVMManager();

    std::optional<std::string> initialize();
    void shutdown();

  public:

    JVMManager(const JVMManager&) = delete;
    JVMManager(JVMManager&&)      = delete;
    JVMManager& operator=(const JVMManager&) = delete;
    JVMManager& operator=(JVMManager&&)      = delete;


    // ------------------------------------------------------------------------
    // supported funcs

    float call_sqrt       (float x)          const { return checked_call(FuncSqrt,       x); }
    float call_invsqrt    (float x)          const { return checked_call(FuncInvsqrt,    x); }
    float call_reciprocal (float x)          const { return checked_call(FuncReciprocal, x); }
    float call_sin        (float x)          const { return checked_call(FuncSin,        x); }
    float call_cos        (float x)          const { return checked_call(FuncCos,        x); }
    float call_acos       (float x)          const { return checked_call(FuncAcos,       x); }
    float call_atan2      (float y, float x) const { return checked_call(FuncAtan2,   y, x); }
    float call_exp        (float x)          const { return checked_call(FuncExp,        x); }
    float call_log        (float x)          const { return checked_call(FuncLog,        x); }

    // RNG: argument is the raw random bits (uint64), result is the generated float.
    float call_gen_float12 (uint64_t rnd) const { return checked_call(FuncGenFloat12,  static_cast<jlong>(rnd)); }
    float call_gen_float01 (uint64_t rnd) const { return checked_call(FuncGenFloat01,  static_cast<jlong>(rnd)); }
    float call_gen_gaussian(uint64_t rnd) const { return checked_call(FuncGenGaussian, static_cast<jlong>(rnd)); }

  private:

    template<typename ... Args>
    float checked_call(std::size_t i, Args ... args) const
    {
        assert(i < NFunc);

        using namespace std::literals::string_literals;
        const std::array<std::string, NFunc> names{{
            "sqrt"s,
            "invsqrt"s,
            "reciprocal"s,
            "sin"s,
            "cos"s,
            "acos"s,
            "atan2"s,
            "exp"s,
            "log"s,
            "gen_float12"s,
            "gen_float01"s,
            "gen_gaussian"s
        }};

        const auto [obj, method] = funcs_.at(i);
        if(obj == nullptr || method == nullptr)
        {
            throw std::runtime_error("JVMManager: "s + names.at(i) +
                    " is not found in "s + std::string(RIAL_JAR_PATH));
        }
        return this->env_->CallFloatMethod(obj, method, args...);
    }

  private:

    JavaVM* jvm_;
    JNIEnv* env_;

    // list of currently supported funcs (to cache java ref)
    constexpr static std::size_t FuncSqrt       = 0;
    constexpr static std::size_t FuncInvsqrt    = 1;
    constexpr static std::size_t FuncReciprocal = 2;
    constexpr static std::size_t FuncSin        = 3;
    constexpr static std::size_t FuncCos        = 4;
    constexpr static std::size_t FuncAcos       = 5;
    constexpr static std::size_t FuncAtan2      = 6;
    constexpr static std::size_t FuncExp        = 7;
    constexpr static std::size_t FuncLog        = 8;
    constexpr static std::size_t FuncGenFloat12  = 9;
    constexpr static std::size_t FuncGenFloat01  = 10;
    constexpr static std::size_t FuncGenGaussian = 11;
    constexpr static std::size_t NFunc          = 12; // size
    std::array<std::pair<jobject, jmethodID>, NFunc> funcs_;

    std::mutex mutex_;
    std::atomic<bool> initialized_{false};
};

} // detail
} // rial
#endif // RIALCPP_JVM_MANAGER_HPP
