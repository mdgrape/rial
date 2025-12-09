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

    float call_sqrt       (float x)          const { return checked_call(FuncSqrt, x); }
    float call_invsqrt    (float x)          const { return checked_call(FuncInvsqrt, x); }
    float call_reciprocal (float x)          const { return checked_call(FuncReciprocal, x); }
    float call_sin        (float x)          const { return checked_call(FuncSin, x); }
    float call_cos        (float x)          const { return checked_call(FuncCos, x); }
    float call_acosPhase1 (float x)          const { return checked_call(FuncAcosPhase1, x); }
    float call_acosPhase2 (float x)          const { return checked_call(FuncAcosPhase2, x); }
    float call_acos       (float x)          const { return call_acosPhase2(call_acosPhase1(x)); }
    float call_atan2Phase1(float y, float x) const { return checked_call(FuncAtan2Phase1, y, x); }
    float call_atan2Phase2(float x)          const { return checked_call(FuncAtan2Phase2, x); }
    float call_atan2      (float y, float x) const { return call_atan2Phase2(call_atan2Phase1(y, x)); }
    float call_exp        (float x)          const { return checked_call(FuncExp, x); }
    float call_log        (float x)          const { return checked_call(FuncLog, x); }

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
            "acosPhase1"s,
            "acosPhase2"s,
            "atan2Phase1"s,
            "atan2Phase2"s,
            "exp"s,
            "log"s
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
    constexpr static std::size_t FuncSqrt        =  0;
    constexpr static std::size_t FuncInvsqrt     =  1;
    constexpr static std::size_t FuncReciprocal  =  2;
    constexpr static std::size_t FuncSin         =  3;
    constexpr static std::size_t FuncCos         =  4;
    constexpr static std::size_t FuncAcosPhase1  =  5;
    constexpr static std::size_t FuncAcosPhase2  =  6;
    constexpr static std::size_t FuncAtan2Phase1 =  7;
    constexpr static std::size_t FuncAtan2Phase2 =  8;
    constexpr static std::size_t FuncExp         =  9;
    constexpr static std::size_t FuncLog         = 10;
    constexpr static std::size_t NFunc           = 11; // size
    std::array<std::pair<jobject, jmethodID>, NFunc> funcs_;

    std::mutex mutex_;
    std::atomic<bool> initialized_{false};
};

} // detail
} // rial
#endif // RIALCPP_JVM_MANAGER_HPP
