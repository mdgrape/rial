#ifndef RIAL_JNI_JVM_MANAGER_HPP
#define RIAL_JNI_JVM_MANAGER_HPP

#include <jni.h>

#include <atomic>
#include <mutex>
#include <string>
#include <optional>
#include <unordered_map>

#ifndef RIAL_JAR_PATH
#error "RIAL_JAR_PATH must be defined"
#endif

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

    JVMManager()
        : jvm_(nullptr), env_(nullptr)
    {}

    ~JVMManager()
    {
        if(initialized_.load())
        {
            cleanup();
        }
    }

    void cleanup();

  public:

    JVMManager(const JVMManager&) = delete;
    JVMManager(JVMManager&&)      = delete;
    JVMManager& operator=(const JVMManager&) = delete;
    JVMManager& operator=(JVMManager&&)      = delete;

    std::optional<std::string> initialize();
    void shutdown();

    float call_sqrt(float x);

  private:

    JavaVM* jvm_;
    JNIEnv* env_;

    std::unordered_map<std::string, std::pair<jobject, jmethodID>> funcs_;

    std::mutex mutex_;
    std::atomic<bool> initialized_{false};
};

} // detail
} // rial
#endif // MATHFUNCTION_H
