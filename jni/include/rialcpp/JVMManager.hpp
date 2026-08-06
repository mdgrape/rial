#ifndef RIALCPP_JVM_MANAGER_HPP
#define RIALCPP_JVM_MANAGER_HPP

#include <jni.h>

#include <atomic>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstdint>

namespace rial
{
namespace detail
{

/** True until the JVM is shut down.
 *
 *  Constant-initialized and never destroyed, so an object holding a global
 *  reference can still ask whether releasing it is safe even if it happens to
 *  outlive the JVMManager singleton during static destruction.
 */
bool jvm_is_alive() noexcept;

/** Owns the JVM and the low-level JNI plumbing.
 *
 *  It knows nothing about the simulated functions themselves: the facades
 *  (rial::MathSim, rial::RngSim, ...) instantiate their own Scala-side object
 *  through `create_object` / `singleton_object` and cache the method ids they
 *  need.  That way one process can hold several instances configured
 *  differently.
 */
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

    JNIEnv* env() const noexcept {return this->env_;}

    /** Instantiates `class_name` (an internal name such as "rial/math/MathSim")
     *  through the constructor taking `args.size()` ints, and returns a global
     *  reference to it.  Release it with `release_object`. */
    jobject create_object(const char* class_name, const std::vector<jint>& args) const;

    /** Same, for a constructor with a mixed signature.  `args` must match
     *  `ctor_sig`; local references inside it stay owned by the caller. */
    jobject create_object(const char* class_name, const char* ctor_sig,
                          const std::vector<jvalue>& args) const;

    /** Builds a java String[]; the returned local reference is the caller's. */
    jobjectArray new_string_array(const std::vector<std::string>& strs) const;

    /** Copies a java String[] out; `arr` stays the caller's. */
    std::vector<std::string> to_string_vector(jobjectArray arr) const;

    /** Returns a global reference to the singleton instance of the Scala
     *  `object class_name` (its `MODULE$` field). */
    jobject singleton_object(const char* class_name) const;

    /** Resolves a method of the runtime class of `obj`. */
    jmethodID method_of(jobject obj, const char* name, const char* sig) const;

    void release_object(jobject obj) const noexcept;

    /** Turns a pending Java exception into a std::runtime_error.  Every call
     *  into Scala goes through this: a misconfigured instance reports the Scala
     *  side's message (bad table widths, unsupported spec, ...) instead of
     *  silently returning garbage. */
    void check_exception(const char* context) const;

  private:

    jclass find_class(const std::string& name) const;

  private:

    JavaVM* jvm_;
    JNIEnv* env_;

    std::mutex mutex_;
    std::atomic<bool> initialized_{false};
};

} // detail
} // rial
#endif // RIALCPP_JVM_MANAGER_HPP
