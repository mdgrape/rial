#include "rial/JVMManager.hpp"

#include <array>
#include <string>
#include <vector>

namespace rial
{
namespace detail
{

std::optional<std::string> JVMManager::initialize()
{
    std::lock_guard<std::mutex> lock(this->mutex_); // lock!

    if (initialized_.load())
    {
        return std::nullopt; // already ok
    }

    JavaVMInitArgs vm_args;
    std::array<JavaVMOption, 1> options;

    const std::string jarpath("-Djava.class.path=" RIAL_JAR_PATH);
    std::vector<char> jarpath_raw(jarpath.size()+1, '\0');
    for(std::size_t i=0; i<jarpath.size(); ++i)
    {
        jarpath_raw.at(i) = jarpath.at(i);
    }
    options.at(0).optionString = jarpath_raw.data();

    vm_args.version  = JNI_VERSION_1_8;
    vm_args.nOptions = 1;
    vm_args.options  = options.data();
    vm_args.ignoreUnrecognized = false;

    jint rc = JNI_CreateJavaVM(&jvm_, (void**)&env_, &vm_args);
    if(rc == JNI_OK)
    {
        initialized_.store(true);
        return std::nullopt;
    }
    else
    {
        initialized_.store(false);
        return std::string("Failed to create JVM (error code: ") +
            std::to_string(rc) + std::string(")");
    }
}

void JVMManager::shutdown()
{
    std::lock_guard<std::mutex> lock(mutex_);

    if ( ! this->initialized_.load())
    {
        return;
    }
    this->cleanup();
    initialized_.store(false);
}

void JVMManager::cleanup()
{
    if (env_)
    {
        for(auto& [k, v] : funcs_)
        {
            auto& [instance, method] = v;
            env_->DeleteGlobalRef(instance);

            instance = nullptr;
            method   = nullptr;
        }
    }
    if (jvm_)
    {
        jvm_->DestroyJavaVM();
        jvm_ = nullptr;
        env_ = nullptr;
    }
}

} // detail
} // rial
