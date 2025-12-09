#include "rialcpp/JVMManager.hpp"

#ifndef RIAL_JAR_PATH
#error "RIAL_JAR_PATH must be defined"
#endif

#include <array>
#include <string>
#include <vector>

namespace rial
{
namespace detail
{

JVMManager::JVMManager()
    : jvm_(nullptr), env_(nullptr)
{
    const auto err = this->initialize();
    if(err.has_value())
    {
        this->shutdown();
        throw std::runtime_error(err.value());
    }
}
JVMManager::~JVMManager()
{
    this->shutdown();
}

std::optional<std::string> JVMManager::initialize()
{
    std::lock_guard<std::mutex> lock(this->mutex_);

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
    }
    else
    {
        initialized_.store(false);
        return std::string("Failed to create JVM (error code: ") +
            std::to_string(rc) + std::string(")");
    }

    // -----------------------------------------------------------------------
    // find funcs from library

    const auto find_func = [this](const std::string& objname,
        const std::string& funcname, const std::string& funcsig,
        std::pair<jobject, jmethodID>& result
    ) -> std::optional<std::string> {
        using namespace std::literals::string_literals;
        const auto objpath  = "rial/math/"s + objname + "$"s;
        const auto Lobjpath = "L"s + objpath + ";"s;

        const auto obj_class = this->env_->FindClass(objpath.c_str());
        if(obj_class == nullptr)
        {
            this->env_->ExceptionClear();
            return "FindClass failed with: "s + objpath;
        }

        const auto module_field = this->env_->GetStaticFieldID(obj_class, "MODULE$", Lobjpath.c_str());
        if(module_field == nullptr)
        {
            this->env_->ExceptionClear();
            return "GetStaticFieldID failed with: "s + Lobjpath;
        }

        const auto instance = this->env_->GetStaticObjectField(obj_class, module_field);
        if(instance == nullptr)
        {
            this->env_->ExceptionClear();
            return "GetStaticObjectField failed with: "s + objname;
        }

        const auto func_method = this->env_->GetMethodID(obj_class, funcname.c_str(), funcsig.c_str());
        if(func_method == nullptr)
        {
            return "GetMethodID failed with: "s + funcname + " : "s + funcsig;
        }

        result = std::make_pair(instance, func_method);
        return std::nullopt;
    };

    if(auto r = find_func("SqrtSim",        "sqrtDefaultFP32",        "(F)F",  this->funcs_.at(FuncSqrt       ))) {return r;}
    if(auto r = find_func("InvSqrtSim",     "invsqrtDefaultFP32",     "(F)F",  this->funcs_.at(FuncInvsqrt    ))) {return r;}
    if(auto r = find_func("ReciprocalSim",  "reciprocalDefaultFP32",  "(F)F",  this->funcs_.at(FuncReciprocal ))) {return r;}
    if(auto r = find_func("SinCosSim",      "sinDefaultFP32",         "(F)F",  this->funcs_.at(FuncSin        ))) {return r;}
    if(auto r = find_func("SinCosSim",      "cosDefaultFP32",         "(F)F",  this->funcs_.at(FuncCos        ))) {return r;}
    if(auto r = find_func("ACosPhase1Sim",  "acosPhase1DefaultFP32",  "(F)F",  this->funcs_.at(FuncAcosPhase1 ))) {return r;}
    if(auto r = find_func("ACosPhase2Sim",  "acosPhase2DefaultFP32",  "(F)F",  this->funcs_.at(FuncAcosPhase2 ))) {return r;}
    if(auto r = find_func("ATan2Phase1Sim", "atan2Phase1DefaultFP32", "(FF)F", this->funcs_.at(FuncAtan2Phase1))) {return r;}
    if(auto r = find_func("ATan2Phase2Sim", "atan2Phase2DefaultFP32", "(F)F",  this->funcs_.at(FuncAtan2Phase2))) {return r;}
    if(auto r = find_func("ExpSim",         "expDefaultFP32",         "(F)F",  this->funcs_.at(FuncExp        ))) {return r;}
    if(auto r = find_func("LogSim",         "logDefaultFP32",         "(F)F",  this->funcs_.at(FuncLog        ))) {return r;}

    return std::nullopt;
}

void JVMManager::shutdown()
{
    std::lock_guard<std::mutex> lock(mutex_);
    // -----------------------------------------------------------------------
    // find funcs from library

    if ( ! this->initialized_.load())
    {
        return;
    }

    if (jvm_)
    {
        jvm_->DestroyJavaVM();
        jvm_ = nullptr;
        env_ = nullptr;
    }
    initialized_.store(false);
}

} // detail
} // rial
