#include "rialcpp/JVMManager.hpp"

#include <dlfcn.h>

#include <array>
#include <string>

#include <cstdlib>

namespace rial
{
namespace detail
{

namespace
{
// Constant-initialized, so it is alive for the whole program -- including the
// static destruction phase, where objects holding global references may still
// be running their destructors after JVMManager is gone.
std::atomic<bool> g_jvm_alive{false};

/** Where rial.jar is, worked out at run time rather than baked in.
 *
 *  The jar and the library are one unit -- the scala classes have to be the
 *  ones this build was written against -- so the library looks for the jar
 *  beside itself:
 *
 *      <directory holding librialcpp>/../share/rialcpp/rial.jar
 *
 *  which is where `cmake --install` puts it, and where the build tree stages it
 *  too, so nothing has to be installed to run out of a build directory. Since
 *  the path is relative to the library rather than absolute, several installs
 *  can sit side by side and each one finds its own jar; nothing is tied to the
 *  prefix the library happened to be configured with.
 *
 *  `RIALCPP_JAR` overrides it, for pointing a build at a jar somewhere else.
 */
std::string resolve_jar_path()
{
    if (const char* env = std::getenv("RIALCPP_JAR"); env != nullptr && *env != '\0')
    {
        return std::string(env);
    }

    // any address inside this library will do; dladdr reports the file it was
    // loaded from.
    Dl_info info;
    if (::dladdr(reinterpret_cast<const void*>(&resolve_jar_path), &info) != 0 &&
        info.dli_fname != nullptr)
    {
        std::string self(info.dli_fname);

        // it may be relative to the working directory, which the caller is free
        // to change later on
        if (char* abs = ::realpath(self.c_str(), nullptr); abs != nullptr)
        {
            self.assign(abs);
            std::free(abs);
        }

        // .../lib/librialcpp.so -> .../lib -> ...
        const auto lib = self.find_last_of('/');
        if (lib != std::string::npos)
        {
            const auto prefix = self.find_last_of('/', lib - 1);
            if (prefix != std::string::npos)
            {
                return self.substr(0, prefix) + "/share/rialcpp/rial.jar";
            }
        }
    }

    return "rial.jar"; // nothing better to go on; let the JVM try the cwd
}

const std::string& jar_path()
{
    static const std::string path = resolve_jar_path();
    return path;
}
} // namespace

bool jvm_is_alive() noexcept
{
    return g_jvm_alive.load();
}

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

    std::string jarpath("-Djava.class.path=" + jar_path());
    options.at(0).optionString = jarpath.data();

    vm_args.version  = JNI_VERSION_1_8;
    vm_args.nOptions = 1;
    vm_args.options  = options.data();
    vm_args.ignoreUnrecognized = JNI_FALSE;

    jint rc = JNI_CreateJavaVM(&jvm_, (void**)&env_, &vm_args);
    if(rc == JNI_OK)
    {
        initialized_.store(true);
        g_jvm_alive.store(true);
    }
    else
    {
        initialized_.store(false);
        return std::string("Failed to create JVM (error code: ") +
            std::to_string(rc) + std::string(")");
    }
    return std::nullopt;
}

void JVMManager::shutdown()
{
    std::lock_guard<std::mutex> lock(mutex_);

    if ( ! this->initialized_.load())
    {
        return;
    }

    g_jvm_alive.store(false);
    if (jvm_)
    {
        jvm_->DestroyJavaVM();
        jvm_ = nullptr;
        env_ = nullptr;
    }
    initialized_.store(false);
}

// ----------------------------------------------------------------------------

void JVMManager::check_exception(const char* context) const
{
    using namespace std::literals::string_literals;

    if(this->env_->ExceptionCheck() == JNI_FALSE)
    {
        return;
    }

    const jthrowable exc = this->env_->ExceptionOccurred();
    this->env_->ExceptionClear();

    // Throwable.toString() gives "java.lang.IllegalArgumentException: <msg>",
    // which is what we want to surface as the C++ error message.
    std::string msg;
    const jclass    cls = this->env_->FindClass("java/lang/Throwable");
    const jmethodID mid = (cls == nullptr) ? nullptr :
        this->env_->GetMethodID(cls, "toString", "()Ljava/lang/String;");
    if(exc != nullptr && mid != nullptr)
    {
        const jstring str = static_cast<jstring>(
                this->env_->CallObjectMethod(exc, mid));
        if(this->env_->ExceptionCheck()) {this->env_->ExceptionClear();}
        if(str != nullptr)
        {
            const char* chars = this->env_->GetStringUTFChars(str, nullptr);
            if(chars != nullptr)
            {
                msg = chars;
                this->env_->ReleaseStringUTFChars(str, chars);
            }
            this->env_->DeleteLocalRef(str);
        }
    }
    if(cls != nullptr) {this->env_->DeleteLocalRef(cls);}
    if(exc != nullptr) {this->env_->DeleteLocalRef(exc);}

    if(msg.empty()) {msg = "(exception without a message)"s;}

    throw std::runtime_error("rial: "s + context + ": "s + msg);
}

jclass JVMManager::find_class(const std::string& name) const
{
    using namespace std::literals::string_literals;

    const jclass cls = this->env_->FindClass(name.c_str());
    if(cls == nullptr)
    {
        this->env_->ExceptionDescribe();
        this->env_->ExceptionClear();
        throw std::runtime_error("rial: class "s + name + " is not found in "s +
                jar_path());
    }
    return cls;
}

jobject JVMManager::create_object(const char* class_name,
        const std::vector<jint>& args) const
{
    using namespace std::literals::string_literals;

    std::vector<jvalue> argv(args.size());
    for(std::size_t i=0; i<args.size(); ++i) {argv.at(i).i = args.at(i);}

    const auto sig = "("s + std::string(args.size(), 'I') + ")V"s;
    return this->create_object(class_name, sig.c_str(), argv);
}

jobject JVMManager::create_object(const char* class_name, const char* ctor_sig,
        const std::vector<jvalue>& argv) const
{
    using namespace std::literals::string_literals;

    const jclass cls = this->find_class(class_name);

    const std::string sig(ctor_sig);
    const jmethodID ctor = this->env_->GetMethodID(cls, "<init>", sig.c_str());
    if(ctor == nullptr)
    {
        this->env_->ExceptionClear();
        this->env_->DeleteLocalRef(cls);
        throw std::runtime_error("rial: "s + class_name +
                " has no constructor "s + sig + " in "s +
                jar_path());
    }

    const jobject local = this->env_->NewObjectA(cls, ctor, argv.data());
    this->check_exception((std::string(class_name) + " constructor"s).c_str());
    if(local == nullptr)
    {
        this->env_->DeleteLocalRef(cls);
        throw std::runtime_error("rial: failed to construct "s + class_name);
    }

    const jobject global = this->env_->NewGlobalRef(local);
    this->env_->DeleteLocalRef(local);
    this->env_->DeleteLocalRef(cls);
    if(global == nullptr)
    {
        throw std::runtime_error("rial: NewGlobalRef failed for "s + class_name);
    }
    return global;
}

jobject JVMManager::singleton_object(const char* class_name) const
{
    using namespace std::literals::string_literals;

    // a scala `object Foo` compiles to the class Foo$ holding its instance in
    // the static field MODULE$.
    const auto objpath  = class_name + "$"s;
    const auto Lobjpath = "L"s + objpath + ";"s;

    const jclass cls = this->find_class(objpath);

    const jfieldID module = this->env_->GetStaticFieldID(cls, "MODULE$", Lobjpath.c_str());
    if(module == nullptr)
    {
        this->env_->ExceptionClear();
        this->env_->DeleteLocalRef(cls);
        throw std::runtime_error("rial: "s + objpath +
                " has no MODULE$ field in "s + jar_path());
    }

    const jobject local = this->env_->GetStaticObjectField(cls, module);
    this->check_exception((objpath + ".MODULE$"s).c_str());
    if(local == nullptr)
    {
        this->env_->DeleteLocalRef(cls);
        throw std::runtime_error("rial: "s + objpath + ".MODULE$ is null"s);
    }

    const jobject global = this->env_->NewGlobalRef(local);
    this->env_->DeleteLocalRef(local);
    this->env_->DeleteLocalRef(cls);
    if(global == nullptr)
    {
        throw std::runtime_error("rial: NewGlobalRef failed for "s + objpath);
    }
    return global;
}

jobjectArray JVMManager::new_string_array(const std::vector<std::string>& strs) const
{
    using namespace std::literals::string_literals;

    const jclass cls = this->find_class("java/lang/String");

    const jobjectArray arr = this->env_->NewObjectArray(
            static_cast<jsize>(strs.size()), cls, nullptr);
    this->env_->DeleteLocalRef(cls);
    if(arr == nullptr)
    {
        this->check_exception("new_string_array");
        throw std::runtime_error("rial: NewObjectArray(String) failed"s);
    }

    for(std::size_t i=0; i<strs.size(); ++i)
    {
        const jstring str = this->env_->NewStringUTF(strs.at(i).c_str());
        if(str == nullptr)
        {
            this->env_->DeleteLocalRef(arr);
            this->check_exception("new_string_array");
            throw std::runtime_error("rial: NewStringUTF failed for "s + strs.at(i));
        }
        this->env_->SetObjectArrayElement(arr, static_cast<jsize>(i), str);
        this->env_->DeleteLocalRef(str);
    }
    return arr;
}

jintArray JVMManager::new_int_array(const std::vector<jint>& values) const
{
    using namespace std::literals::string_literals;

    const jintArray arr =
        this->env_->NewIntArray(static_cast<jsize>(values.size()));
    if(arr == nullptr)
    {
        this->check_exception("new_int_array");
        throw std::runtime_error("rial: NewIntArray failed"s);
    }
    this->env_->SetIntArrayRegion(arr, 0, static_cast<jsize>(values.size()),
            values.data());
    return arr;
}

std::vector<std::string> JVMManager::to_string_vector(jobjectArray arr) const
{
    std::vector<std::string> strs;
    if(arr == nullptr) {return strs;}

    const jsize n = this->env_->GetArrayLength(arr);
    strs.reserve(static_cast<std::size_t>(n));
    for(jsize i=0; i<n; ++i)
    {
        const jstring str = static_cast<jstring>(
                this->env_->GetObjectArrayElement(arr, i));
        if(str == nullptr) {continue;}

        const char* chars = this->env_->GetStringUTFChars(str, nullptr);
        if(chars != nullptr)
        {
            strs.emplace_back(chars);
            this->env_->ReleaseStringUTFChars(str, chars);
        }
        this->env_->DeleteLocalRef(str);
    }
    return strs;
}

jmethodID JVMManager::method_of(jobject obj, const char* name, const char* sig) const
{
    using namespace std::literals::string_literals;

    const jclass cls = this->env_->GetObjectClass(obj);
    if(cls == nullptr)
    {
        this->env_->ExceptionClear();
        throw std::runtime_error("rial: GetObjectClass failed for "s + name);
    }

    const jmethodID method = this->env_->GetMethodID(cls, name, sig);
    this->env_->DeleteLocalRef(cls);
    if(method == nullptr)
    {
        this->env_->ExceptionClear();
        throw std::runtime_error("rial: method "s + name + " : "s + sig +
                " is not found in "s + jar_path());
    }
    return method;
}

void JVMManager::release_object(jobject obj) const noexcept
{
    if(obj != nullptr && this->env_ != nullptr)
    {
        this->env_->DeleteGlobalRef(obj);
    }
}

} // detail
} // rial
