// The entry point of the test binary, with doctest's signal handling turned
// off.
//
// A JVM raises SIGSEGV as part of running normally. Once a method is hot enough
// to be JIT compiled, the null check on a field access is left out of the
// generated code and the MMU is allowed to trap instead; the JVM's own SIGSEGV
// handler then resumes execution on the deoptimization path. Nothing is wrong
// when this happens -- it is how the null check is made free.
//
// doctest installs a SIGSEGV handler of its own around every test case, and it
// reports and re-raises rather than resuming, so that ordinary event surfaces
// as `test case CRASHED: SIGSEGV` and takes the test down. It is not reliably
// reproducible either: it needs the JIT to have compiled the code being called,
// so only the tests that go through JNI often enough hit it.
//
// The signal belongs to the JVM, so doctest is kept away from it.
#define DOCTEST_CONFIG_NO_POSIX_SIGNALS
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
