CPMAddPackage(
    NAME doctest
    GITHUB_REPOSITORY "doctest/doctest"
    VERSION 2.5.3
)

# ------------------------------------------------------------------------------
# Reference implementation of Threefry

include(FetchContent)
FetchContent_Declare(
    random123
    GIT_REPOSITORY https://github.com/DEShawResearch/random123.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(random123)
FetchContent_GetProperties(random123)

add_library(random123_lib INTERFACE)
target_include_directories(random123_lib INTERFACE
    ${random123_SOURCE_DIR}/include
)

# ------------------------------------------------------------------------------

add_executable(main main.cpp)
target_compile_features(main PRIVATE cxx_std_20)
target_compile_options (main PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries  (main PRIVATE rialcpp)

# ------------------------------------------------------------------------------

# The entry point, shared by every test binary. It carries doctest's
# implementation and `main`, so each executable below only compiles its own
# tests. See doctest_main.cpp for why it is ours rather than
# doctest::doctest_with_main.
add_library(test_main STATIC doctest_main.cpp)
target_compile_features(test_main PUBLIC cxx_std_20)
target_link_libraries  (test_main PUBLIC doctest::doctest)

include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)

# One binary per source file: each gets its own JVM, so a failure in one cannot
# take the others with it.
foreach(name bits math rng threefry)
    add_executable(test_${name} test_${name}.cpp)
    target_compile_options (test_${name} PRIVATE -Wall -Wextra -Wpedantic)
    target_link_libraries  (test_${name} PRIVATE rialcpp test_main)
    doctest_discover_tests(test_${name})
endforeach()

# the reference implementation, used by test_threefry only
target_link_libraries(test_threefry PRIVATE random123_lib)
