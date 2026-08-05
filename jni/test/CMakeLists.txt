cmake_minimum_required(VERSION 3.14)
project(rialcpp-test)

find_package(rialcpp REQUIRED)

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
target_link_libraries  (main PRIVATE rialcpp::rialcpp random123_lib)
