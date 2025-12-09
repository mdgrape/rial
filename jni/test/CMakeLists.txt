cmake_minimum_required(VERSION 3.10)
project(rialcpp-test)

find_package(rialcpp REQUIRED)

add_executable(main main.cpp)
target_compile_features(main PRIVATE cxx_std_20)
target_compile_options (main PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries  (main PRIVATE rialcpp::rialcpp)
