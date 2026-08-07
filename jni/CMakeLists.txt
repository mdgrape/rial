cmake_minimum_required(VERSION 3.10)
project(rialcpp VERSION 0.1.0)

find_package(JNI REQUIRED)

include(cmake/CPM.cmake)
include(CTest) # enable_testing() must run in the top level to reach `test/`

set(RIAL_JAR_PATH "${CMAKE_SOURCE_DIR}/../target/scala-2.13/rial.jar")

# ------------------------------------------------------------------------------

add_subdirectory(src)
add_subdirectory(test)
