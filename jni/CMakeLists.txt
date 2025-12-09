cmake_minimum_required(VERSION 3.10)
project(rialcpp VERSION 0.1.0)

find_package(JNI REQUIRED)

set(RIAL_JAR_PATH "${CMAKE_SOURCE_DIR}/../target/scala-2.12/rial.jar")

# ------------------------------------------------------------------------------

add_subdirectory(src)
