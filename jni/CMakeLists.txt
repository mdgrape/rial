cmake_minimum_required(VERSION 3.10)
project(rial-jni)

find_package(JNI REQUIRED)

set(RIAL_JAR_PATH "${CMAKE_SOURCE_DIR}/../target/scala-2.12/rial.jar")

# ------------------------------------------------------------------------------

add_subdirectory(src)
