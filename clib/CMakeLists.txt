cmake_minimum_required(VERSION 3.10)
project(crial)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

# ------------------------------------------------------------------------------

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

add_subdirectory(arith)
add_subdirectory(table)

add_library(crial SHARED)
target_link_libraries(crial PUBLIC random123_lib crial_arith crial_table)
target_include_directories(crial PUBLIC "${PROJECT_SOURCE_DIR}/include")

INSTALL(TARGETS crial
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/crial)
