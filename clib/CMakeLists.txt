cmake_minimum_required(VERSION 3.10)
project(crial)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

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
message("include dir: ${random123_SOURCE_DIR}/include")

add_subdirectory(arith)
add_subdirectory(table)

add_library(crial SHARED
  $<TARGET_OBJECTS:arith_objs>
  $<TARGET_OBJECTS:table_objs> )
target_link_libraries(crial PUBLIC random123_lib)

set_target_properties(crial PROPERTIES PUBLIC_HEADER
  "include/defs.h;include/md.h;include/mdg4a_atom.h;include/r2sqrt.h;include/rinv.h;include/rinv3.h")

INSTALL(TARGETS crial
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/crial)


