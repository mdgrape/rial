cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

add_subdirectory(arith)
add_subdirectory(table)

add_library(crial SHARED
  $<TARGET_OBJECTS:arith_objs>
  $<TARGET_OBJECTS:table_objs> )

set_target_properties(crial PROPERTIES PUBLIC_HEADER
  "include/defs.h;include/md.h;include/mdg4a_atom.h;include/r2sqrt.h;include/rinv.h;include/rinv3.h")

INSTALL(TARGETS crial
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/crial)


