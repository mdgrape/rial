cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(table_objs OBJECT chebyshev.c gentable.c gen_funcev.c nonbond_functions.c)

target_link_libraries(table_objs PUBLIC random123_lib)
target_include_directories(table_objs PRIVATE . ../include ../arith)

install(FILES chebyshev.h gentable.h gen_funcev.h nonbond_functions.h
  DESTINATION include/mdg)
