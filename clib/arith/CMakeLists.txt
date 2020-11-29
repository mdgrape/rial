cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(arith_objs OBJECT fpmd.c mt.c random.c scale.c mdg_threefry.c)

target_include_directories(arith_objs PRIVATE . ../include ../table)

