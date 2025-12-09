add_library(crial_arith OBJECT fpmd.c mt.c random.c scale.c mdg_threefry.c)
target_link_libraries(crial_arith PUBLIC random123_lib)
target_include_directories(crial_arith PUBLIC "${PROJECT_SOURCE_DIR}/include")
