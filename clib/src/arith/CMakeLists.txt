add_library(crial_arith OBJECT fpmd.c mt.c random.c scale.c mdg_threefry.c)
target_link_libraries(crial_arith PUBLIC random123_lib)
target_include_directories(crial_arith PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/arith"
    "${PROJECT_SOURCE_DIR}/table")
