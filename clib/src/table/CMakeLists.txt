add_library(crial_table OBJECT chebyshev.c gentable.c gen_funcev.c nonbond_functions.c)
target_link_libraries(crial_table PUBLIC random123_lib)
target_include_directories(crial_table PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/arith"
    "${PROJECT_SOURCE_DIR}/table")
# install(FILES chebyshev.h gentable.h gen_funcev.h nonbond_functions.h
#   DESTINATION include/mdg)
