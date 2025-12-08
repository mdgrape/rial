add_library(crial_table OBJECT chebyshev.c gentable.c gen_funcev.c nonbond_functions.c)
target_include_directories(crial_table PUBLIC
    "${PROJECT_SOURCE_DIR}/include")

# install(FILES chebyshev.h gentable.h gen_funcev.h nonbond_functions.h
#   DESTINATION include/mdg)
