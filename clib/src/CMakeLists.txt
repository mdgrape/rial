add_subdirectory(arith)
add_subdirectory(table)

add_library(crial SHARED)
target_link_libraries(crial PUBLIC random123_lib crial_arith crial_table)
target_include_directories(crial PUBLIC "${PROJECT_SOURCE_DIR}/include")

INSTALL(TARGETS crial
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/crial)
