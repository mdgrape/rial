add_library(rialcpp SHARED
    JVMManager.cpp real.cpp math.cpp rng.cpp
)

target_compile_features(rialcpp PRIVATE cxx_std_20)
target_compile_options(rialcpp PRIVATE -Wall -Wextra -Wpedantic)

target_link_libraries(rialcpp PRIVATE ${JNI_LIBRARIES} ${CMAKE_DL_LIBS})
target_include_directories(rialcpp
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${JNI_INCLUDE_DIRS}
)

# rialcpp finds rial.jar beside itself at run time, at
# `<dir of librialcpp>/../share/rialcpp/rial.jar` (see resolve_jar_path in
# JVMManager.cpp).  Staging the jar at that spot in the build tree as well means
# the build tree runs without being installed.
add_custom_command(TARGET rialcpp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_BINARY_DIR}/share/rialcpp"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${RIAL_JAR_PATH}" "${CMAKE_BINARY_DIR}/share/rialcpp/rial.jar"
    COMMENT "staging rial.jar for the build tree"
    VERBATIM
)

get_filename_component(JVM_LIB_DIR "${JAVA_JVM_LIBRARY}" DIRECTORY)
set_target_properties(rialcpp PROPERTIES
    BUILD_RPATH   "${JVM_LIB_DIR}"
    INSTALL_RPATH "${JVM_LIB_DIR}"
)
message(STATUS "JVM_LIBRARY = ${JAVA_JVM_LIBRARY}")
message(STATUS "JVM_LIB_DIR = ${JVM_LIB_DIR}")

# -----------------------------------------------------------------------------

install(TARGETS rialcpp
    EXPORT rialcppTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/rialcpp DESTINATION include)
install(FILES ${RIAL_JAR_PATH} DESTINATION share/rialcpp)
install(EXPORT rialcppTargets
    FILE        rialcppTargets.cmake
    NAMESPACE   rialcpp::
    DESTINATION lib/cmake/rialcpp
)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/rialcppConfig.cmake.in
    ${PROJECT_BINARY_DIR}/rialcppConfig.cmake
    INSTALL_DESTINATION lib/cmake/rialcpp
)
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/rialcppConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${PROJECT_BINARY_DIR}/rialcppConfig.cmake
    ${PROJECT_BINARY_DIR}/rialcppConfigVersion.cmake
    DESTINATION lib/cmake/rialcpp
)
