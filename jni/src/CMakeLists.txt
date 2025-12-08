add_library(rialcpp SHARED
    JVMManager.cpp
)

target_compile_features(rialcpp PRIVATE cxx_std_20)
target_compile_options(rialcpp PRIVATE -Wall -Wextra -Wpedantic)

target_link_libraries(rialcpp PRIVATE ${JNI_LIBRARIES})
target_include_directories(rialcpp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${JNI_INCLUDE_DIRS}
)
target_compile_definitions(rialcpp PRIVATE
    RIAL_JAR_PATH="${CMAKE_INSTALL_PREFIX}/share/rialcpp/rial.jar"
)

install(TARGETS rialcpp
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORIES ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
install(FILES       ${RIAL_JAR_PATH}             DESTINATION share/rialcpp)
