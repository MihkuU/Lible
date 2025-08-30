cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Setting up some external libraries
# List of libraries:
#   - armadillo, https://arma.sourceforge.net/

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)
execute_process(COMMAND "${CMAKE_COMMAND}" . WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/)
execute_process(COMMAND "${CMAKE_COMMAND}" --build . WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/)

# arma
target_compile_definitions(lible PRIVATE ARMA_DONT_USE_WRAPPER)
target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/include)
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib")
    target_link_libraries(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib/libarmadillo.a)
elseif(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib64")
    target_link_libraries(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib64/libarmadillo.a)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_definitions(lible PRIVATE ARMA_NO_DEBUG)
endif()
