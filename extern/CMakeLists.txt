cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Setting up some external libraries
# List of libraries:
#   - armadillo, https://arma.sourceforge.net/
#   - fmt, https://github.com/fmtlib/fmt
#   - json, https://github.com/nlohmann/json

configure_file(${CMAKE_SOURCE_DIR}/extern/CMakeLists.txt.in ${CMAKE_BINARY_DIR}/extern/CMakeLists.txt @ONLY)
execute_process(COMMAND "${CMAKE_COMMAND}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/extern)
execute_process(COMMAND "${CMAKE_COMMAND}" --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/extern)

target_include_directories(lible PUBLIC ${CMAKE_BINARY_DIR}/extern/arma/build/include)
target_link_libraries(lible ${CMAKE_BINARY_DIR}/extern/arma/build/lib/libarmadillo.so)
target_include_directories(lible PUBLIC ${CMAKE_BINARY_DIR}/extern/fmt/include)
target_include_directories(lible PUBLIC ${CMAKE_BINARY_DIR}/extern/json/include)