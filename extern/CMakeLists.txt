cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Setting up some external libraries
# List of libraries:
#   - armadillo, https://arma.sourceforge.net/
#   - ankerl, https://github.com/martinus/unordered_dense
#   - fmt, https://github.com/fmtlib/fmt
#   - json, https://github.com/nlohmann/json
#   - mpl, https://github.com/rabauke/mpl

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)
execute_process(COMMAND "${CMAKE_COMMAND}" . WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/)
execute_process(COMMAND "${CMAKE_COMMAND}" --build . WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/)

target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/include)
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib")
    target_link_libraries(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib/libarmadillo.a)
elseif(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib64")
    target_link_libraries(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib64/libarmadillo.a)
endif()

# fmt
target_compile_definitions(lible PRIVATE FMT_HEADER_ONLY)
target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/fmt/include)

# json
target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/json/include)

# ankerl-unordered dense
find_package(unordered_dense CONFIG REQUIRED PATHS ${CMAKE_CURRENT_BINARY_DIR}/unordered_dense/build/lib/cmake/unordered_dense)
target_link_libraries(lible PRIVATE unordered_dense::unordered_dense)

# mpl
if(USE_MPI)
	target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mpl/build/include)
endif()
