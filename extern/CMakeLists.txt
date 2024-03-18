cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Setting up some external libraries
# List of libraries:
#   - armadillo, https://arma.sourceforge.net/
#   - fmt, https://github.com/fmtlib/fmt
#   - json, https://github.com/nlohmann/json
#   - mpl, https://github.com/rabauke/mpl

configure_file(${CMAKE_SOURCE_DIR}/extern/CMakeLists.txt.in ${CMAKE_BINARY_DIR}/extern/CMakeLists.txt @ONLY)
execute_process(COMMAND "${CMAKE_COMMAND}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/extern)
execute_process(COMMAND "${CMAKE_COMMAND}" --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/extern)

target_include_directories(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/arma/build/include)
if(EXISTS "${CMAKE_BINARY_DIR}/extern/arma/build/lib")
    target_link_libraries(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/arma/build/lib/libarmadillo.a)
elseif(EXISTS "${CMAKE_BINARY_DIR}/extern/arma/build/lib64")
    target_link_libraries(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/arma/build/lib64/libarmadillo.a)
endif()


# fmt
target_compile_definitions(lible PRIVATE FMT_HEADER_ONLY)
target_include_directories(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/fmt/include)

# json
target_include_directories(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/json/include)

# TODO: figure what to do about this one
#target_include_directories(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/robinhood/src/include/)
#target_include_directories(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/)
find_package(unordered_dense CONFIG REQUIRED PATHS ${CMAKE_BINARY_DIR}/extern/unordered_dense/build/lib/cmake/unordered_dense)
target_link_libraries(lible PRIVATE unordered_dense::unordered_dense)

# mpl
if(USE_MPI)
	target_include_directories(lible PRIVATE ${CMAKE_BINARY_DIR}/extern/mpl/build/include)
endif()
