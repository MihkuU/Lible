cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(Lible VERSION 0.0.0 LANGUAGES CXX)

if (PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
   message(FATAL_ERROR "In-source builds are not allowed, specify another build directory")
endif()

# TODO: compile which submodule of the library.
# option(BUILD_GEOMOPT)
# option(BUILD_INTS)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER "clang++")

# set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "-g")

add_library(lible SHARED)
add_subdirectory(extern)
add_subdirectory(src)

# miscellaneous
add_definitions(-DBASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/ints/basis_sets")

message("WTH IS GOING ON")
