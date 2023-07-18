cmake_minimum_required(VERSION 3.2)

project(Lible VERSION 0.0.0 LANGUAGES CXX)

if (PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
   message(FATAL_ERROR "In-source builds are not allowed, specify another build directory")
endif()

# TODO: compile which submodule of the library.
# option(BUILD_GEOMOPT)
# option(BUILD_INTS)

# Specify C++ version and compiler.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER "clang++")

# Specify some C++ flags.
# set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "-g")

### The main target of Lible project.
add_library(lible SHARED)

# Add sources from the subdirectories based on which Lible-options are enabled.
add_subdirectory(src/geomopt) # TODO: reduce it to only src/
add_subdirectory(src/ints) # TODO: reduce it to only src/

# Incorporating the dependencies.
add_subdirectory(dep)

# miscellaneous
add_definitions(-DBASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/ints/basis_sets")