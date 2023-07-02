cmake_minimum_required(VERSION 3.2)

project(Lible
		VERSION 0.0.0
		LANGUAGES CXX)

# TODO: compile which submodule of the library.

# Specify C++ version and compiler.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER "clang++")

# Specify some C++ flags.
set(CMAKE_CXX_FLAGS "-O3")

# The main target of Lible project.
add_library(lible SHARED)

# Incorporating some dependencies.
target_include_directories(lible PUBLIC "${PROJECT_SOURCE_DIR}/dep/json-nlohman-develop/include")
target_include_directories(lible PUBLIC "${PROJECT_SOURCE_DIR}/dep/fmt-master/include")

# Add sources from the subdirectories based on which Lible-options are enabled.
add_subdirectory(src/geomopt)
add_subdirectory(src/ints)

# miscellaneous
add_definitions(-DBASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/ints/basis_sets")