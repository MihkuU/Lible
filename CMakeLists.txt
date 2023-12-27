cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(Lible VERSION 0.0.0 LANGUAGES CXX)

if (PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
   message(FATAL_ERROR "In-source builds are not allowed, specify another build directory")
endif()

if(NOT((CMAKE_BUILD_TYPE STREQUAL "Debug") OR (CMAKE_BUILD_TYPE STREQUAL "Release")))
    message(FATAL_ERROR "Please specify either 'Debug' or 'Release' build type: \
                         -DCMAKE_BUILD_TYPE=<build type>")
endif()

# TODO: compile which submodule of the library.
option(BUILD_GEOMOPT "" ON)
option(BUILD_INTS    "" ON)
option(BUILD_GPUINTS "" OFF)

option(USE_MPI "" ON)

if(USE_MPI)
	find_package(MPI REQUIRED)
	#set(CMAKE_CXX_COMPILER mpic++)
	#execute_process(COMMAND which mpic++ OUTPUT_VARIABLE FOO)
	#set(CMAKE_CXX_COMPILER "${FOO}")
	set(CMAKE_CXX_COMPILER /home/alastor/libs/openmpi-4.1.1/build/bin/mpic++)
	add_definitions(-D_USE_MPI_)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_COMPILER "hipcc") # TODO: use hipcc only if "ints"-module is built
#set(CMAKE_CXX_COMPILER mpic++)

add_library(le SHARED) # lib-prefix will be prepended
add_subdirectory(extern)
add_subdirectory(src)

# miscellaneous
add_definitions(-DBASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/ints/basis_sets")


# installation
include(GNUInstallDirs)
message("CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
message("CMAKE_INSTALL_INCLUDEDIR = ${CMAKE_INSTALL_INCLUDEDIR}/lible")
install(TARGETS le DESTINATION "${CMAKE_INSTALL_LIBDIR}/lible")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/include/lible.h"
		DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/lible")

get_target_property(LE_PUBLIC_HEADERS le INTERFACE_SOURCES)
message("LE_PUBLIC_HEADERS: ${LE_PUBLIC_HEADERS}")
