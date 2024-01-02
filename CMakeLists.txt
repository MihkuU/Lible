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
### options
option(BUILD_DAVIDSON "" ON)
option(BUILD_GEOMOPT  "" ON)
option(BUILD_GUGASCI  "" ON)
option(BUILD_INTS     "" ON)
option(BUILD_GPUINTS  "" OFF)

option(USE_MPI "" ON)

if(USE_MPI)
	find_package(MPI REQUIRED)
	set(CMAKE_CXX_COMPILER mpic++)
	add_definitions(-D_USE_MPI_)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_COMPILER "hipcc") # TODO: use hipcc only if "ints"-module is built
#set(CMAKE_CXX_COMPILER mpic++)

### the main target
add_library(lible SHARED)
set_target_properties(lible PROPERTIES PREFIX "")

add_subdirectory(extern)
add_subdirectory(src)

find_package(OpenMP REQUIRED)
target_link_libraries(lible PUBLIC ${OpenMP_CXX_LIBRARIES})

# miscellaneous
add_definitions(-DBASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/ints/basis_sets")

### installation
set(lible_public_headers "src/util/lible/util.h")
list(APPEND lible_public_headers "src/util/lible/types.h")
list(APPEND lible_public_headers "src/gugaci/lible/gci.h")

set_target_properties(lible PROPERTIES
	PUBLIC_HEADER "${lible_public_headers}"
)

include(GNUInstallDirs)
install(TARGETS lible
	EXPORT LibleTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible
)

install(EXPORT LibleTargets
	FILE LibleTargets.cmake
	NAMESPACE Lible::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible
)

install(FILES "src/LibleConfig.cmake"
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible
)
