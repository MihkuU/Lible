cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(Lible VERSION 0.0.0 LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
   message(FATAL_ERROR "In-source builds are not allowed, specify another build directory")
endif()

if(NOT((CMAKE_BUILD_TYPE STREQUAL "Debug") OR (CMAKE_BUILD_TYPE STREQUAL "Release")))
    message(FATAL_ERROR "Please specify either 'Debug' or 'Release' build type: \
                         -DCMAKE_BUILD_TYPE=<build type>")
endif()

### TMP
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/extern/arma/build/lib/") 
# ...SOB that was hard. Tried by following: https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
# TODO: figure out a better organization for it in the future.
### TMP

### options
option(BUILD_DAVIDSON "" ON)
option(BUILD_GEOMOPT  "" ON)
option(BUILD_GCI      "" ON)
option(BUILD_INTS     "" ON)
option(BUILD_GPUINTS  "" OFF)

option(USE_MPI "" OFF)

if(USE_MPI)
	find_package(MPI REQUIRED)
	set(CMAKE_CXX_COMPILER mpic++)
	add_definitions(-D_USE_MPI_)
else()
	set(CMAKE_CXX_COMPILER g++)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_COMPILER "hipcc") # TODO: use hipcc only if "ints"-module is built
#set(CMAKE_CXX_COMPILER mpic++)

### the main target
add_library(lible SHARED)
set_target_properties(lible PROPERTIES PREFIX "")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/extern")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")

find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
target_link_libraries(lible PUBLIC ${OpenMP_CXX_LIBRARIES})

# find_package(unordered_dense CONFIG REQUIRED)
# target_link_libraries(lible PRIVATE unordered_dense::unordered_dense)

# miscellaneous
add_definitions(-DBASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/ints/basis_sets")

### installation
##### TODO: put the installation details into an appropriate file, e.g. lible_installation.cmake
set(lible_public_headers "src/util/lible/util.h")
list(APPEND lible_public_headers "src/util/lible/types.h")
list(APPEND lible_public_headers "src/util/lible/vec2d.h")
list(APPEND lible_public_headers "src/util/lible/vec3d.h")
list(APPEND lible_public_headers "src/util/lible/vec4d.h")

if(BUILD_GCI)
	list(APPEND lible_public_headers "src/gugaci/lible/gci.hpp")
	list(APPEND lible_public_headers "src/gugaci/lible/gci_settings.hpp")
	list(APPEND lible_public_headers "src/davidson/lible/davidson_settings.h")
endif()

if(BUILD_INTS)
	list(APPEND lible_public_headers "src/ints/lible/ints.hpp")
	list(APPEND lible_public_headers "src/ints/lible/shell.hpp")
	list(APPEND lible_public_headers "src/ints/lible/structure.hpp")
endif()

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
