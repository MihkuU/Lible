cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(lible
	VERSION      0.0.0 
	LANGUAGES    CXX
	DESCRIPTION  "A C++ quantum chemistry library with various usable gadgets."
	HOMEPAGE_URL "https://github.com/MihkuU/Lible"
)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
   message(FATAL_ERROR "In-source builds are not allowed, specify another build directory!")
endif()

### options
option(BUILD_DAVIDSON "" ON)
option(BUILD_GEOMOPT  "" OFF)
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

### the main target
add_library(lible SHARED)
add_library(lible::lible ALIAS lible)

set_target_properties(lible PROPERTIES PREFIX "")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/extern")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")

# OMP
find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
target_link_libraries(lible PUBLIC ${OpenMP_CXX_LIBRARIES})

# BLAS/LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

target_link_libraries(lible PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# miscellaneous # TODO: manage this
add_definitions(-DBASIS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/src/ints/basis_sets")

### installation
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	set(lible_is_top_level TRUE)
else()
	set(lible_is_top_level FALSE)
endif()

if(lible_is_top_level)
	set(lible_public_headers "src/util/lible/util.hpp")
	list(APPEND lible_public_headers "src/util/lible/types.hpp")
	list(APPEND lible_public_headers "src/util/lible/vec2d.hpp") # TODO: put in one header file
	list(APPEND lible_public_headers "src/util/lible/vec3d.hpp") # TODO: put in one header file 
	list(APPEND lible_public_headers "src/util/lible/vec4d.hpp") # TODO: put in one header file

	if(BUILD_DAVIDSON)
	endif()

	if(BUILD_GCI)
		list(APPEND lible_public_headers "src/gugaci/lible/gci.hpp")
		list(APPEND lible_public_headers "src/gugaci/lible/gci_settings.hpp")
		list(APPEND lible_public_headers "src/davidson/lible/davidson_settings.hpp")
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
		EXPORT libleTargets
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
		PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible
	)

	install(EXPORT libleTargets
		FILE libleTargets.cmake
		NAMESPACE lible::
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible
	)

	install(FILES "src/libleConfig.cmake"
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible
	)
endif()
