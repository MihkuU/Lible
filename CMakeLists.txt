cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(lible
	VERSION      0.0.0 
	LANGUAGES    CXX
	DESCRIPTION  "A C++ quantum chemistry library with various usable gadgets."
	HOMEPAGE_URL "https://github.com/MihkuU/Lible")

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
   message(FATAL_ERROR "In-source builds are not allowed, specify another build directory!")
endif()

### Check whether it's the main project for later usage
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
	set(lible_is_top_level TRUE)
else()
	set(lible_is_top_level FALSE)
endif()

### Options

option(LIBLE_BUILD_DOCS "Enables building sphinx-documentation" OFF)
option(LIBLE_USE_MPI    "Enables compiling with the MPI wrapper" ON)
option(LIBLE_USE_LIBXSMM "Enables using libxsmm for fast small matrix multiplication" ON)

### Modules
option(LIBLE_BUILD_DAVIDSON "Build and install the Davidson module" ON)
option(LIBLE_BUILD_GEOMOPT  "Build and install the Geometry optimization module" ON)
option(LIBLE_BUILD_GCI      "Build and install the GUGA-SCI module" ON)
option(LIBLE_BUILD_INTS     "Build and install the molecular integrals module" ON)
option(LIBLE_BUILD_GPUINTS  "Build and install the GPU molecular integrals module" OFF)

if(LIBLE_USE_MPI)
	find_package(MPI REQUIRED)
	set(CMAKE_CXX_COMPILER mpic++)
	add_definitions(-D_LIBLE_USE_MPI_)
else()
	set(CMAKE_CXX_COMPILER g++)
endif()

if(LIBLE_USE_LIBXSMM)
	add_definitions(-D_LIBLE_USE_LIBXSMM)
endif()

if(LIBLE_BUILD_GPUINTS)
	add_definitions(-D_LIBLE_USE_HIP_)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

### the main target
add_library(lible SHARED)
add_library(lible::lible ALIAS lible)

set_target_properties(lible PROPERTIES PREFIX "")

set(src_dir ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory("${src_dir}/extern")
add_subdirectory("${src_dir}/src")

# OMP
find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
target_link_libraries(lible PRIVATE ${OpenMP_CXX_LIBRARIES})

# BLAS/LAPACK
find_package(MKL)
if (MKL_FOUND)
	add_definitions(-D_LIBLE_USE_MKL_)
    target_include_directories(lible PRIVATE ${MKL_ROOT}/include)
    target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_gf_lp64.so)
    target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_core.so)
    target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_sequential.so)
else()
	set(BLA_VENDOR OpenBLAS)
	find_package(BLAS REQUIRED)
	target_link_libraries(lible PRIVATE ${BLAS_LIBRARIES})
endif()

# miscellaneous # TODO: manage this
add_definitions(-DBASIS_DIR="${src_dir}/src/lible/ints/basis")

### Documentation
if(lible_is_top_level AND LIBLE_BUILD_DOCS)
	add_subdirectory("${src_dir}/docs")
endif()

### installation
if(lible_is_top_level)

	# TODO: put in a file and include for better organization?
	set(lible_basic_headers "${src_dir}/src/lible/util.hpp")
	list(APPEND lible_basic_headers "${src_dir}/src/lible/types.hpp")
	list(APPEND lible_basic_headers "${src_dir}/src/lible/vec2d.hpp") 
	list(APPEND lible_basic_headers "${src_dir}/src/lible/vec3d.hpp") 
	list(APPEND lible_basic_headers "${src_dir}/src/lible/vec4d.hpp") 

	if(LIBLE_BUILD_DAVIDSON)
	endif()

	if(LIBLE_BUILD_GCI)
		list(APPEND lible_gci_headers "${src_dir}/src/lible/gugaci/gci.hpp")
		list(APPEND lible_gci_headers "${src_dir}/src/lible/gugaci/gci_settings.hpp")
		list(APPEND lible_gci_headers "${src_dir}/src/lible/davidson/davidson_settings.hpp")
	endif()

	if(LIBLE_BUILD_INTS)
		list(APPEND lible_ints_headers "${src_dir}/src/lible/ints/ints.hpp")
		list(APPEND lible_ints_headers "${src_dir}/src/lible/ints/shell.hpp")
		list(APPEND lible_ints_headers "${src_dir}/src/lible/ints/structure.hpp")
	endif()

	include(GNUInstallDirs)
	install(TARGETS lible
		EXPORT libleTargets
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

	install(FILES ${lible_basic_headers}
			DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible)	

	if(LIBLE_BUILD_INTS)
		install(FILES ${lible_ints_headers}
				DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible/ints)
	endif()

	if(LIBLE_BUILD_GCI)
		install(FILES ${lible_gci_headers}
				DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible/gci)
	endif()

	install(EXPORT libleTargets
		FILE libleTargets.cmake
		NAMESPACE lible::
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible)

	install(FILES "src/libleConfig.cmake"
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible)
endif()

