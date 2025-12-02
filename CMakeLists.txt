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
option(LIBLE_BUILD_DOCS     "Enables building sphinx-documentation" OFF)
option(LIBLE_BUILD_INTS     "Build and install the molecular integrals module" ON)
option(LIBLE_BUILD_STATIC "Enables building the static Lible library" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

### the main target
if(LIBLE_BUILD_STATIC)
	add_library(lible STATIC)
else()
	add_library(lible SHARED)
endif()

add_library(lible::lible ALIAS lible)
target_compile_options(lible PRIVATE -Wall -Wno-unused-local-typedefs)

set_target_properties(lible PROPERTIES PREFIX "")

### Dependencies/Sources 
set(src_dir ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory("${src_dir}/extern")
add_subdirectory("${src_dir}/src")

# OMP
find_package(OpenMP REQUIRED)
target_link_libraries(lible PRIVATE OpenMP::OpenMP_CXX)

# BLAS/LAPACK
find_package(MKL QUIET)
if (MKL_FOUND)
	target_compile_definitions(lible PRIVATE _LIBLE_USE_MKL_)
    target_include_directories(lible PRIVATE ${MKL_ROOT}/include)
	if (LIBLE_BUILD_STATIC)
		target_link_libraries(lible PRIVATE "-m64  -Wl,--start-group ${MKL_ROOT}/lib/libmkl_intel_ilp64.a \
											${MKL_ROOT}/lib/libmkl_sequential.a ${MKL_ROOT}/lib/libmkl_core.a \
											-Wl,--end-group -lpthread -lm -ldl")
	else()
	    target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_gf_lp64.so)
    	target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_core.so)
    	target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_sequential.so)
	endif()
else()
	set(BLA_VENDOR OpenBLAS)
	find_package(BLAS REQUIRED)
	target_link_libraries(lible PRIVATE BLAS::BLAS)
endif()

# Basis sets
## TODO: figure out how to make this stuff work with installation?

if(LIBLE_BUILD_INTS)
	set(lible_main_basis_srcdir "${src_dir}/src/lible/ints/lible_basis")
	target_compile_definitions(lible PUBLIC LIBLE_MAIN_BASIS_DIR="${CMAKE_BINARY_DIR}/lible_basis")
	file(COPY ${lible_main_basis_srcdir} DESTINATION ${CMAKE_BINARY_DIR}/)

	set(lible_aux_basis_srcdir "${src_dir}/src/lible/ints/lible_aux_basis")
	target_compile_definitions(lible PUBLIC LIBLE_AUX_BASIS_DIR="${CMAKE_BINARY_DIR}/lible_aux_basis")
	file(COPY ${lible_aux_basis_srcdir} DESTINATION ${CMAKE_BINARY_DIR}/)
endif()

### Documentation
if(lible_is_top_level AND LIBLE_BUILD_DOCS)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/docs")
endif()

### Installation
if(lible_is_top_level)

	set(lible_basic_headers "${src_dir}/src/lible/utils.hpp")
	list(APPEND lible_basic_headers "${src_dir}/src/lible/types.hpp")
	list(APPEND lible_basic_headers "${src_dir}/src/lible/vectormd.hpp") 

	if(LIBLE_BUILD_INTS)
		file(GLOB lible_ints_headers1 "${src_dir}/src/lible/ints/*.hpp")
		file(GLOB lible_ints_headers2 "${src_dir}/src/lible/ints/oneel/*.hpp")
		file(GLOB lible_ints_headers3 "${src_dir}/src/lible/ints/twoel/*.hpp")
	endif()

	include(GNUInstallDirs)
	install(TARGETS lible
		EXPORT libleTargets
		LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

	install(FILES ${lible_basic_headers}
			DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible)	

	if(LIBLE_BUILD_INTS)	
		install(FILES ${lible_ints_headers1} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible/ints)
		install(FILES ${lible_ints_headers2} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible/ints/oneel)
		install(FILES ${lible_ints_headers3} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lible/ints/twoel)
	endif()

	install(EXPORT libleTargets
		FILE libleTargets.cmake
		NAMESPACE lible::
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible)

	install(FILES "src/libleConfig.cmake"
		DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lible)
endif()

### Testing
if(lible_is_top_level)
	include(CTest)
	enable_testing()
	add_executable(testlible)

	add_subdirectory("${src_dir}/test")

	target_link_libraries(testlible PRIVATE OpenMP::OpenMP_CXX)
endif()
