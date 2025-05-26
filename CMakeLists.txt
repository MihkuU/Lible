cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(lible
	VERSION      0.0.0 
	#LANGUAGES    CXX
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

# Modules
option(LIBLE_BUILD_DAVIDSON "Build and install the Davidson module" OFF)
option(LIBLE_BUILD_DOCS     "Enables building sphinx-documentation" OFF)
option(LIBLE_BUILD_GEOMOPT  "Build and install the Geometry optimization module" OFF)
option(LIBLE_BUILD_GCI      "Build and install the GUGA-SCI module" OFF)
option(LIBLE_BUILD_INTS     "Build and install the molecular integrals module" ON)
option(LIBLE_BUILD_GPUINTS  "Build and install the GPU molecular integrals module" OFF)

### Options
option(LIBLE_USE_MPI    "Enables compiling with the MPI wrapper" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(LIBLE_BUILD_GPUINTS)
    # Based on the boilerplate that I found from AMDs HIP course
    if(NOT CMAKE_GPU_RUNTIME)
        set(GPU_RUNTIME "ROCM" CACHE STRING "Switches between ROCM and CUDA")
    else()
        set(GPU_RUNTIME "${CMAKE_GPU_RUNTIME}" CACHE STRING "Switches between ROCM and CUDA")
    endif()

    # Should only be ROCM or CUDA, but allowing HIP because it is the currently built-in option
    # Select with e.g., -DGPU_RUNTIME=ROCM
    set(GPU_RUNTIMES "ROCM" "CUDA" "HIP")
    if(NOT "${GPU_RUNTIME}" IN_LIST GPU_RUNTIMES)
        set(ERROR_MESSAGE
            "GPU_RUNTIME is set to \"${GPU_RUNTIME}\".\nGPU_RUNTIME must be either HIP, ROCM, or CUDA.")
        message(FATAL_ERROR ${ERROR_MESSAGE})
    endif()

    # GPU_RUNTIME should really be ROCM for AMD GPUs, so manually resetting to HIP if ROCM is selected
    if(${GPU_RUNTIME} MATCHES "ROCM")
        set(GPU_RUNTIME "HIP")
    endif()
    set_property(CACHE GPU_RUNTIME PROPERTY STRINGS ${GPU_RUNTIMES})
    enable_language(${GPU_RUNTIME})
    set(CMAKE_${GPU_RUNTIME}_EXTENSIONS OFF)
    set(CMAKE_${GPU_RUNTIME}_STANDARD_REQUIRED ON)

    set(ROCMCC_FLAGS "${ROCMCC_FLAGS} ")
    set(CUDACC_FLAGS "${CUDACC_FLAGS} ")
    if(${GPU_RUNTIME} MATCHES "HIP")
        set(HIPCC_FLAGS "${ROCMCC_FLAGS}")
    else()
        set(HIPCC_FLAGS "${CUDACC_FLAGS}")
    endif()
endif()

### the main target
add_library(lible SHARED)
add_library(lible::lible ALIAS lible)
target_compile_options(lible PRIVATE -Wall -Wno-unused-local-typedefs)

set_target_properties(lible PROPERTIES PREFIX "")

### Dependencies/Sources 
set(src_dir ${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory("${src_dir}/extern")
add_subdirectory("${src_dir}/src")

if(LIBLE_USE_MPI)
	find_package(MPI REQUIRED)
	target_compile_definitions(lible PRIVATE _LIBLE_USE_MPI_)
	target_link_libraries(lible PRIVATE MPI::MPI_CXX)
endif()

if(LIBLE_BUILD_GPUINTS)
	target_compile_definitions(lible PRIVATE _LIBLE_USE_HIP_)
	target_link_libraries(lible PRIVATE lible_gpuints)
endif()

# OMP
find_package(OpenMP REQUIRED)
target_compile_options(lible PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(lible PRIVATE ${OpenMP_CXX_FLAGS})

# BLAS/LAPACK
find_package(MKL QUIET)
if (MKL_FOUND)
	target_compile_definitions(lible PRIVATE _LIBLE_USE_MKL_)
    target_include_directories(lible PRIVATE ${MKL_ROOT}/include)
    target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_gf_lp64.so)
    target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_core.so)
    target_link_libraries(lible PRIVATE ${MKL_ROOT}/lib/intel64/libmkl_sequential.so)
else()
	set(BLA_VENDOR OpenBLAS)
	find_package(BLAS REQUIRED)
	target_link_libraries(lible PRIVATE BLAS::BLAS)
endif()

# miscellaneous # TODO: manage this
target_compile_definitions(lible PRIVATE BASIS_DIR="${src_dir}/src/lible/ints/basis")
target_compile_definitions(lible PRIVATE AUX_BASIS_DIR="${src_dir}/src/lible/ints/aux_basis")

### Documentation
if(lible_is_top_level AND LIBLE_BUILD_DOCS)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/docs")
endif()

### Installation
if(lible_is_top_level)

	set(lible_basic_headers "${src_dir}/src/lible/utils.hpp")
	list(APPEND lible_basic_headers "${src_dir}/src/lible/types.hpp")
	list(APPEND lible_basic_headers "${src_dir}/src/lible/vectormd.hpp") 

	if(LIBLE_BUILD_DAVIDSON)
	endif()

	if(LIBLE_BUILD_GCI)
		list(APPEND lible_gci_headers "${src_dir}/src/lible/gugaci/gci.hpp")
		list(APPEND lible_gci_headers "${src_dir}/src/lible/gugaci/gci_settings.hpp")
		list(APPEND lible_gci_headers "${src_dir}/src/lible/davidson/davidson_settings.hpp")
	endif()

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

### Testing
include(CTest)
enable_testing()
add_executable(testlible)

add_subdirectory("${src_dir}/test")

