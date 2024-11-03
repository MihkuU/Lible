cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

set(src_dir "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB basic_sources CONFIGURE_DEPENDS
	"${src_dir}/lible/*.cpp"
	"${src_dir}/lible/*.hpp")

target_sources(lible PRIVATE ${basic_sources})

target_include_directories(lible PUBLIC
	$<BUILD_INTERFACE:${src_dir}/>
	$<INSTALL_INTERFACE:include>)

#[=[
if(LIBLE_BUILD_DAVIDSON)
	file(GLOB davidson_sources "davidson/*.cpp")
	target_sources(lible PRIVATE ${davidson_sources})
	target_include_directories(lible PUBLIC 
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/davidson>
		$<INSTALL_INTERFACE:include>)
endif()
#]=]

if(LIBLE_BUILD_GEOMOPT)
	file(GLOB geomopt_sources CONFIGURE_DEPENDS
		"geomopt/*.cpp"
		"geomopt/*.hpp")

	target_sources(lible PRIVATE ${geomopt_sources})

	target_include_directories(lible PUBLIC 
		$<BUILD_INTERFACE:${src_dir}>
		$<INSTALL_INTERFACE:include>)
endif()

if(LIBLE_BUILD_GCI)
	file(GLOB gugaci_sources CONFIGURE_DEPENDS 
		"${src_dir}/lible/gugaci/*.cpp"
		"${src_dir}/lible/gugaci/*.hpp"
		"${src_dir}/lible/davidson/*.cpp"
		"${src_dir}/lible/davidson/*.hpp")

	target_sources(lible PRIVATE ${gugaci_sources})

	target_include_directories(lible PUBLIC 
		$<BUILD_INTERFACE:${src_dir}>
		$<INSTALL_INTERFACE:include>)
endif()

if(LIBLE_BUILD_INTS)
	target_compile_options(lible PRIVATE -fconstexpr-ops-limit=2147483647)
	file(GLOB ints_sources CONFIGURE_DEPENDS
		"${src_dir}/lible/ints/*.cpp"
		"${src_dir}/lible/ints/*.hpp"
		"${src_dir}/lible/ints/oneel/*.cpp"
		"${src_dir}/lible/ints/oneel/*.hpp"
		"${src_dir}/lible/ints/twoel/*.cpp"
		"${src_dir}/lible/ints/twoel/eri4_kernels/*.cpp"
		"${src_dir}/lible/ints/twoel/*.cpp")

	#set_source_files_properties(${ints_sources} PROPERTIES LANGUAGE ${GPU_RUNTIME})
	target_sources(lible PRIVATE ${ints_sources})
	
	if(LIBLE_BUILD_GPUINTS)
   		#target_compile_definitions(lible PRIVATE _LIBLE_USE_HIP_) # Does it work here???

		set(CMAKE_CXX_STANDARD 17)
		set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

		#file(GLOB ints_sources_gpu CONFIGURE_DEPENDS
		#	"${src_dir}/lible/ints/gpu/*.cpp"
		#	"${src_dir}/lible/ints/gpu/*.hpp")
		file(GLOB ints_sources_gpu CONFIGURE_DEPENDS
			"${src_dir}/lible/ints/gpu/*.cpp")

		add_library(lible_gpuints OBJECT)
		target_compile_options(lible_gpuints PRIVATE -fPIC) # Customize for AMD?
		target_include_directories(lible_gpuints PRIVATE ${src_dir})

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

    	set_source_files_properties(${ints_sources_gpu} PROPERTIES LANGUAGE ${GPU_RUNTIME})
	    set_source_files_properties(${ints_sources_gpu} PROPERTIES COMPILE_FLAGS ${HIPCC_FLAGS})

		target_sources(lible_gpuints PRIVATE ${ints_sources_gpu})

		# fmt
		target_compile_definitions(lible_gpuints PRIVATE FMT_HEADER_ONLY)
		target_include_directories(lible_gpuints PRIVATE ${CMAKE_BINARY_DIR}/extern/fmt/include)

	endif()

	target_include_directories(lible PUBLIC 
		$<BUILD_INTERFACE:${src_dir}>
		$<INSTALL_INTERFACE:include>)
endif()
