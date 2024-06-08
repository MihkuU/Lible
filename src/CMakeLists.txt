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
	file(GLOB ints_sources CONFIGURE_DEPENDS
		"${src_dir}/lible/ints/*.cpp"
		"${src_dir}/lible/ints/*.hpp"
		"${src_dir}/lible/ints/oneel/*.cpp"
		"${src_dir}/lible/ints/oneel/*.hpp"
		"${src_dir}/lible/ints/twoel/*.cpp"
		"${src_dir}/lible/ints/twoel/*.hpp")
	
	if(LIBLE_BUILD_GPUINTS)
		file(GLOB ints_sources_gpu CONFIGURE_DEPENDS
			"${src_dir}/lible/ints/experimental/*.cpp"
			"${src_dir}/lible/ints/experimental/*.hpp")
		list(APPEND ints_sources ${ints_sources_gpu})
	endif()

	target_sources(lible PRIVATE ${ints_sources})

	target_include_directories(lible PUBLIC 
		$<BUILD_INTERFACE:${src_dir}>
		$<INSTALL_INTERFACE:include>)
endif()
