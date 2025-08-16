cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

set(src_dir "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB basic_sources CONFIGURE_DEPENDS
	"${src_dir}/lible/*.cpp"
	"${src_dir}/lible/*.hpp")

target_sources(lible PRIVATE ${basic_sources})

target_include_directories(lible PUBLIC
	$<BUILD_INTERFACE:${src_dir}/>
	$<INSTALL_INTERFACE:include>)

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
		"${src_dir}/lible/ints/twoel/*.hpp"
		"${src_dir}/lible/ints/twoel/eri_kernels/*.cpp"
		"${src_dir}/lible/ints/twoel/rints_kernels/*.cpp"
		"${src_dir}/lible/ints/twoel/shark_mm_kernels/*.cpp")

	target_sources(lible PRIVATE ${ints_sources})

	target_include_directories(lible PUBLIC 
		$<BUILD_INTERFACE:${src_dir}>
		$<INSTALL_INTERFACE:include>)
endif()
