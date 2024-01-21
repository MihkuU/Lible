cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

file(GLOB util_sources "util/*.cpp")
target_sources(lible PRIVATE ${util_sources})
target_include_directories(lible PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/util>
	$<INSTALL_INTERFACE:include>
)

file(GLOB davidson_sources "davidson/*.cpp")
target_sources(lible PRIVATE ${davidson_sources})
target_include_directories(lible PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/davidson>
	$<INSTALL_INTERFACE:include>
)

file(GLOB geomopt_sources "geomopt/*.cpp")
target_sources(lible PRIVATE ${geomopt_sources})
target_include_directories(lible PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/geomopt>
	$<INSTALL_INTERFACE:include>
)

file(GLOB gugaci_sources "gugaci/*.cpp")
target_sources(lible PRIVATE ${gugaci_sources})
target_include_directories(lible PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gugaci>
	$<INSTALL_INTERFACE:include>
)

file(GLOB ints_sources "ints/*.cpp")
file(GLOB ints_1el_sources "ints/one_el_ints/*.cpp")
target_sources(lible PRIVATE ${ints_sources})
target_sources(lible PRIVATE ${ints_1el_sources})
target_include_directories(lible PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ints>
	$<INSTALL_INTERFACE:include>
)

#target_sources(le PRIVATE experimental/overlap_ints_gpu.cpp)
#target_sources(le PRIVATE experimental/structure_gpu.cpp)
