set(tests_dir ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(testlible PRIVATE 
	${tests_dir}/main.cpp 
	${tests_dir}/tests.cpp)

target_include_directories(testlible PRIVATE ${tests_dir})
target_include_directories(testlible PRIVATE ${src_dir}/src/)

add_test(NAME "lible::ints::numCartesians()" COMMAND ${CMAKE_BINARY_DIR}/testlible "numCartesians")
add_test(NAME "lible::ints::numSphericals()" COMMAND ${CMAKE_BINARY_DIR}/testlible "numSphericals")
