#[=[ 
	Based on the blogpost:
	1) https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/
	2) https://coderefinery.github.io/documentation/
#]=]

# The main header corresponding to the lible-ints module
set(lible_ints_main_header "${PROJECT_SOURCE_DIR}/src/ints/lible/ints.hpp")

### Handle the Doxygen part
find_package(Doxygen REQUIRED)

set(doxygen_input_dir ${PROJECT_SOURCE_DIR}/src/ints/lible/ints.hpp)
set(doxygen_output_dir ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(doxygen_index_file ${CMAKE_CURRENT_SOURCE_DIR}/html/index.html)
set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(doxyfile_out ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

configure_file(${doxyfile_in} ${doxyfile_out} @ONLY)

file(MAKE_DIRECTORY ${doxygen_output_dir})
add_custom_command(OUTPUT ${doxygen_index_file}
	DEPENDS               ${lible_ints_main_header}
	COMMAND               ${DOXYGEN_EXECUTABLE} ${doxyfile_out}
	WORKING_DIRECTORY     ${CMAKE_CURRENT_SOURCE_DIR}
	MAIN_DEPENDENCY       ${doxyfile_out} ${doxyfile_in}
	COMMENT               "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${doxygen_index_file})

### Handle the Sphinx part
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(Sphinx REQUIRED)

set(sphinx_source ${CMAKE_CURRENT_SOURCE_DIR})
set(sphinx_build ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(sphinx_index_file ${sphinx_build}/index.html)

# Only regenerate Sphinx when:
# - Doxygen has rerun
# - Our doc files have been updated
# - The Sphinx config has been updated
add_custom_command(OUTPUT ${sphinx_index_file}
	COMMAND
		${SPHINX_EXECUTABLE} -b html
		-Dbreathe_projects.Lible=${doxygen_output_dir}/xml
		${sphinx_source} ${sphinx_build}
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	DEPENDS
		${CMAKE_CURRENT_SOURCE_DIR}/index.rst
		${doxygen_index_file}
	MAIN_DEPENDENCY   ${sphinx_source}/conf.py
	COMMENT           "Generating documentation with Sphinx")

add_custom_target(Sphinx ALL DEPENDS ${sphinx_index_file})

# Add an install target to install the docs
include(GNUInstallDirs)
install(DIRECTORY ${SPHINX_BUILD}
DESTINATION ${CMAKE_INSTALL_DOCDIR})

#[=[
add_custom_target(Sphinx ALL
	COMMAND           
		${SPHINX_EXECUTABLE} -b html
		-Dbreathe_projects.Lible=${doxygen_output_dir}/xml # Tell Breathe where to find the Doxygen output
		${sphinx_source} ${sphinx_build}
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT           "Generating docs with Sphinx")
#]=]
