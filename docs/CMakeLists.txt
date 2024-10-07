#[=[ 
	Based on:
	1) https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/
	2) https://coderefinery.github.io/documentation/

    Some other examples:
	1) http://ceres-solver.org/index.html
    2) https://rabauke.github.io/mpl/html/
    3) https://xtb-docs.readthedocs.io/en/latest/
	4) https://github.com/roman-ellerbrock/QuTree
#]=]

set(lible_header_dir "${PROJECT_SOURCE_DIR}/src/lible/ints")
set(lible_ints_headers "${lible_header_dir}/ints.hpp" 
					   "${lible_header_dir}/shell.hpp"
					   "${lible_header_dir}/shell_pair_data.hpp"
					   "${lible_header_dir}/structure.hpp")

### The Doxygen part
find_package(Doxygen REQUIRED)

set(doxygen_input_dir ${PROJECT_SOURCE_DIR}/src/lible/ints/)
set(doxygen_output_dir ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(doxygen_index_file ${doxygen_output_dir}/xml/index.xml)
set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(doxyfile_out ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

configure_file(${doxyfile_in} ${doxyfile_out} @ONLY)

file(MAKE_DIRECTORY ${doxygen_output_dir})
add_custom_command(OUTPUT ${doxygen_index_file}
	DEPENDS               ${lible_ints_headers}
	COMMAND               ${DOXYGEN_EXECUTABLE} ${doxyfile_out}
	MAIN_DEPENDENCY       ${doxyfile_out} ${doxyfile_in}
	WORKING_DIRECTORY     ${doxygen_output_dir}
	COMMENT               "Generating docs"
	VERBATIM)

add_custom_target(Doxygen ALL DEPENDS ${doxygen_index_file})

### The Sphinx part
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
		${CMAKE_CURRENT_SOURCE_DIR}/integrals.rst
		${CMAKE_CURRENT_SOURCE_DIR}/ints/definitions.rst
		${CMAKE_CURRENT_SOURCE_DIR}/ints/ints.rst
		${CMAKE_CURRENT_SOURCE_DIR}/ints/shell.rst
		${CMAKE_CURRENT_SOURCE_DIR}/ints/spdata.rst
		${CMAKE_CURRENT_SOURCE_DIR}/ints/structure.rst
		${CMAKE_CURRENT_SOURCE_DIR}/ints/tutorials.rst
		${doxygen_index_file}
	MAIN_DEPENDENCY   ${sphinx_source}/conf.py
	COMMENT           "Generating documentation with Sphinx")

add_custom_target(Sphinx ALL DEPENDS ${sphinx_index_file})

#[=[
# Add an install target to install the docs
include(GNUInstallDirs)
install(DIRECTORY ${SPHINX_BUILD}
DESTINATION ${CMAKE_INSTALL_DOCDIR})
#]=]
