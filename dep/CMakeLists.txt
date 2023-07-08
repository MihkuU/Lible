# Setting up some external libraries
# List of libraries:
#   - armadillo, https://arma.sourceforge.net/
#   - fmt, https://github.com/fmtlib/fmt
#   - json, https://github.com/nlohmann/json

set(src_dep_dir "${CMAKE_SOURCE_DIR}/dep")
set(build_dep_dir "${CMAKE_BINARY_DIR}/dep")
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${build_dep_dir})

# armadillo
message("1.  Building armadillo library")
set(arma_archive "armadillo-12.4.1.tar.xz")
set(arma_extract "armadillo-12.4.1")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${src_dep_dir}/${arma_archive} ${build_dep_dir}/${arma_archive})
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${arma_archive} WORKING_DIRECTORY ${build_dep_dir})
execute_process(COMMAND ${CMAKE_COMMAND} -E rm ${arma_archive} WORKING_DIRECTORY ${build_dep_dir})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${build_dep_dir}/${arma_extract}/build)
message("   ... Configuring")
execute_process(COMMAND ${CMAKE_COMMAND} . -DCMAKE_INSTALL_PREFIX:PATH=${build_dep_dir}/${arma_extract}/build 
                WORKING_DIRECTORY ${build_dep_dir}/${arma_extract} OUTPUT_FILE arma_conf.log)
message("   ... Building")
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${build_dep_dir}/${arma_extract} OUTPUT_FILE arma_build.log)
message("   ... Installing")
execute_process(COMMAND ${CMAKE_COMMAND} --install . WORKING_DIRECTORY ${build_dep_dir}/${arma_extract} OUTPUT_FILE arma_install.log)
target_include_directories(lible PUBLIC ${build_dep_dir}/${arma_extract}/build/include)

# fmt
message("2.  Setting up fmt-library (header-only)")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${src_dep_dir}/fmt-master.zip ${build_dep_dir}/fmt-master.zip)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf fmt-master.zip WORKING_DIRECTORY ${build_dep_dir})
execute_process(COMMAND ${CMAKE_COMMAND} -E rm fmt-master.zip WORKING_DIRECTORY ${build_dep_dir})
target_include_directories(lible PUBLIC ${build_dep_dir}/fmt-master/include)

# json
message("3.  Setting for json-library (header-only)")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${src_dep_dir}/json-develop.zip ${build_dep_dir}/json-develop.zip)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf json-develop.zip WORKING_DIRECTORY ${build_dep_dir})
execute_process(COMMAND ${CMAKE_COMMAND} -E rm json-develop.zip WORKING_DIRECTORY ${build_dep_dir})
target_include_directories(lible PUBLIC ${build_dep_dir}/json-develop/include)