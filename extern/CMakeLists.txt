cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Setting up some external libraries
# List of libraries:
#   - armadillo, https://arma.sourceforge.net/
#   - ankerl, https://github.com/martinus/unordered_dense
#   - fmt, https://github.com/fmtlib/fmt
#   - libxsmml, https://github.com/libxsmm/libxsmm
#   - mpl, https://github.com/rabauke/mpl

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)
execute_process(COMMAND "${CMAKE_COMMAND}" . WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/)
execute_process(COMMAND "${CMAKE_COMMAND}" --build . WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/)

# arma
target_compile_definitions(lible PRIVATE ARMA_DONT_USE_WRAPPER)
target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/include)
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib")
    target_link_libraries(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib/libarmadillo.a)
elseif(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib64")
    target_link_libraries(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/arma/build/lib64/libarmadillo.a)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_definitions(lible PRIVATE ARMA_NO_DEBUG)
endif()

# fmt
target_compile_definitions(lible PRIVATE FMT_HEADER_ONLY)
target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/fmt/include)

# ankerl-unordered dense
find_package(unordered_dense CONFIG REQUIRED PATHS ${CMAKE_CURRENT_BINARY_DIR}/unordered_dense/build/lib/cmake/unordered_dense)
target_link_libraries(lible PRIVATE unordered_dense::unordered_dense)

# mpl
if(LIBLE_USE_MPI)
	target_include_directories(lible PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mpl/build/include)
endif()

# libxsmml
target_include_directories(lible PRIVATE /home/alastor/libs/libxsmm-main/build/include)
target_include_directories(lible PRIVATE /home/alastor/libs/libxsmm-main/include)
#[=[
include(FetchContent)
FetchContent_Declare(
  xsmm
  GIT_REPOSITORY https://github.com/libxsmm/libxsmm
  GIT_TAG b8b085e
)
#FetchContent_GetProperties(xsmm)
#if(NOT xsmm_POPULATED)
  #FetchContent_Populate(xsmm)
#endif()
FetchContent_MakeAvailable(xsmm)

set(LIBXSMMROOT ${xsmm_SOURCE_DIR})
file(GLOB _GLOB_XSMM_SRCS LIST_DIRECTORIES false CONFIGURE_DEPENDS ${LIBXSMMROOT}/src/*.c)
list(REMOVE_ITEM _GLOB_XSMM_SRCS ${LIBXSMMROOT}/src/libxsmm_generator_gemm_driver.c)
set(XSMM_INCLUDE_DIRS ${LIBXSMMROOT}/include)

add_library(xsmm STATIC ${_GLOB_XSMM_SRCS})
set_target_properties(xsmm PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(xsmm PUBLIC ${XSMM_INCLUDE_DIRS})
target_compile_definitions(xsmm PUBLIC
  LIBXSMM_DEFAULT_CONFIG
)

target_link_libraries(lible PRIVATE xsmm)
#]=]
