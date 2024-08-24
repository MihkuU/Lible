cmake_minimum_required(VERSION 3.2)

project(LibleTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER "clang++")

# set(CMAKE_CXX_FLAGS "-g")
set(CMAKE_CXX_FLAGS "-O3")

add_executable(testlible testlible.cpp)

target_include_directories(testlible PUBLIC "/home/alastor/stuff/Lible/include")
target_include_directories(testlible PUBLIC "/home/alastor/stuff/Lible/build/dep/armadillo-12.4.1/build/include")
target_link_libraries(testlible PUBLIC "/home/alastor/stuff/Lible/build/dep/armadillo-12.4.1/build/lib/libarmadillo.so")
target_link_libraries(testlible PUBLIC "/home/alastor/stuff/Lible/build/liblible.so")
