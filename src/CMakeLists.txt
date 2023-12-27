cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

add_subdirectory(davidson)
add_subdirectory(geomopt)
add_subdirectory(gugasci)
add_subdirectory(ints)
add_subdirectory(util)

target_include_directories(le PUBLIC util)
#target_include_directories(le PUBLIC davidson)
target_include_directories(le PRIVATE davidson)
