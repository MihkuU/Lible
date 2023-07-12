cmake_minimum_required(VERSION 3.2)

project(LibleTest)

# set(CMAKE_CXX_FLAGS "-g")
set(CMAKE_CXX_FLAGS "-O3")

add_executable(testlible testlible.cpp)
target_include_directories(testlible PUBLIC "/home/alastor/stuff/Lible/include")
target_include_directories(testlible PUBLIC ../build/dep/armadillo-12.4.1/build/include)
target_link_libraries(testlible PUBLIC /home/alastor/stuff/Lible/build/liblible.so)