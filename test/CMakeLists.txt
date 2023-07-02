cmake_minimum_required(VERSION 3.2)

project(LibleTest)

add_executable(testlible testlible.cpp)
target_include_directories(testlible PUBLIC "/home/alastor/stuff/Lible/include")
target_link_libraries(testlible PUBLIC /home/alastor/stuff/Lible/build/liblible.so)