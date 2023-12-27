file(GLOB ints_sources "*.cpp")
target_sources(le PRIVATE ${ints_sources})
target_sources(le PRIVATE experimental/overlap_ints_gpu.cpp)
target_sources(le PRIVATE experimental/structure_gpu.cpp)
