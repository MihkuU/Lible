# Lible

## Set up 

### Preqrequisites
- Compiler that enables C++20 features such as `std::format` or `consteval`. For example, GCC-13 or newer.
A comprehensive overview of appropriate compilers can be seen at [C++ compiler support](https://en.cppreference.com/w/cpp/compiler_support.html).
- Some form of BLAS/LAPACK. Recommended either [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) or
[IntelMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html).

### Installation  
Building and installing Lible is necessary primarily if you intend to work with the source code, for example, to 
add new features. If you want to just integrate the library in your code and use its features, then it is advisable 
to use the CMake `FetchContent` feature that is explained below. For the separate build, follow the steps:
  - Clone the repo as usual.  
  - Then run:
  ```bash
    1. cmake -S . -B build -DCMAKE_BUILD_TYPE=<Specify Debug or Release>
    2. cmake --build build/ -j <nprocs>
    3. cmake --install build/ --prefix "<full path to the chosen installdir, can be build/>"
   ```

### CMake integration
  1. The 3. step in Installation ensures that Lible can be incorporated in your CMake project using the 'find_package()' command. 
  To do that, write in your 'CMakeLists.txt' file:
  ```
    find_package(lible REQUIRED)
    target_link_libraries(<your target> PRIVATE lible::lible)
  ```
  For the `find_package()` command to work, your CMake configuration has to find Lible. This can be facilitated by pointing the
  CMake search path to the Lible installation location:
  ```
    -DCMAKE_PREFIX_PATH=<full path to the Lible installdir>
  ```
  
  2. Lible can be incorporated in a CMake project directly using `FetchContent`. Write in your `CMakeLists.txt` file:
  ```
    include(FetchContent)

    FetchContent_Declare(lible
      GIT_REPOSITORY https://github.com/MihkuU/Lible
      GIT_TAG <XXX>)

    FetchContent_MakeAvailable(lible)
    target_link_libraries(<your target> PRIVATE lible::lible)
  ```
  This approach downloads the library and integrates it directly in your build. For more details see 
  [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html).

### Tests
Lible builds a test program whenever it is the [top level project](https://cmake.org/cmake/help/latest/variable/PROJECT_IS_TOP_LEVEL.html), 
i.e., it is built separately. The test cases are implemented via the [CTest](https://cmake.org/cmake/help/latest/module/CTest.html) module. 
To run the tests, simply run
```
ctest --test-dir build/
```
