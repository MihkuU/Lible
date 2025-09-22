# Lible

## Set up 

### Preqrequisites
- A compiler that enables C++20 features such as `std::format` or `consteval`. For example, GNU-13 or newer.
For the appropriate compilers, see [C++ compiler support](https://en.cppreference.com/w/cpp/compiler_support.html).
- A form of BLAS/LAPACK. Recommended either [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) or
[IntelMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html).

### Installation  
Building and installing Lible is necessary primarily if you intend to work with the source code, for example, to 
add new features. If you want to just integrate the library in your code and use its features, then it is advisable 
to use the CMake `FetchContent` feature that is explained below. For the separate build, follow the steps:
  - Clone the repo as usual.  
  - Then run:
  ```
    1. cmake -S . -B build -DCMAKE_BUILD_TYPE=<Specify Debug or Release>
    2. cmake --build build/ -j <nprocs>
    3. cmake --install build/ --prefix "<full path to the chosen installdir, can be build/"
   ```

### CMake integration

** Using Lible **

The 3. step in Installation ensures that the Lible can be conveniently incorporated in your CMake project using the 'find_package()' function call. In your 'CMakeLists.txt' file you can write:
  ```
  find_package(Lible REQUIRED)
  target_link_libraries(YourProject PRIVATE Lible::lible)
  ```
For the find_package() to work, your CMake project configuration has to find Lible. That means, you need to provide the path to the directory, where it was installed:
```
  -DCMAKE_PREFIX_PATH=<Path to where you installed Lible in 'Installation'>
```
