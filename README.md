# Lible
Lible is a toolbox containing multiple libraries for implementing quantum chemical methods. The focus is on providing a clear and modern C++ programming interface while being as efficient as possible in the backend. The name LIBLE corresponds to Library for Linear Efficiency. It doesn't mean everything here is asymptotically linear scaling, the name just sounded good.

** Prerequisites **
  - A C++17-capable compiler
  - A form of BLAS/LAPACK is probably required. IntelMKL is preferable, otherwise get OpenBLAS.

** Installation **
  - Pull the repo as usual
  - Then run:
  ```
    1. cmake -S . -B build -DCMAKE_BUILD_TYPE=<Specify Debug or Release>
    2. cmake --build build/ -j <nprocs>
    3. cmake --install build/ --prefix "<full path to the build/-directory>"
  ```
** Using Lible **

The 3. step in the Installation ensures that the Lible can be conveniently incorporated in your CMake project using the 'find_package()' function call. In your 'CMakeLists.txt' file you could write:
  ```
  find_package(Lible REQUIRED)
  target_link_libraries(YourProject PRIVATE Lible::lible)
  ```
For the find_package to work, your CMake project configuration has to find Lible. That means, you need to provide the path to the build directory of Lible, where it was installed:
```
  -DCMAKE_PREFIX_PATH=<Path to where you built Lible during 'Installation'>
```
