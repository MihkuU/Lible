# Lible

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

The 3. step in Installation ensures that the Lible can be conveniently incorporated in your CMake project using the 'find_package()' function call. In your 'CMakeLists.txt' file you can write:
  ```
  find_package(Lible REQUIRED)
  target_link_libraries(YourProject PRIVATE Lible::lible)
  ```
For the find_package() to work, your CMake project configuration has to find Lible. That means, you need to provide the path to the directory, where it was installed:
```
  -DCMAKE_PREFIX_PATH=<Path to where you installed Lible in 'Installation'>
```
<details>
<summary>
Short summary
</summary>

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod
tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam,
quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse
cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non
proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
</details>
