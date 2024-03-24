.. _Installation:

Installation
============

Prerequisites
-------------

At least C++17-capable compiler. To check whether your compiler is sufficiently new, see the
`C++ compiler support <https://en.cppreference.com/w/cpp/compiler_support>`_. A version of BLAS
and LAPACK must be available to compile the library. Some possible BLAS/LAPACK implementations are:

* `OpenBLAS <https://www.openblas.net/>`_
* `Intel oneAPI MKL <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html>`_

Optionally, if distributed memory parallelization (e.g., across multiple nodes) is desired, a version
of MPI must be installed and accessible. An `OpenMPI <https://www.open-mpi.org/>`_ version since 3.0 
is recommended. To toggle the use of MPI, please refer to the **Compile options** section below.

CMake integration
-----------------

If you are using CMake, which is recommended, Lible can be incorporated into your C++ project with 
ease using the `FetchContent <https://cmake.org/cmake/help/latest/module/FetchContent.html>`_ module
in modern CMake. All you have to do is write the following lines in the CMakeLists.txt file of your 
project::

   # Enable FetchContent
   include(FetchContent)

   # Declare the Lible dependencies to be fetched
   FetchContent_Declare(lible
       GIT_REPOSITORY https://github.com/MihkuU/Lible
       GIT_TAG        xxxxxx)

   # Make the dependencies available
   FetchContent_MakeAvailable(lible)

   # Link against your target. Necessary to actually build Lible
   target_link_libraries(your_target PRIVATE lible::lible)

This automatically downloads and compiles Lible, making its sources available for use in your 
project. No additional effort is required, and you may proceed with using what the library offers.
Note that a specific git tag for ``xxxxxx`` has to be provided to choose a particular version of a
commit from the repo, e.g.

.. note::
   This approach supports setting or disabiling options that enable certain features from the library,
   such as using MPI. An example will be provided below.

Manual installation
-------------------

This option is suitable if you want to work with Lible separately, e.g., as a developer

.. _My target:

Compile options
---------------
