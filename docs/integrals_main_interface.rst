
Main Interface
==============

Assume that Lible is properly built/installed and linked against your code. To calculate 
integrals all you have to then do is include the main header ``#include <lible/ints/ints.hpp>``
in your source code. The contents of this header file are documented below.

.. doxygengroup:: IntsMainInterface

.. cpp:function:: bool myMethod(int arg1, std::string arg2)

   A typedef-like declaration of a type. Also,
   eat shit.
.. math::

   \frac{ \sum_{t=0}^{N}f(t,k) }{N} 

.. cpp:type:: overlap_kernel_t = Kernel1<Kernel1Type::overlap, vec2d>

.. cpp:function:: vec2d overlap_kernel_t::operator()(int ipair, const ShellPairData &sp_data)
