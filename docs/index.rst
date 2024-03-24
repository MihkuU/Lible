.. Lible documentation master file, created by
   sphinx-quickstart on Fri Mar 22 14:36:43 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Lible overview
==============

Lible is a modern C++ library that aims to enhance the development of quantum chemical methods. 
Its primary goal is to provide thin-grained utilities in the form of functions and types that cover
various concepts in quantum chemistry:

* calculation of molecular integrals
* geometry optimization
* configuration interaction
* etc...

Therefore, Lible is designed in a modular manner, where each of these ideas is represented by
a separate programming interface, or *namespace*, to be specific. 

Another goal is to make the esoteric knowledge of quantum chemists more accessible to the people,
and especially especially newcomers, who are working in the this field. To achieve that, the 
documentation of the project aims to be `tutorial-driven <https://coderefinery.github.io/documentation/summary/>`_ 
as much as possible. A mere "Wham bam, thank you, Mam!" automatically generated doxygen C++ API is not 
enough here. Thus, the overall documentation is split into two parts:

1. User guide
2. Developer guide

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Introduction

   installation
   basic_usage
   modules

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: User guide

   dox_test

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Developer guide

..
   Docs
   ====
   .. doxygenstruct:: lible::ints::test
   :members:
