.. Lible documentation master file, created by
   sphinx-quickstart on Fri Mar 22 14:36:43 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====
Lible
=====

`Lible <https://github.com/MihkuU/Lible/>`_ is a modern C++ library that aims to enhance the 
development of quantum chemical methods. Its primary goal is to provide thin-grained utilities
in the form of functions and types that cover various concepts in quantum chemistry:

* calculation of molecular integrals
* geometry optimization
* configuration interaction
* etc.

Therefore, Lible is designed in a modular manner, where each of these ideas is represented by
a separate programming interface, or *namespace*, to be specific. 

Another goal is to make the esoteric knowledge of quantum chemists more accessible to the people,
and especially newcomers, who are working in this field. To achieve that, the documentation of the
project aims to be `tutorial-driven <https://coderefinery.github.io/documentation/summary/>`_ 
as much as possible. A bare "Wham bam, thank you, Mam!" automatically generated doxygen C++ API is 
not enough here. Thus, the overall documentation is split into three parts:

1. Introduction. Gives a basic overview of setting up the library, its usage and the available features.
2. User guide. Explains the theoretical background behind the modules and how to use them.
3. Developer API. Provides the full documentation of the library utilities.

.. toctree::
   :hidden:
   :maxdepth: 1

   integrals
