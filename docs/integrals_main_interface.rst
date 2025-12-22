
Main Interface
==============

Assume that Lible is properly built/installed and linked against your code. To use the library for 
calculating integrals, the main header ``#include <lible/ints/ints.hpp>`` has to be then included. 
This header file constitutes the main interface. The main interface contains inclusions of other 
header files as well and may be illustrated diagrammatically:

.. graphviz::

   digraph {
      rankdir="LR";
      graph [fontname="Verdana", fontsize="12"];
      node [fontname="Verdana", fontsize="12"];
      edge [fontname="Sans", fontsize="9"];

      vectormd [label="<lible/vectormd.hpp>"]
      types [label="<lible/types.hpp>"]
      ints [label="<lible/ints/ints.hpp>"]
      shell [label="<lible/ints/shell.hpp>"]
      spdata [label="<lible/ints/shell_pair_data.hpp>"]
      boysfun [label="<lible/ints/boys_function.hpp>"]
      utils [label="<lible/ints/utils.hpp>"]
      structure [label="<lible/ints/structure.hpp>"]
      erikernels [label="<lible/ints/eri_kernels.hpp>"]

      vectormd -> types 
      types -> ints 
      shell -> spdata
      spdata -> ints
      boysfun -> ints
      utils -> ints
      structure -> ints
      erikernels -> ints
   }

Below are is described the programming utilities from these files. But before that, another diagram 
is useful for a conceptual understanding of how to calculate integrals with Lible:

.. graphviz::

   digraph {
      rankdir="LR"
      graph [fontname="Verdana", fontsize="12"];
      node [fontname="Verdana", fontsize="12"];
      edge [fontname="Sans", fontsize="9"];

      geom
      basisset
      shells
      spdata
      integrals
      
      geom -> shells
      basisset -> shells
      shells -> spdata
      spdata -> integrals
   }

A simplified path to calculating integrals, that hides dealing with shells and the shell pair data,
is via the structure object:

.. graphviz::

   digraph {
      rankdir="LR"
      graph [fontname="Verdana", fontsize="12"];
      node [fontname="Verdana", fontsize="12"];
      edge [fontname="Sans", fontsize="9"];

      geom
      basisset
      structure
      integrals 
     
      geom -> structure
      basisset -> structure
      structure -> integrals
   }

For prototyping and simple usage this scheme is to be preferrer. For more control and large scale 
application, the previous scheme is the way to go.

Basis and Shells
----------------

The information content of the previous diagrams is a simplicistic representation of what happens 
in quantum chemical calculations and more specifically, in calculating integrals which is a central
task in the former. Having chosen a molecular geometry and a basis set, the basic building blocks, 
i.e., the shells, can be constructed. This section provides an overview of the contents in 
``<lible/ints/shell.hpp>``.

.. cpp:struct:: lible::ints::BasisShell
   
   Structure representing the basis set of a specific atomic orbital shell.

   .. cpp:var:: int l_; 

      Angular momentum of the shell.

   .. cpp:var:: std::vector<double> exps_;

      Gaussian primitive exponents. 

   .. cpp:var:: std::vector<double> coeffs_;

      Contraction coefficients of the Gaussian primitives.      

.. cpp:type:: basis_shells_t = std::vector<lible::ints::BasisShell>

   Type representing a list of basis sets on shells.

.. cpp:struct:: BasisAtom 

   Structure representing the basis set on a specific atom.

   .. cpp:var:: int atomic_nr_; 

      Atomic number of the atom.

   .. cpp:var:: basis_shells_t basis_shells_;

      List of shell basis sets on the atom.

.. cpp:type:: basis_atoms_t = std::vector<BasisAtom>

   Type representing a list of basis sets on atoms. This object can be used to represent the entire 
   molecular basis set (main or auxiliary).


.. cpp:struct:: lible::ints::Shell

   Structure representing an atomic orbital shell. Contains essential data for calculating 
   integrals. rambleramble.

   .. cpp:var:: int l_; 

      Angular momentum of the shell.

   .. cpp:var:: int z_;
      
      Atomic number of the shell's atom.

   .. cpp:var:: size_t dim_cart_;

      Number of Cartesian atomic orbitals.

   .. cpp:var:: size_t dim_sph_;

      Number of spherical atomic orbitals.

   .. cpp:var:: size_t ofs_cart_;

      Starting position of the current shell's Cartesian atomic orbitals in the list of all AOs.

   .. cpp:var:: size_t ofs_sph_;

      Starting position of the current shell's spherical atomic orbitals in the list of all AOs.

   .. cpp:var:: size_t idx_;

      Index of the current shell in the list of all shells.

   .. cpp:var:: size_t idx_atom_;

      Index of the current shell's atom in the list of all atoms.

   .. cpp:var:: std::array<double, 3> xyz_coords_;

      Cartesian coordinates of the shell's atom.

   .. cpp:var:: std::vector<double> exps_;

      Gaussian primitive exponents.

   .. cpp:var:: std::vector<double> coeffs_;

      Contraction coefficients of the Gaussian primitives.      

   .. cpp:var:: std::vector<double> norms_;

      Normalization coefficients of the spherical atomic orbitals.

   .. cpp:var:: std::vector<double> norms_prim_;

      Normalization coefficients of the Gaussian primitives.
  