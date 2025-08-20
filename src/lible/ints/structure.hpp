#pragma once

#include <map>
#include <string>

#include <lible/ints/shell.hpp>

namespace lible
{
    namespace ints
    {
        /** \defgroup structure */

        typedef std::pair<std::vector<double>, std::vector<double>> shell_exps_coeffs_t;

        /**
         * \ingroup structure
         * Typedef for bundling shell Gaussian exponents and contraction coefficients for each
         * angular momentum per atom.
         */
        typedef std::map<int, std::map<int, std::vector<shell_exps_coeffs_t>>> basis_atoms_t;

        /**
         * \struct Structure structure.hpp <lible/ints/structure.hpp>
         * \ingroup structure
         *
         * Structure for representing the atomic structure. The main data carrier for calculating
         * the integrals. Note that to use the RI approximation, the Structure object has to be
         * initialized by providing an auxiliary basis set name. By default, RI is disabled.
         */
        // TODO: add a flattened list of shells
        // TODO: replace the 'coords_angstrom' with 'coords_bohr'        
        struct Structure
        {
            /** Constructur for initializing a structure object. */
            Structure();

            /** Constructor for initializing a Structure object. */
            Structure(const std::string &basis_set, const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            /** Constructor for initializing a Structure object. RI approximation is disabled. */
            Structure(const std::string &basis_set, const std::vector<std::string> &elements,
                      const std::vector<double> &coords_angstrom);

            /** Constructor for initializing a Structure object with enabling the RI approximation. */
            Structure(const std::string &basis_set, const std::string &basis_set_aux,
                      const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            /** Constructor for initializing a Structure object with enabling the RI approximation. */
            Structure(const std::string &basis_set, const std::string &basis_set_aux,
                      const std::vector<std::string> &elements,
                      const std::vector<double> &coords_angstrom);

            /** Constructor for initializing a structure object with a given custom basis set. */
            Structure(const basis_atoms_t &basis_set_custom, const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            /**
             * Constructor for initializing a structure object with a given main basis set
             * and custom auxiliary basis set.
             */
            Structure(const std::string &basis_set, const basis_atoms_t &basis_set_custom_aux,
                      const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            /**
             * Constructor for initializing a structure object with a custom main basis set
             * and a given auxiliary basis set.
             */
            Structure(const basis_atoms_t &basis_set_custom, const std::string &basis_set_aux,
                      const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            /**
             * Constructor for initializing a structure object with custom main and auxiliary
             * basis sets.
             */
            Structure(const basis_atoms_t &basis_set_custom, const basis_atoms_t &basis_set_custom_aux,
                      const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            // Some getters

            bool getUseRI() const;

            int getMaxL() const;

            int getMaxLAux() const;

            int getZ(const int iatom) const;

            int getDimAO() const;

            int getDimAOAux() const;

            int getDimAOCart() const;

            int getDimAOCartAux() const;

            int getNAtoms() const;

            std::array<double, 3> getCoordsAtom(const int iatom) const;

            /** Returns for every atom {x, y, z, charge} */
            std::vector<std::array<double, 4>> getZs() const;

            std::vector<Shell> getShellsL(const int l) const;

            std::vector<Shell> getShellsLAux(const int l) const;

            std::map<int, std::vector<Shell>> getShells() const;

            std::map<int, std::vector<Shell>> getShellsAux() const;

        private:
            bool use_ri{false}; /** Flag for using the RI-approximation. */

            int max_l{};     /** Highest angular momentum among the shells of atomic orbitals. */
            int max_l_aux{}; /** Highest angular momentum among the shells of auxiliary atomic orbitals. */

            int dim_ao{};          /** Number of atomic orbitals. */
            int dim_ao_aux{};      /** Number of auxiliary atomic orbitals. */
            int dim_ao_cart{};     /** Number of atomic orbitals in Cartesian Gaussian basis. */
            int dim_ao_cart_aux{}; /** Number of auxiliary atomic orbitals in Cartesian Gaussian basis. */
            int n_atoms{};         /** Number of atoms in the structure. */

            std::string basis_set;     /** Name of the main basis set. */
            std::string basis_set_aux; /** Name of the auxiliary basis set. */

            std::vector<double> coords;                    /** Rolled out coordinates of atoms.*/
            std::vector<int> atomic_nrs;                   /** Atomic numbers. */
            std::vector<std::array<double, 3>> coords_xyz; /** Coordinates of the atoms. */
            std::vector<std::string> elements;             /** Symbols of the atoms. */

            std::map<int, std::vector<Shell>> shells;     /** Shells corresponding to the main basis set for each angular momentum. */
            std::map<int, std::vector<Shell>> shells_aux; /** Shells corresponding to the auxiliary basis set for each angular momentum. */

            /**
             * Constructs the shells for the given basis. Used for both the main and auxiliary
             * basis.
             */
            void constructShells(const basis_atoms_t &basis_atoms, int &max_l, int &dim_ao,
                                 int &dim_ao_cart, std::map<int, std::vector<Shell>> &shells);
        };
    }
}