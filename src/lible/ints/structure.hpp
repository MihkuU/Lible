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

            /** Constructor for initializing a structure object with a given custom basis sets. */
            Structure(const basis_atoms_t& basis_set_custom, const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            /**
             * Constructor for initializing a structure object with a given custom main basis sets
             * and auxiliary.
             */
            Structure(const basis_atoms_t& basis_set_custom, const basis_atoms_t& basis_set_custom_aux,
                      const std::vector<int> &atomic_nrs,
                      const std::vector<double> &coords_angstrom);

            bool getUseRI() const;

            int getMaxL() const;

            int getMaxLAux() const;

            int getZ(const size_t iatom) const;

            size_t getDimAO() const;

            size_t getDimAOAux() const;

            size_t getDimAOCart() const;

            size_t getDimAOCartAux() const;

            size_t getNAtoms() const;

            std::array<double, 3> getCoordsAtom(const size_t iatom) const;

            std::vector<Shell> getShellsL(const int l) const;

            std::vector<Shell> getShellsLAux(const int l) const;

            std::map<int, std::vector<Shell>> getShells() const;

            std::map<int, std::vector<Shell>> getShellsAux() const;

        private:
            /** Flag for using the RI-approximation. */
            bool use_ri{false};

            /** */
            bool prescreen_built{false};

            /** */
            bool prescreen_built_ri{false};

            /** Highest angular momentum among the shells of atomic orbitals. */
            int max_l{};

            /** Highest angular momentum among the shells of auxiliary atomic orbitals. */
            int max_l_aux{};

            /** Number of atomic orbitals. */
            size_t dim_ao{};

            /** Number of auxiliary atomic orbitals. */
            size_t dim_ao_aux{};

            /** Number of atomic orbitals in Cartesian Gaussian basis. */
            size_t dim_ao_cart{};

            /** Number of auxiliary atomic orbitals in Cartesian Gaussian basis. */
            size_t dim_ao_cart_aux{};

            /** Number of atoms in the structure. */
            size_t n_atoms{};

            /** Name of the main basis set. */
            std::string basis_set;

            /** Name of the auxiliary basis set. */
            std::string basis_set_aux;            

            /** Rolled out coordinates of atoms.*/
            std::vector<double> coords;

            /** Atomic numbers. */
            std::vector<int> atomic_nrs;

            /** Coordinates of the atoms. */
            std::vector<std::array<double, 3>> coords_xyz;

            /** Symbols of the atoms. */
            std::vector<std::string> elements;

            /** Shells corresponding to the main basis set for each atomic number. */
            std::map<int, std::vector<Shell>> shells;

            /** Shells corresponding to the auxiliary basis set for each atomic number. */
            std::map<int, std::vector<Shell>> shells_aux;

            /**
             * Constructs the shells for the given basis. Used for both the main and auxiliary
             * basis.
             */
            void constructShells(const basis_atoms_t &basis_atoms, int &max_l, size_t &dim_ao,
                                 size_t &dim_ao_cart, std::map<int, std::vector<Shell>> &shells);
        };
    }
}