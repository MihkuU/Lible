#pragma once

#include <map>
#include <string>

#include <lible/ints/shell.hpp>

namespace lible::ints
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
        Structure() = default;

        /** Constructor for initializing a Structure object. */
        Structure(const std::string &basis_set, const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        /** Constructor for initializing a Structure object with enabling the RI approximation. */
        Structure(const std::string &basis_set, const std::string &basis_set_aux,
                  const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        /** Constructor for initializing a structure object with a given custom basis set. */
        Structure(const basis_atoms_t &basis_set, const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        /**
         * Constructor for initializing a structure object with custom main and auxiliary
         * basis sets.
         */
        Structure(const basis_atoms_t &basis_set, const basis_atoms_t &basis_set_aux,
                  const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        // Some getters

        bool getUseRI() const;

        int getMaxL() const;

        int getMaxLAux() const;

        int getZ(int iatom) const;

        int getDimAO() const;

        int getDimAOAux() const;

        int getDimAOCart() const;

        int getDimAOCartAux() const;

        int getNAtoms() const;

        std::array<double, 3> getCoordsAtom(int iatom) const;

        /** Returns for every atom {x, y, z, charge} */
        std::vector<std::array<double, 4>> getZs() const;

        std::vector<Shell> getShells() const;

        std::vector<Shell> getShellsAux() const;

        std::vector<Shell> getShellsL(int l) const;

        std::vector<Shell> getShellsLAux(int l) const;

        std::map<int, std::vector<Shell>> getShellsMap() const;

        std::map<int, std::vector<Shell>> getShellsMapAux() const;

    private:
        bool use_ri{false}; /** Flag for using the RI-approximation. */

        int max_l{}; /** Highest angular momentum among the shells of atomic orbitals. */
        int max_l_aux{}; /** Highest angular momentum among the shells of auxiliary atomic orbitals. */

        int dim_ao{}; /** Number of atomic orbitals. */
        int dim_ao_aux{}; /** Number of auxiliary atomic orbitals. */
        int dim_ao_cart{}; /** Number of atomic orbitals in Cartesian Gaussian basis. */
        int dim_ao_cart_aux{}; /** Number of auxiliary atomic orbitals in Cartesian Gaussian basis. */
        int n_atoms{}; /** Number of atoms in the structure. */

        std::string basis_set; /** Name of the main basis set. */
        std::string basis_set_aux; /** Name of the auxiliary basis set. */

        std::vector<std::array<double, 3>> coords; /** Coordinates of the atoms in Bohr (a.u.). */
        std::vector<int> atomic_nrs; /** Atomic numbers. */

        std::vector<Shell> shells; /** All of the shells corresponding to the main basis set. */
        std::vector<Shell> shells_aux; /** All of the shells corresponding to the auxiliary basis set. */

        std::map<int, std::vector<Shell>> shells_map;
        /** Shells corresponding to the main basis set for each angular momentum. */
        std::map<int, std::vector<Shell>> shells_map_aux;
        /** Shells corresponding to the auxiliary basis set for each angular momentum. */
    };
}
