#pragma once

#include <map>
#include <string>

#include <lible/ints/shell.hpp>

namespace lible::ints
{
    /// Structure for representing the atomic structure and various data for calculating intgrals.
    /// Contains the xyz-coordinates, shells, etc.
    struct Structure
    {
        /// Default constructor.
        Structure() = default;

        /// Constructor for initializing the object without the RI approximiation.
        Structure(const std::string &basis_set, const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        /// Constructor for initializing the object with the RI approximation.
        Structure(const std::string &basis_set, const std::string &basis_set_aux,
                  const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        /// Constructor for initializing the object without the RI approximiation.
        Structure(const basis_atoms_t &basis_set, const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        /// Constructor for initializing the object with the RI approximation.
        Structure(const basis_atoms_t &basis_set, const basis_atoms_t &basis_set_aux,
                  const std::vector<int> &atomic_nrs,
                  const std::vector<std::array<double, 3>> &coords_angstrom);

        // Some getters

        /// Returns the flag for using RI.
        bool getUseRI() const;

        /// Returns the maximum angular momentum in the main basis set shells.
        int getMaxL() const;

        /// Returns the maximum angular momentum in the auxiliary basis set shells.
        int getMaxLAux() const;

        /// Returns the atomic number (nuclear charge) at atom `iatom`.
        int getZ(size_t iatom) const;

        /// Returns the number of (spherical) atomic orbitals from the main basis set.
        size_t getDimAO() const;

        /// Returns the number of (spherical) atomic orbitals from the auxiliary basis set.
        size_t getDimAOAux() const;

        /// Returns the number of Cartesian atomic orbitals from the main basis set.
        size_t getDimAOCart() const;

        /// Returns the number of Cartesian atomic orbitals from the auxiliary basis set.
        size_t getDimAOCartAux() const;

        /// Returns the number of atoms.
        size_t getNAtoms() const;

        /// Returns the xyz-coordinates (in a.u.) at atom `iatom`.
        std::array<double, 3> getCoordsAtom(size_t iatom) const;

        /// Returns for every atom {x, y, z, charge}.
        std::vector<std::array<double, 4>> getZs() const;

        /// Returns all the main basis set shells.
        std::vector<Shell> getShells() const;

        /// Returns all the auxiliary basis set shells.
        std::vector<Shell> getShellsAux() const;

        /// Returns all the main basis shells with angular momentum `l`.
        std::vector<Shell> getShellsL(int l) const;

        /// Returns all the auxiliary basis shells with angular momentum `l`.
        std::vector<Shell> getShellsLAux(int l) const;

        /// Returns the main basis shells map, {l, {shells}}.
        std::map<int, std::vector<Shell>> getShellsMap() const;

        /// Returns the auxiliary basis shells map, {l, {shells}}.
        std::map<int, std::vector<Shell>> getShellsMapAux() const;

    private:
        /// Flag for using the RI approximation.
        bool use_ri_{false};

        /// Highest angular momentum among the shells of atomic orbitals.
        int max_l_{};
        /// Highest angular momentum among the shells of auxiliary atomic orbitals.
        int max_l_aux_{};

        /// Number of atomic orbitals.
        size_t dim_ao_{};
        /// Number of auxiliary atomic orbitals.
        size_t dim_ao_aux_{};
        /// Number of atomic orbitals in the Cartesian Gaussian basis.
        size_t dim_ao_cart_{};
        /// Number of auxiliary atomic orbitals in the Cartesian Gaussian basis.
        size_t dim_ao_cart_aux_{};
        /// Number of atoms in the structure.
        size_t n_atoms_{};

        /// Coordinates of the atoms in Bohr (a.u.).
        std::vector<std::array<double, 3>> coords_;
        /// Atomic numbers.
        std::vector<int> atomic_nrs_;

        /// Name of the main basis set.
        std::string basis_set_;
        /// Name of the auxiliary basis set.
        std::string basis_set_aux_;
        /// All the shells corresponding to the main basis set.
        std::vector<Shell> shells_;
        /// All the shells corresponding to the auxiliary basis set.
        std::vector<Shell> shells_aux_;
        /// Shells corresponding to the main basis set for each angular momentum.
        std::map<int, std::vector<Shell>> shells_map_;
        /// Shells corresponding to the auxiliary basis set for each angular momentum.
        std::map<int, std::vector<Shell>> shells_map_aux_;
    };
}
