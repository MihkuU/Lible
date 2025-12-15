#pragma once

#include <array>
#include <vector>

namespace lible::ints
{
    /// Structure containing angular momentum, Gaussian primitive exponents and contraction
    /// coefficients in one shell.
    struct BasisShell
    {
        int l_;
        std::vector<double> exps_;
        std::vector<double> coeffs_;
    };

    /// Type alias for representing a basis set on an atom.
    using basis_shells_t = std::vector<BasisShell>;

    /// Structure containing an atomic number and the basis sets of atomic orbital shells.
    struct BasisAtom
    {
        int atomic_nr_;
        basis_shells_t basis_shells_;
    };

    /// Type alias for representing the basis set on all atoms.
    using basis_atoms_t = std::vector<BasisAtom>;

    /// Structure for representing a shell of atomic orbitals.
    struct Shell
    {
        Shell() = default;

        Shell(const int l, const int z, const size_t idx_atom, const size_t dim_cart,
              const size_t dim_sph, const size_t ofs_cart, const size_t ofs_sph, const size_t idx,
              const std::array<double, 3> &xyz_coords, const std::vector<double> &exps,
              const std::vector<double> &coeffs, const std::vector<double> &norms,
              const std::vector<double> &norms_prim)
            : l_(l), z_(z), dim_cart_(dim_cart), dim_sph_(dim_sph), ofs_cart_(ofs_cart),
              ofs_sph_(ofs_sph), idx_(idx), idx_atom_(idx_atom), xyz_coords_(xyz_coords),
              exps_(exps), coeffs_(coeffs), norms_(norms), norms_prim_(norms_prim)
        {
        }

        /// Angular momentum.
        int l_{};
        /// Atomic number.
        int z_{};

        /// Number of atomic orbitals in a Cartesian basis.
        size_t dim_cart_{};
        /// Number of atomic orbitals in a spherical basis.
        size_t dim_sph_{};
        /// Starting position in the list of atomic orbitals in a Cartesian basis.
        size_t ofs_cart_{};
        /// Starting position in the list of atomic orbitals.
        size_t ofs_sph_{};

        /// Index of the shell in the list of all shells.
        size_t idx_{};
        /// Index of the shell atom in the list of all atoms.
        size_t idx_atom_{};

        /// Coordinates of the atom corresponding to the shell.
        std::array<double, 3> xyz_coords_{};
        /// Exponents of the Gaussian primitives.
        std::vector<double> exps_;
        /// Contraction coefficients of the Gaussian primitives.
        std::vector<double> coeffs_;
        /// Normalization constants of the atomic orbitals in a spherical basis.
        std::vector<double> norms_;
        /// Normalization constants of the Gaussian primitives.
        std::vector<double> norms_prim_;
    };
}
