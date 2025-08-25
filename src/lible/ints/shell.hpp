#pragma once

#include <array>
#include <map>
#include <utility>
#include <vector>

namespace lible
{
    //  TODO: clean it up!
    //  Here we deal with atomic orbital shells. We follow conventions outlined in https://iodata.readthedocs.io/en/latest/basis.html.
    //  Definitions:
    //    Shell - set of basis functions with the same angular momentum, contraction coefficients and
    //    exponents of gaussian primitives.
    namespace ints
    {
        /** \defgroup shell */

        /**
         * \struct Shell shell.hpp <lible/ints/shell.hpp>
         * \ingroup shell
         *
         * Structure for representing a shell of atomic orbitals. Groups together various
         * data that defines a shell: angular momentum, contraction coefficients, contraction
         * exponents etc.
         */
        struct Shell
        {
            /** The constructor. */
            Shell(const int l, const int z, const int atom_idx, const int dim_cart,
                  const int dim_sph, const int ofs_cart, const int ofs_sph, const int idx,
                  const std::array<double, 3> &xyz_coords, const std::vector<double> &exps,
                  const std::vector<double> &coeffs, const std::vector<double> &norms,
                  const std::vector<double> &norms_prim)
                : l(l), z(z), atom_idx(atom_idx), dim_cart(dim_cart), dim_sph(dim_sph),
                  ofs_cart(ofs_cart), ofs_sph(ofs_sph), idx(idx), xyz_coords(xyz_coords),
                  exps(exps), coeffs(coeffs), norms(norms), norms_prim(norms_prim)
            {
            }

            int l; /** Angular momentum. */
            int z; /** Atomic number. */

            int atom_idx; /** Index of the shell atom in the list of all atoms. */
            int dim_cart; /** Number of atomic orbitals in Cartesian basis. */
            int dim_sph;  /** Number of atomic orbitals in spherical basis. */
            int ofs_cart; /** Starting position in the list of atomic orbitals in Cartesian basis. */
            int ofs_sph;  /** Starting position in the list of atomic orbitals. */            

            int idx; /** Index of the shell in the list of all shells. */

            std::array<double, 3> xyz_coords; /** Coordinates of the atom corresponding to the shell. */

            std::vector<double> exps;       /** Exponents of the Gaussian primitives. */
            std::vector<double> coeffs;     /** Contraction coefficients of the Gaussian primitives. */            
            std::vector<double> norms;      /** Normalization constants of the atomic orbitals in spherical basis. */
            std::vector<double> norms_prim; /** Normalization constants of the Gaussian primitives. */
        };

        using shell_exps_coeffs_t = std::pair<std::vector<double>, std::vector<double>>;
        using basis_atom_t = std::map<int, std::vector<shell_exps_coeffs_t>>;
        using basis_atoms_t = std::map<int, basis_atom_t>;

        /**
         * \ingroup shell
         *
         * Calculates the norms of the shell atomic orbitals in spherical basis.
         */
        std::vector<double> calcShellNorms(const int l, const std::vector<double> &coeffs,
                                           const std::vector<double> &exps,
                                           const std::vector<double> &primitive_norms);

        /** */
        std::vector<Shell> constructShells(const basis_atoms_t &basis_atoms,
                                           const std::vector<int> &atomic_nrs,
                                           const std::vector<std::array<double, 3>> &coords_atoms);
    }
}
