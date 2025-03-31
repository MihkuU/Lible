#pragma once

#include <array>
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
            Shell(const int l, const int z, const size_t dim_cart, const size_t dim_sph,
                  const size_t pos, const size_t pos_cart, const std::array<double, 3> &xyz_coords,
                  const std::vector<double> &coeffs, const std::vector<double> &coeffs_raw,
                  const std::vector<double> &exps, const std::vector<double> &norms)
                : l(l), z(z), dim_cart(dim_cart), dim_sph(dim_sph), pos(pos), pos_cart(pos_cart),
                  xyz_coords(xyz_coords), coeffs(coeffs), coeffs_raw(coeffs_raw), exps(exps),
                  norms(norms)
            {
            }

            int l; /** Angular momentum. */
            int z; /** Atomic number. */

            size_t dim_cart; /** Number of atomic orbitals in Cartesian basis. */
            size_t dim_sph;  /** Number of atomic orbitals in spherical basis. */
            size_t pos;      /** Starting position in the list of atomic orbitals. */
            size_t pos_cart; /** Starting position in the list of atomic orbitals in Cartesian basis. */

            std::array<double, 3> xyz_coords; /** Coordinates of the atom corresponding to the shell. */

            std::vector<double> coeffs;     /** Contraction coefficients with primitive norms multiplied into. */
            std::vector<double> coeffs_raw; /** Contraction coefficients without primitive norms. */
            std::vector<double> exps;       /** Exponents of the Gaussian primitives. */
            std::vector<double> norms;      /** Normalization constants of the atomic orbitals in spherical basis. */
        };

        /**
         * \ingroup shell
         *
         * Calculates the norms of the shell atomic orbitals in spherical basis.
         */
        std::vector<double> calcShellNorms(const int l, const std::vector<double> &coeffs,
                                           const std::vector<double> &exps);
    }
}