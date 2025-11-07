#pragma once

#include <lible/ints/structure.hpp>

namespace lible
{
    namespace ints
    {
        /** \defgroup spdata */    

        /**
         * \ingroup spdata
         *
         * Structure for representing rolled out data of shells that is utilized in the
         * calculation of integrals. By rolled out it is meant that the data in the shells
         * is placed into memory contiguously. Desired data can be accessed with offset indices
         * and dimensions.
         */
        struct ShellData
        {
            /** */
            ShellData(int l, const std::vector<Shell> &shells);

            int l{};            /** Angular momentum. */
            int n_shells{};     /** Number of shells. */
            int n_primitives{}; /** Total number of Gaussian primitives. */

            std::vector<double> coeffs; /** Contraction coefficients with the primitive norms multiplied into. */
            std::vector<double> coords; /** Coordinates of atoms corresponding to shells. */
            std::vector<double> exps;   /** Exponents of the Gaussian primitives. */
            std::vector<double> norms;  /** Normalization constants of the atomic orbitals in spherical basis. */

            std::vector<int> atomic_idxs;     /** Indices of atoms involved in the shells. */
            std::vector<int> cdepths;         /** Contraction depths corresponding to each shell. */
            std::vector<int> coffsets;        /** Offsets of the contraction data (coeffs and exps) for each shell.*/
            std::vector<int> offsets_ecoeffs; /** Offsets of the spherical Hermite expansion coefficients. */
            std::vector<int> offsets_norms;   /** Offsets of the shell atomic orbital norms. */
            std::vector<int> offsets_sph;     /** Offsets of the atomic orbital positions in the list of all atomic orbitals. */
            std::vector<int> shell_idxs;      /** Indices of the involved shells in the list of all shells. */
        };

        /**
         * \ingroup spdata
         *
         * Structure for representing contiguously rolled out data of shell pairs. The data of
         * each shell in a pair is placed side-by-side to each other in a corresponding vector.
         * If \f$l_a \neq l_b\f$, all of the shells are taken. When \f$l_a = l_b\f$, we enforce
         * \f$i \geq j\f$ where \f$i, j\f$ refer to the shell indices. Contains a few additional
         * gadgets that might be useful.
         */
        struct ShellPairData
        {
            ShellPairData(bool use_symm, int la, int lb, const std::vector<Shell> &shells_a,
                          const std::vector<Shell> &shells_b);

            bool uses_symm{}; /** Flag indicating whether symmetry is used or not. */

            int la{};           /** Angular momentum of the left shell. */
            int lb{};           /** Angular momentum of the right shell. */
            int n_pairs{};      /** Total number of shell pairs. */
            int n_prim_pairs{}; /** Total number of Gaussian primitive pairs. */

            std::vector<double> coeffs; /** Contraction coefficients including the primitive norms. */
            std::vector<double> coords; /** Coordinates of shells in the shell pairs. */
            std::vector<double> exps;   /** Exponents of the Gaussian primitives. */
            std::vector<double> norms;  /** Normalization consts of the atomic orbitals in spherical basis. */

            std::vector<int> atomic_idxs;            /** Indices of atoms involved in the shell pairs. */
            std::vector<int> cdepths;                /** Contraction depths corresponding to each shell in a shell pair. */
            std::vector<int> coffsets;               /** Offsets of the contraction data (coeffs and exps) for each shell in a shell pair.*/
            std::vector<int> offsets_cart;           /** Offsets of the atomic orbital (cartesian) positions in the list of all atomic orbitals. */
            std::vector<int> offsets_ecoeffs;        /** Offsets of the spherical Hermite expansion coefficients. */
            std::vector<int> offsets_ecoeffs_deriv1; /** Offsets of the 1st derivative spherical Hermite expansion coefficients. */
            std::vector<int> offsets_ecoeffs_deriv2; /** Offsets of the 2nd derivative spherical Hermite expansion coefficients. */
            std::vector<int> offsets_norms;          /** Offsets of the shell atomic orbital norms. */
            std::vector<int> offsets_sph;            /** Offsets of the atomic orbital (spherical) positions in the list of all atomic orbitals. */
            std::vector<int> shell_idxs;             /** Indices of the involved shells in the list of all shells. */

            /** Returns the angular momentum pair. */
            std::pair<int, int> getLPair() const
            {
                return {la, lb};
            }
        };

        /**
         * Constructs the shell datas for the auxiliary basis set, up to max L.
         */
        std::vector<ShellData> shellDataAux(const Structure &structure);

        /**
         *
         */
        std::vector<ShellPairData> shellPairData(bool use_symm, const Structure &structure);

        std::vector<ShellPairData> shellPairData(const std::vector<Shell> &shells_a,
                                                 const std::vector<Shell> &shells_b);
    }
}