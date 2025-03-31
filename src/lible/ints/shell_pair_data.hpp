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
            /** The constructor. */
            ShellData(const int l, const int n_shells, const int n_primitives,
                      const std::vector<double> &coeffs, const std::vector<double> &coords,
                      const std::vector<double> &exps, const std::vector<double> &norms,
                      const std::vector<int> &cdepths, const std::vector<int> &coffsets,
                      const std::vector<int> &offsets_ecoeffs, const std::vector<int> &offsets_norms,
                      const std::vector<int> &offsets_sph)
                : l(l), n_shells(n_shells), n_primitives(n_primitives),
                  coeffs(coeffs), coords(coords), exps(exps), norms(norms), cdepths(cdepths),
                  coffsets(coffsets), offsets_ecoeffs(offsets_ecoeffs),
                  offsets_norms(offsets_norms), offsets_sph(offsets_sph)
            {
            }

            int l;            /** Angular momentum. */
            int n_shells;     /** Number of shells. */
            int n_primitives; /** Total number of Gaussian primitives. */

            std::vector<double> coeffs; /** Contraction coefficients including the primitive norms. */
            std::vector<double> coords; /** Coordinates of atoms corresponding to shells. */
            std::vector<double> exps;   /** Exponents of the Gaussian primitives. */
            std::vector<double> norms;  /** Normalization constants of the atomic orbitals in spherical basis. */

            std::vector<int> cdepths;         /** Contraction depths corresponding to each shell. */
            std::vector<int> coffsets;        /** Offsets of the contraction data (coeffs and exps) for each shell.*/
            std::vector<int> offsets_ecoeffs; /** Offsets of the spherical Hermite expansion coefficients. */
            std::vector<int> offsets_norms;   /** Offsets of the shell atomic orbital norms. */
            std::vector<int> offsets_sph;     /** Offsets of the atomic orbital positions in the list of all atomic orbitals. */
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
            /** The constructor. */
            ShellPairData(const int la, const int lb, const int n_atoms, const int n_pairs,
                          const int n_prim_pairs, const std::vector<double> &atomic_coords,
                          const std::vector<double> &coeffs, const std::vector<double> &coords,
                          const std::vector<double> &exps, const std::vector<double> &norms,
                          const std::vector<int> &atomic_nrs, const std::vector<int> &cdepths,
                          const std::vector<int> &coffsets, const std::vector<int> &offsets_cart,
                          const std::vector<int> &offsets_ecoeffs,
                          const std::vector<int> &offsets_norms,
                          const std::vector<int> &offsets_sph)
                : la(la), lb(lb), n_atoms(n_atoms), n_pairs(n_pairs), n_prim_pairs(n_prim_pairs),
                  atomic_coords(atomic_coords), coeffs(coeffs), coords(coords), exps(exps),
                  norms(norms), atomic_nrs(atomic_nrs), cdepths(cdepths), coffsets(coffsets),
                  offsets_cart(offsets_cart), offsets_ecoeffs(offsets_ecoeffs),
                  offsets_norms(offsets_norms), offsets_sph(offsets_sph)
            {
            }

            int la; /** Angular momentum of the left shell. */
            int lb; /** Angular momentum of the right shell. */
            int n_atoms;
            int n_pairs; /** Total number of atoms in the system. */
            int n_prim_pairs;

            std::vector<double> atomic_coords; /** Coordinates of each atom in the structure. */
            std::vector<double> coeffs;
            std::vector<double> coords; /** Coordinates of shells in the shell pairs. */
            std::vector<double> exps;
            std::vector<double> norms;

            std::vector<int> atomic_nrs; /** Atomic numbers of all the atoms in the system. */
            std::vector<int> cdepths;
            std::vector<int> coffsets;
            std::vector<int> offsets_cart;
            std::vector<int> offsets_ecoeffs;
            std::vector<int> offsets_norms;
            std::vector<int> offsets_sph;
        };

        /**
         * \ingroup spdata
         *
         * Constructs the shell data corresponding to the auxilary basis set.
         */
        ShellData constructShellDataAux(const int l, const Structure &structure);

        /**
         * \ingroup spdata
         *
         *  Constructs the shell pair data corresponding to the main basis set.
         */
        ShellPairData constructShellPairData(const int la, const int lb,
                                             const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell datas for the auxiliary basis set, up to l_max.
         */
        std::vector<ShellData>
        constructShellDatasAux(const int l_max, const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell pair datas for the given l-pairs.
         */
        std::vector<ShellPairData>
        constructShellPairDatas(const std::vector<std::pair<int, int>> &l_pairs,
                                const Structure &structure);
    }
}