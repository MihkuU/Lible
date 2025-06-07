#pragma once

#include <lible/ints/structure.hpp>

namespace lible
{
    namespace ints
    {
        /** \defgroup spdata */

        // TODO: make variables private, add getters!

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
                      const std::vector<int> atomic_idxs, const std::vector<int> &cdepths, 
                      const std::vector<int> &coffsets, const std::vector<int> &offsets_ecoeffs, 
                      const std::vector<int> &offsets_norms, const std::vector<int> &offsets_sph)
                : l(l), n_shells(n_shells), n_primitives(n_primitives), coeffs(coeffs),
                  coords(coords), exps(exps), norms(norms),  atomic_idxs(atomic_idxs),
                  cdepths(cdepths), coffsets(coffsets), offsets_ecoeffs(offsets_ecoeffs),
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

            std::vector<int> atomic_idxs;     /** Indices of atoms involved in the shells. */
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
                          const std::vector<int> &atomic_idxs, const std::vector<int> &atomic_nrs,
                          const std::vector<int> &cdepths, const std::vector<int> &coffsets,
                          const std::vector<int> &offsets_cart,
                          const std::vector<int> &offsets_ecoeffs,
                          const std::vector<int> &offsets_ecoeffs_deriv1,
                          const std::vector<int> &offsets_norms,
                          const std::vector<int> &offsets_sph)
                : la(la), lb(lb), n_atoms(n_atoms), n_pairs(n_pairs), n_prim_pairs(n_prim_pairs),
                  atomic_coords(atomic_coords), coeffs(coeffs), coords(coords), exps(exps),
                  norms(norms), atomic_idxs(atomic_idxs), atomic_nrs(atomic_nrs), cdepths(cdepths),
                  coffsets(coffsets), offsets_cart(offsets_cart), offsets_ecoeffs(offsets_ecoeffs),
                  offsets_ecoeffs_deriv1(offsets_ecoeffs_deriv1), offsets_norms(offsets_norms), 
                  offsets_sph(offsets_sph)
            {
            }

            int la;           /** Angular momentum of the left shell. */
            int lb;           /** Angular momentum of the right shell. */
            int n_atoms;      /** Total numer of atoms in the system. */
            int n_pairs;      /** Total number of shell pairs. */
            int n_prim_pairs; /** Total number of Gaussian primitive pairs. */

            std::vector<double> atomic_coords; /** Coordinates of each atom in the structure. */
            std::vector<double> coeffs;        /** Contraction coefficients including the primitive norms. */
            std::vector<double> coords;        /** Coordinates of shells in the shell pairs. */
            std::vector<double> exps;          /** Exponents of the Gaussian primitives. */
            std::vector<double> norms;         /** Normalization consts of the atomic orbitals in spherical basis. */

            std::vector<int> atomic_idxs;            /** Indices of atoms involved in the shell pairs. */
            std::vector<int> atomic_nrs;             /** Atomic numbers of all the atoms in the system. */
            std::vector<int> cdepths;                /** Contraction depths corresponding to each shell in a shell pair. */
            std::vector<int> coffsets;               /** Offsets of the contraction data (coeffs and exps) for each shell in a shell pair.*/
            std::vector<int> offsets_cart;           /** Offsets of the atomic orbital (cartesian) positions in the list of all atomic orbitals. */
            std::vector<int> offsets_ecoeffs;        /** Offsets of the spherical Hermite expansion coefficients. */
            std::vector<int> offsets_ecoeffs_deriv1; /** Offsets of the 1st derivative spherical Hermite expansion coefficients. */
            std::vector<int> offsets_norms;          /** Offsets of the shell atomic orbital norms. */
            std::vector<int> offsets_sph;            /** Offsets of the atomic orbital (spherical) positions in the list of all atomic orbitals. */
        };

        /**
         * \ingroup spdata
         *
         * Constructs the shell data corresponding to the auxilary basis set.
         */
        ShellData shellDataAux(const int l, const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell pair data corresponding to the main basis set. It is assumed
         * that the integrals involving shells A and B are symmetric w.r.t. interchanging A and B.
         * If la != lb, shell pair data is created for pairs (ishellA, ishellB). When la == lb,
         * the data is create for (ishellA, ishellA') pairs such that ishellA >= ishellA'.
         */
        ShellPairData shellPairDataSymm(const int la, const int lb, const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell pair data corresponding to the main basis set. No symmetries are
         * being used in this version.
         */
        ShellPairData shellPairDataNoSymm(const int la, const int lb, const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell pair data corresponding to the main basis set. If symmetry is
         * enabled the data is created for (ishellA, ishellA') pairs such that ishellA >= ishellA'
         * when la == lb. If la != lb, shell pair data is created for pairs (ishellA, ishellB)
         * irrespective of whether symmetry is used or not.
         */
        ShellPairData constructShellPairData(const bool use_symm, const int la, const int lb,
                                             const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell datas for the auxiliary basis set, up to l_max.
         */
        std::vector<ShellData>
        shellDatasAux(const int l_max, const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell pair datas for the given l-pairs. The shell pair datas for each
         * (la, lb)-pair assume symmetries.
         */
        std::vector<ShellPairData>
        shellPairDatasSymm(const std::vector<std::pair<int, int>> &l_pairs,
                           const Structure &structure);

        /**
         * \ingroup spdata
         *
         * Constructs the shell pair datas for the given l-pairs. The shell pair datas for each
         * (la, lb)-pair assume symmetries.
         */
        std::vector<ShellPairData>
        shellPairDatasNoSymm(const std::vector<std::pair<int, int>> &l_pairs,
                                      const Structure &structure);
    }
}