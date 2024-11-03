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

            /** Angular momentum. */
            const int l;

            /** Number of shells. */
            const int n_shells;

            /** Total number of Gaussian primitives. */
            const int n_primitives;

            /** Contraction coefficients including the primitive norms. */
            const std::vector<double> coeffs;

            /** Coordinates of atoms corresponding to shells. */
            const std::vector<double> coords;

            /** Exponents of the Gaussian primitives. */
            const std::vector<double> exps;

            /** Normalization constants of the atomic orbitals in spherical basis. */
            const std::vector<double> norms;

            /** Contraction depths corresponding to each shell. */
            const std::vector<int> cdepths;

            /** Offsets of the contraction data (coeffs and exps) for each shell.*/
            const std::vector<int> coffsets;

            /** Offsets of the spherical Hermite expansion coefficients. */
            const std::vector<int> offsets_ecoeffs;
            
            /** Offsets of the shell atomic orbital norms. */
            const std::vector<int> offsets_norms;

            /** Offsets of the atomic orbital positions in the list of all atomic orbitals. */
            const std::vector<int> offsets_sph;
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

            /** Angular momentum of the left shell. */
            const int la;

            /** Angular momentum of the right shell. */
            const int lb;

            /** Total number of atoms in the system. */
            const int n_atoms;
            const int n_pairs;
            const int n_prim_pairs;

            /** Coordinates of each atom in the structure. */
            const std::vector<double> atomic_coords;

            const std::vector<double> coeffs;

            /** Coordinates of shells in the shell pairs. */
            const std::vector<double> coords;

            const std::vector<double> exps;
            const std::vector<double> norms;

            /** Atomic numbers of all the atoms in the system. */
            const std::vector<int> atomic_nrs;

            const std::vector<int> cdepths;
            const std::vector<int> coffsets;

            const std::vector<int> offsets_cart;
            const std::vector<int> offsets_ecoeffs;
            const std::vector<int> offsets_norms;
            const std::vector<int> offsets_sph;
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
    }
}