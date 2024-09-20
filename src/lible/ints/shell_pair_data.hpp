#pragma once

#include <lible/ints/structure.hpp>

namespace lible
{
    namespace ints
    {
        /** */
        struct ShellData
        {
            ShellData(const int l, const int n_shells, const int n_primitives,
                      const std::vector<double> &coeffs, const std::vector<double> &exps,
                      const std::vector<double> &norms, const std::vector<int> &cdepths,
                      const std::vector<int> &coffsets, const std::vector<int> &offsets_ecoeffs,
                      const std::vector<int> &offsets_norms, const std::vector<int> &offsets_sph)
                : l(l), n_shells(n_shells), n_primitives(n_primitives),
                  coeffs(coeffs), exps(exps), norms(norms), cdepths(cdepths), coffsets(coffsets),
                  offsets_ecoeffs(offsets_ecoeffs), offsets_norms(offsets_norms),
                  offsets_sph(offsets_sph)
            {
            }

            const int l;
            const int n_shells;
            const int n_primitives;

            const std::vector<double> coeffs;
            const std::vector<double> coords;
            const std::vector<double> exps;
            const std::vector<double> norms;

            const std::vector<int> cdepths;
            const std::vector<int> coffsets;

            const std::vector<int> offsets_ecoeffs;
            const std::vector<int> offsets_norms;
            const std::vector<int> offsets_sph;
        };

        /** */
        struct ShellPairData
        {
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

            const int la;
            const int lb;
            const int n_atoms;
            const int n_pairs;
            const int n_prim_pairs;

            const std::vector<double> atomic_coords;

            const std::vector<double> coeffs;
            const std::vector<double> coords;
            const std::vector<double> exps;
            const std::vector<double> norms;

            const std::vector<int> atomic_nrs;

            const std::vector<int> cdepths;
            const std::vector<int> coffsets;

            const std::vector<int> offsets_cart;
            const std::vector<int> offsets_ecoeffs;
            const std::vector<int> offsets_norms;
            const std::vector<int> offsets_sph;
        };

        /** */
        ShellData constructShellDataAux(const int l, const Structure &structure);

        /** */
        ShellPairData constructShellPairData(const int la, const int lb,
                                             const Structure &structure);
    }
}