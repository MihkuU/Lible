#pragma once

#include <lible/ints/structure.hpp>

namespace lible
{
    namespace ints
    {
        // TODO: redesign the shelldata so that everything is flattened.
        // For example, the contraction coefficients for each shell-pair are aligned
        // side-by-side. This will require addition of offset/dimension vectors.
        struct ShellPairData
        {
            ShellPairData();

            ShellPairData(const int la, const int lb, const Structure &structure);

            const Structure *structure;

            int la;
            int lb;
            size_t n_pairs;
            size_t n_prim_pairs;

            std::vector<size_t> offsets_ecoeffs;
            std::vector<size_t> offsets_prims;
            std::vector<std::pair<size_t, size_t>> offsets;
            std::vector<std::pair<size_t, size_t>> offsets_cart;
            std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> coords;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> ccoeffs;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> exps;
            std::vector<std::pair<std::vector<double>, std::vector<double>>> norms;
        };

        /** */
        struct ShellPairData_new
        {
            ShellPairData_new(const int la, const int lb, const int n_atoms, const int n_pairs,
                              const int n_prim_pairs, const std::vector<double> &atomic_coords,
                              const std::vector<double> &coeffs, const std::vector<double> &coords,
                              const std::vector<double> &exps, const std::vector<double> &norms,
                              const std::vector<int> &atomic_nrs, const std::vector<int> &cdepths,
                              const std::vector<int> &coffsets,
                              const std::vector<int> &offsets_cart,
                              const std::vector<int> &offsets_ecoeffs,
                              const std::vector<int> &offsets_norms,
                              const std::vector<int> &offsets_primitives,
                              const std::vector<int> &offsets_sph)
                : la(la), lb(lb), n_atoms(n_atoms), n_pairs(n_pairs), n_prim_pairs(n_prim_pairs),
                  atomic_coords(atomic_coords), coeffs(coeffs), coords(coords), exps(exps),
                  norms(norms), atomic_nrs(atomic_nrs), cdepths(cdepths), coffsets(coffsets),
                  offsets_cart(offsets_cart), offsets_ecoeffs(offsets_ecoeffs),
                  offsets_norms(offsets_norms), offsets_primitives(offsets_primitives),
                  offsets_sph(offsets_sph)
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
            const std::vector<int> offsets_primitives;
            const std::vector<int> offsets_sph;
        };

        ShellPairData_new constructShellPairData(const int la, const int lb,
                                                 const Structure &structure);
    }
}