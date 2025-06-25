#pragma once

#include <lible/ints/utils.hpp>

#include <array>
#include <vector>

namespace lible
{
    namespace ints
    {
        /**
         * Calculates the Hermite Coulomb integrals as a flattened matrix,
         *   R(t + t', u + u', v + v') -> R(tuv, t'u'v').
         * This is used in the SHARK method.
         */
        std::vector<double> calcRIntsMatrix(const int l, const double fac, const double p,
                                            const double *xyz_pq, const double *fnx,
                                            const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                            const std::vector<std::array<int, 3>> &hermite_idxs_b);

        /** Calculates the Hermite Coulomb integrals as a 3D array R(t, u, v). */
        vec3d calcRInts3D(const int l, const double p, const double *xyz_ab, const double *fnx);

        std::vector<double> calcRInts_ERI2D1(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_ab,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_b);

        std::vector<double> calcRInts_ERI2D2(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_ab,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_a,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_b);

        std::vector<double> calcRInts_ERI3D1(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_pc,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_bra,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_ket);

        std::vector<double> calcRInts_ERI3D2(const int l, const double alpha, const double fac,
                                             const double *fnx, const double *xyz_pc,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_bra,
                                             const std::vector<std::array<int, 3>> &hermite_idxs_ket);
    }
}