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
                                            const std::vector<std::array<int, 3>> &tuv_idxs_a,
                                            const std::vector<std::array<int, 3>> &tuv_idxs_b);

        /** Calculates the Hermite Coulomb integrals as a 3D array R(t, u, v). */
        vec3d calcRInts3D(const int l, const double p, const double *xyz_ab, const double *fnx);

        /** TODO: dox */
        std::vector<double> calcRInts_ERI4_Deriv1(const int l, const double fac, const double p,
                                                  const double *xyz_pq, const double *fnx,
                                                  const std::vector<std::array<int, 3>> &tuv_idxs_a,
                                                  const std::vector<std::array<int, 3>> &tuv_idxs_b);


        /** TODO: dox */
        std::vector<double> calcRInts_ERI2_deriv1(const int l, const double fac, const double p,
                                                  const double *xyz_ab, const double *fnx,
                                                  const std::vector<std::array<int, 3>> &tuv_idxs_a,
                                                  const std::vector<std::array<int, 3>> &tuv_idxs_b);                                                  
    }
}