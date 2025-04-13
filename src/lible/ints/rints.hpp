#pragma once

#include <lible/ints/utils.hpp>

#include <array>
#include <vector>

namespace lible
{
    namespace ints
    {
        /** */
        vec3d calcRInts(const int l, const double p, const double *xyz_ab, const double *fnx);

        /** Calculates the Hermite Coulomb integrals as a 3D-array. */
        void calcRInts_(const int la, const int lb, const double p,
                        const std::array<double, 3> &xyz_ab, const std::vector<double> &fnx,
                        vec4d &rints_tmp, vec3d &rints_out);

        /** Calculates the Hermite Coulomb integrals as a flattened vector. */
        void calcRInts_(const int la, const int lb, const double fac, const double p,
                        const std::array<double, 3> &xyz_pq, const std::vector<double> &fnx,
                        const std::vector<std::array<int, 3>> &tuv_idxs_a,
                        const std::vector<std::array<int, 3>> &tuv_idxs_b,
                        vec4d &rints_tmp, std::vector<double> &rints_out);
    }
}
