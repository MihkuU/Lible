#pragma once

#include <lible/ints/utils.hpp>

#include <array>
#include <vector>

#include <armadillo> // TODO: remove

namespace lible
{
    namespace ints
    {
        // TODO: use std::array<double, 3> instead of arma::vec::fixed<3>

        /** Calculates the Hermite Coulomb integrals as a 3D-array. */
        void calcRInts(const int la, const int lb, const double p,
                       const arma::vec::fixed<3> &xyz_ab, const std::vector<double> &fnx,
                       vec4d &rints_tmp, vec3d &rints_out);

        /** Calculates the Hermite Coulomb integrals as a flattened vector. */
        void calcRInts(const int la, const int lb, const double fac, const double p,
                       const arma::vec::fixed<3> &xyz_pq, const std::vector<double> &fnx,
                       const std::vector<std::array<int, 3>> &tuv_idxs_a,
                       const std::vector<std::array<int, 3>> &tuv_idxs_b,
                       vec4d &rints_tmp, std::vector<double> &rints_out);
    }
}
