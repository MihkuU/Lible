#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <vector>

#include <armadillo>

namespace lible
{
    namespace ints
    {
        /**
         *
         */
        void coeffs(const double a, const double b, const double PA, const double PB,
                    const double one_o_2p, const int la, const int lb, vec3d &E);

        /**
         *
         */
        void coeffs(const double a, const double b, const int la, const int lb,
                    const std::array<double, 3> &A, const std::array<double, 3> &B,
                    const std::array<double, 3> &Kab, vec3d &Ex, vec3d &Ey, vec3d &Ez);

        /**
         *
         */
        void calcECoeffs(const int l, const std::vector<double> &exps,
                         std::vector<arma::dmat> &ecoeffs_out);

        /**
         *
         */
        void calcECoeffs(const int la, const int lb, const ShellPairData &shell_pair_data,
                         std::vector<std::vector<arma::dmat>> &ecoeffs_out);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, 0) | r = x,y,z } for each
         * pair of Gaussian primitives for each shell pair.
         */
        void calcECoeffs(const int la, const int lb, const ShellPairData &shell_pair_data,
                         std::vector<std::vector<vec3d>> &ecoeffs_out);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, t) | r = x,y,z && t <= la + lb}
         * for each pair of Gaussian primitives for each shell pair.
         */
        void calcECoeffs(const int la, const int lb, const ShellPairData &shell_pair_data,
                         std::vector<std::vector<vec4d>> &ecoeffs_out);

        /**
         *
         */
        void calcECoeffsSpherical(const int la, const int lb,
                                  const ShellPairData &shell_pair_data,
                                  std::vector<std::vector<arma::dmat>> &ecoeffs_out);

        /**
         *
         */
        void calcECoeffsSpherical(const int la, const int lb,
                                  const ShellPairData &shell_pair_data,
                                  std::vector<double> &ecoeffs_out,
                                  std::vector<double> &ecoeffs_tsp_out);
    }
}