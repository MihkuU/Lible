#pragma once

#include <lible/shell_pair_data.hpp>
#include <lible/types.hpp>

#include <armadillo>

namespace lible
{
    namespace ints
    {
        namespace MD
        {
            /**
             *
             */
            struct IdxsTUV
            {
                int t, u, v;
            };

            struct IdxsCart
            {
                int i, j, k;
            };

            /**
             *
             */
            vec3i returnTUVPoss(const int l);

            /**
             *
             */
            std::vector<IdxsTUV> returnIdxsTUV(const int l);

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
            void calcRInts(const int la, const int lb, const double p,
                           const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                           vec4d &rints_tmp, vec3d &rints_out);

            /**
             *
             */
            void calcRInts(const int la, const int lb, const double p,
                           const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,                           
                           const std::vector<IdxsTUV> &tuv_idxs_a,
                           const std::vector<IdxsTUV> &tuv_idxs_b,
                           vec4d &rints_tmp, arma::dmat &rints_out);

            void coeffs(const double a, const double b, const double PA, const double PB,
                        const double one_o_2p, const int la, const int lb, vec3d &E);

            void coeffs(const double a, const double b, const int la, const int lb,
                        const std::array<double, 3> &A, const std::array<double, 3> &B,
                        const std::array<double, 3> &Kab, vec3d &Ex, vec3d &Ey, vec3d &Ez);

            // void coeffs(const double one_o_2p, const int la, const int lb,
            //             const arma::vec::fixed<3> &A, const arma::vec::fixed<3> &B,
            //             const std::array<double, 3> &Kab,
            //             vec4d &E);

            void coeffsSpherical(); // TODO: future

            void test(); // TODO: remove
        }
    }
}