#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/ints_util.hpp>

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

            /** Calculates the Hermite expansion coefficients in a column-major order. */
            void calcECoeffsSphericalCM(const int la, const int lb,
                                        const ShellPairData &shell_pair_data,
                                        std::vector<double> &ecoeffs_out,
                                        std::vector<double> &ecoeffs_tsp_out);

            /** Calculates the Hermite expansion coefficients ina row-major order. */
            void calcECoeffsSphericalRM(const int la, const int lb,
                                        const ShellPairData &shell_pair_data,
                                        std::vector<double> &ecoeffs_out,
                                        std::vector<double> &ecoeffs_tsp_out);

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

            /**
             *
             */
            void calcRInts(const int la, const int lb, const double p,
                           const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                           const std::vector<IdxsTUV> &tuv_idxs_a,
                           const std::vector<IdxsTUV> &tuv_idxs_b,
                           vec4d &rints_tmp, std::vector<double> &rints_out);

            /** Calculates the Hermite Coulomb integrals in a column-major order. */
            void calcRIntsCM(const int la, const int lb, const double p, const double fac,
                             const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                             const std::vector<IdxsTUV> &tuv_idxs_a,
                             const std::vector<IdxsTUV> &tuv_idxs_b,
                             vec4d &rints_tmp, std::vector<double> &rints_out);

            /** Calculates the Hermite Coulomb integrals in a row-major order. */
            void calcRIntsRM(const int la, const int lb, const double p, const double fac,
                             const arma::vec::fixed<3> &RPC, const std::vector<double> &fnx,
                             const std::vector<IdxsTUV> &tuv_idxs_a,
                             const std::vector<IdxsTUV> &tuv_idxs_b,
                             vec4d &rints_tmp, std::vector<double> &rints_out);

            void coeffs(const double a, const double b, const double PA, const double PB,
                        const double one_o_2p, const int la, const int lb, vec3d &E);

            void coeffs(const double a, const double b, const int la, const int lb,
                        const std::array<double, 3> &A, const std::array<double, 3> &B,
                        const std::array<double, 3> &Kab, vec3d &Ex, vec3d &Ey, vec3d &Ez);
        }
    }
}