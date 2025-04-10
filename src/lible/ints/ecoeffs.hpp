#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <vector>

namespace lible
{
    namespace ints
    {
        /**
         * Calculates the Hermite expansion coefficients for the given l-pairs and shell-pair
         * datas. The function assumes the given shell-pair datas correspond to the l-pairs.
         * The E-coefficients in the ket are transposed.
         */
        std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
        ecoeffsFromSPDatas(const std::vector<std::pair<int, int>> &l_pairs,
                           const std::vector<ShellPairData> &sp_datas);

        /**
         * Calculates the Hermite expansion coefficients for the given l-value. It is assumed
         * that the shell datas are ordered as 0,...,l_max_aux. The function assumes shell data
         * for the auxiliary basis set.
         */
        std::vector<std::vector<double>>
        ecoeffsFromShellDatasAux(const int l_max_aux, const std::vector<ShellData> &sh_datas);

        /**
         * Calculates the Hermite expansion coefficients for the given l value. It is assumed 
         * that the shell datas are ordered as 0,...,l_max_aux. The function assume shell data
         * for the auxiliary basis set. The E-coefficients in the ket are transposed.
         */
        std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
        ecoeffsFromShellDatasAuxPairs(const int l_max_aux, const std::vector<ShellData> &sh_datas);

        /**
         * Calculates Hermite expansion coefficients for a diagonal shell pair located on the same
         * atom.
         */
        void ecoeffsShell(const int l, const std::vector<double> &exps,
                          std::vector<std::vector<double>> &ecoeffs_out);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, 0) | r = x,y,z } for each
         * pair of Gaussian primitives for each shell pair.
         */
        void ecoeffsShellPairs3D(const int la, const int lb, const ShellPairData &shell_pair_data,
                                 std::vector<std::vector<vec3d>> &ecoeffs_out);

        /**
         * Calculates the Hermite expansion coefficients { E(r, i, j, t) | r = x,y,z && t <= la + lb}
         * for each pair of Gaussian primitives for each shell pair.
         */
        void ecoeffsShellPairs4D(const int la, const int lb, const ShellPairData &sp_data,
                                 std::vector<std::vector<vec4d>> &ecoeffs_out);

        /**
         * Calculates the Hermite expanesion coefficients { E(r, i, j, t) | r = x,y,z && t <= t_max}
         * for each pair of Gaussian primitives for each shell pair.
         */
        void ecoeffsShellPairs4D(const int la, const int lb, const int t_max,
                                 const ShellPairData &sp_data,
                                 std::vector<std::vector<vec4d>> &ecoeffs_out);

        /**
         * Calculates the Cartesian-to-spherical-transformed Hermite expansion coefficients for
         * given shells.
         */
        void ecoeffsShellsSpherical(const int l, const ShellData &sh_data,
                                    std::vector<double> &ecoeffs_out);

        /**
         * Calculates the Cartesian-to-spherical-transformed Hermite expansion coefficients for
         * given shells. Calculates also the transposed coefficients
         */
        void ecoeffsShellsSpherical(const int l, const ShellData &sh_data,
                                    std::vector<double> &ecoeffs_out,
                                    std::vector<double> &ecoeffs_tsp_out);

        /**
         * Calculates the Cartesian-to-spherical-transformed Hermite expansion coefficients for
         * given shell pairs.
         */
        void ecoeffsSPsSpherical(const int la, const int lb, const ShellPairData &sp_data,
                                 std::vector<double> &ecoeffs_out);
                                 
        /**
         * Calculates the Cartesian-to-spherical-transformed Hermite expansion coefficients for
         * given shell pairs. Calculates also the transposed coefficients.
         */
        void ecoeffsSPsSpherical(const int la, const int lb, const ShellPairData &sp_data,
                                 std::vector<double> &ecoeffs_out,
                                 std::vector<double> &ecoeffs_tsp_out);

        /**
         * Calculates the Hermite expansion coefficients for a single primitive Gaussian along
         * the three Cartesian directions.
         */
        void ecoeffsPrimitive(const double a, const int l, vec2d &ecoeffs_x, vec2d &ecoeffs_y,
                              vec2d &ecoeffs_z);

        /**
         * Calculates the Hermite expansion coefficients for a pair of primitive Gaussians along
         * the three Cartesian directions.
         */
        void ecoeffsPrimitivePair(const double a, const double b, const int la, const int lb,
                                  const std::array<double, 3> &xyz_a,
                                  const std::array<double, 3> &xyz_b,
                                  const std::array<double, 3> &Kab, vec3d &ecoeffs_x,
                                  vec3d &ecoeffs_y, vec3d &ecoeffs_z);

        /**
         * Calculates the Hermite expansion coefficents of a Cartesian Gaussian along one
         * Cartesian direction.
         */
        void ecoeffsRecurrence1(const double one_o_2a, const int l, vec2d &ecoeffs);

        /**
         * Calculates the Hermite expansion coefficients of a Cartesian Gaussian product along
         * one Cartesian direction.
         */
        void ecoeffsRecurrence2(const double a, const double b, const double PA, const double PB,
                                const double one_o_2p, const int la, const int lb, vec3d &ecoeffs);
    }
}