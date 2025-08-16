#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <vector>

namespace lible
{
    namespace ints
    {
        /** */
        vec2d ecoeffsRecurrence1(const double one_o_2a, const int l);

        /** */
        vec3d ecoeffsRecurrence2(const double a, const double b, const int la, const int lb,
                                 const double PA, const double PB, const double Kab);

        /** */
        vec3d ecoeffsRecurrence2_n1(const double a, const double b, const int la, const int lb,
                                    const double A, const double B, const vec3d &ecoeffs0);

        /** */
        std::array<lible::vec2d, 3> ecoeffsPrimitive(const double a, const int l);

        /** */
        std::array<vec3d, 3> ecoeffsPrimitivePair(const double a, const double b, const int la,
                                                  const int lb, const double *xyz_a,
                                                  const double *xyz_b);

        /** */
        std::array<vec3d, 3> ecoeffsPrimitivePair_n1(const double a, const double b, const int la,
                                                     const int lb, const double *xyz_a,
                                                     const double *xyz_b,
                                                     const std::array<lible::vec3d, 3> &ecoeffs0);
        /** */
        std::vector<std::vector<double>> ecoeffsShell(const int l, const std::vector<double> &exps);

        /** */
        std::vector<double> ecoeffsSHARK(const ShellPairData &sp_data, const bool transpose = false);

        /** */
        std::vector<double> ecoeffsSHARK(const ShellData &sh_data, const bool transpose = false);

        /** */
        std::vector<double> ecoeffsD1SHARK(const ShellPairData &sp_data, const bool transpose = false);

        // SHARK E-coeffs contain norms and contraction coeffs TODO: mention that.

        // /** */
        // vec3d ecoeffsRecurrence2_n2(const double a, const double b, const int la, const int lb,
        //                             const double A, const double B, const vec3d &ecoeffs1);
                
        /** */
        // std::array<vec3d, 3> ecoeffsPrimitivePair_n2(const double a, const double b, const int la,
        //                                              const int lb, const double *xyz_a,
        //                                              const double *xyz_b,
        //                                              const std::array<lible::vec3d, 3> &ecoeffs1);

        // /** */
        // std::vector<double> ecoeffsD2SHARK(const ShellPairData &sp_data, const bool transpose = false);
    }
}