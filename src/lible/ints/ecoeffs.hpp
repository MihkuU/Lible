#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <vector>

namespace lible::ints
{
    /** */
    vec2d ecoeffsRecurrence1(double one_o_2a, int l);

    /** */
    vec3d ecoeffsRecurrence2(double a, double b, int la, int lb, double PA, double PB,
                             double Kab);

    /** */
    vec3d ecoeffsRecurrence2_n1(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0);

    /** */
    vec3d ecoeffsRecurrence2_n2(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0, const vec3d &ecoeffs1);

    /** */
    std::array<lible::vec2d, 3> ecoeffsPrimitive(double a, int l);

    /** */
    std::array<vec3d, 3> ecoeffsPrimitivePair(double a, double b, int la, int lb,
                                              const double *xyz_a, const double *xyz_b);

    /** */
    std::array<vec3d, 3> ecoeffsPrimitivePair_n1(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0);

    /** */
    std::array<vec3d, 3> ecoeffsPrimitivePair_n2(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0,
                                                 const std::array<vec3d, 3> &ecoeffs1);

    /** */
    std::vector<std::vector<double>> ecoeffsShell(int l, const std::vector<double> &exps);

    /** */
    std::vector<double> ecoeffsSHARK(const ShellPairData &sp_data, bool transpose = false);

    /** */
    std::vector<double> ecoeffsSHARK(const ShellData &sh_data, bool transpose = false);

    /** */
    std::vector<double> ecoeffsD1SHARK(const ShellPairData &sp_data, bool transpose = false);

    // SHARK E-coeffs contain norms and contraction coeffs TODO: mention that.

    // /** */
    // vec3d ecoeffsRecurrence2_n2(const double a, const double b, const int la, const int lb,
    //                             const double A, const double B, const vec3d &ecoeffs1);



    // /** */
    // std::vector<double> ecoeffsD2SHARK(const ShellPairData &sp_data, const bool transpose = false);
}
