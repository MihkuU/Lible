#pragma once

#include <lible/types.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <vector>

namespace lible::ints
{
    /// Calculates the Hermite expansion coefficients for a single primitive Gaussian function
    /// in one Cartesian direction. Based on eq. (8) from https://doi.org/10.1063/5.0217001.
    /// Returned as E^i_t -> (i, t).
    vec2d ecoeffsRecurrence1(double one_o_2a, int l);

    /// Calculates the Hermite expansion coefficients for a primitive Gaussian function product in
    /// one Cartesian direction. Based on eqs. (9.5.6) and (9.5.7) from
    /// https://doi.org/10.1002/9781119019572. Returned as E^{ij}_t -> (i, j, t).
    vec3d ecoeffsRecurrence2(double a, double b, int la, int lb, double PA, double PB,
                             double Kab);

    /// Calculates the first derivative of the Hermite expansion coefficients. Based on eq. (20)
    /// from https://doi.org/10.1007/BF01132826.
    vec3d ecoeffsRecurrence2_n1(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0);

    /// Calculates the second derivative of the Hermite expansion coefficients. Based on eq. (23)
    /// from https://doi.org/10.1007/BF01132826
    vec3d ecoeffsRecurrence2_n2(double a, double b, int la, int lb, double A, double B,
                                const vec3d &ecoeffs0, const vec3d &ecoeffs1);

    /// Calculates the Hermite expansion coefficients for a single primitive Gaussian function
    /// in the three Cartesian directions.
    std::array<vec2d, 3> ecoeffsPrimitive(double a, int l);

    /// Calculates the Hermite expansion coefficients for a single primitive Gaussian function
    /// product in three Cartesian directions.
    std::array<vec3d, 3> ecoeffsPrimitivePair(double a, double b, int la, int lb,
                                              const double *xyz_a, const double *xyz_b);

    /// Calculates the first derivative of the Hermite expansion coefficients for a single
    /// primitive Gaussian function product in three Cartesian directions.
    std::array<vec3d, 3> ecoeffsPrimitivePair_n1(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0);

    /// Calculates the second derivative of the Hermite expansion coefficients for a single
    /// primitive Gaussian function product in three Cartesian directions.
    std::array<vec3d, 3> ecoeffsPrimitivePair_n2(double a, double b, int la, int lb,
                                                 const double *xyz_a, const double *xyz_b,
                                                 const std::array<vec3d, 3> &ecoeffs0,
                                                 const std::array<vec3d, 3> &ecoeffs1);

    /// Calculates the Hermite expansion coefficients for a single shell with given primitive
    /// Gaussian exponents. Returns E^{ii'}_0 E^{jj'}_0 E^{kk'}_0 for each primitive pair.
    std::vector<std::vector<double>> ecoeffsShell(int l, const std::vector<double> &exps);

    /// Calculates the Hermite expansion coefficients, E_{munu, tuv}, for each primitive Gaussian
    /// pair in each shell pair, where \mu \nu refer to the AO indices. The coefficients have AO
    /// norms and primitive contraction coefficients multiplied into.
    std::vector<double> ecoeffsSHARK(const ShellPairData &sp_data, bool transpose = false);

    /// Calculates the Hermite expansion coefficients, E_{mu, tuv}, for each primitive Gaussian
    /// pair in each shell pair. The coefficients have AO norms and primitive contraction
    /// coefficients multiplied into.
    std::vector<double> ecoeffsSHARK(const ShellData &sh_data, bool transpose = false);

    /// Calculates the first derivate Hermite expansion coefficients,
    /// The coefficients have AO norms and primitive contraction coefficients multiplied into.
    std::vector<double> ecoeffsD1SHARK(const ShellPairData &sp_data, bool transpose = false);
}
