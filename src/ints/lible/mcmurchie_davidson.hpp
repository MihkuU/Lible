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
            void calcHCoeffs(const int &l, const std::vector<double> &exps,
                             std::vector<arma::dmat> &h_coeffs_out);

            void calcHCoeffs(const int &la, const int &lb, const ShellPairData &shell_pair_data,
                             std::vector<std::vector<arma::dmat>> &h_coeffs_out);

            void coeffs(const double &a, const double &b, const double &PA, const double &PB,
                        const double one_o_2p, const int &la, const int &lb, vec3d &E);

            void coeffs(const double &a, const double &b, const int &la, const int &lb,
                        const std::array<double, 3> &A, const std::array<double, 3> &B,
                        const std::array<double, 3> &Kab,
                        vec3d &E_x, vec3d &E_y, vec3d &E_z);

            void coeffsSpherical(); // TODO: future

            void test(); // TODO: remove
        }
    }
}