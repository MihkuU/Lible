#pragma once

#define _USE_MATH_DEFINES

#include <array>
#include <cmath>
#include <utility>
#include <vector>

#include "ints_util.h"

namespace lible
{
    /*
     * Here we deal with atomic orbital shells. We follow conventions outlined in https://iodata.readthedocs.io/en/latest/basis.html.
     * Definitions:
     *   Shell - set of basis functions with the same angular momentum, contraction coefficients and
     *   exponents of gaussian primitives.
     *
     */
    // TODO: Shove it inside 'ints' namespace
    namespace shells
    {
        inline double calcPureGaussianPrimitiveNorm(const int angular_momentum,
                                                    const double &exponent)
        {
            double norm = std::sqrt(std::pow(2 * exponent / M_PI, 1.5) * std::pow(4 * exponent, angular_momentum) /
                                    IntsUtil::doubleFactorial(2 * angular_momentum - 1));
            return norm;
        }

        inline std::size_t calcShellDimCartesian(const int &angular_momentum)
        {
            return (angular_momentum + 1) * (angular_momentum + 2) / 2;
        }

        inline std::size_t calcShellDimSpherical(const int &angular_momentum)
        {
            return 2 * angular_momentum + 1;
        }

        inline std::vector<std::array<int, 3>> calcShellCartesianExps(const int &angular_momentum)
        {
            std::size_t dim_cartesian = calcShellDimCartesian(angular_momentum);

            std::vector<std::array<int, 3>> cartesian_exps(dim_cartesian);
            for (int x = angular_momentum, pos = 0; x >= 0; x--)
                for (int y = angular_momentum; y >= 0; y--)
                    for (int z = angular_momentum; z >= 0; z--)
                        if (x + y + z == angular_momentum)
                        {
                            cartesian_exps[pos] = {x, y, z};
                            pos++;
                        }

            return cartesian_exps;
        }

        std::vector<double> calcShellNormalization(const int &angmom,
                                                   const std::vector<double> &contraction_coeffs,
                                                   const std::vector<double> &contraction_exps,
                                                   const std::vector<std::array<int, 3>> &cartesian_exps);

        struct Shell
        {
            Shell(const int &angular_momentum,
                  const int &atomic_number,
                  const std::size_t &dim_cartesian,
                  const std::size_t &dim_spherical,
                  const std::size_t &pos,
                  const std::array<double, 3> &xyz_coordinates,
                  const std::vector<double> &contraction_coeffs,
                  const std::vector<double> &contraction_exps,
                  const std::vector<double> &normalization,
                  const std::vector<std::array<int, 3>> cartesian_exps)
                : angular_momentum(angular_momentum),
                  atomic_number(atomic_number),
                  dim_cartesian(dim_cartesian),
                  dim_spherical(dim_spherical),
                  pos(pos),
                  xyz_coordinates(xyz_coordinates),
                  contraction_coeffs(contraction_coeffs),
                  contraction_exps(contraction_exps),
                  normalization(normalization),
                  cartesian_exps(cartesian_exps)
            {
            }

            Shell(const int &angular_momentum,
                  const std::size_t &dim_cartesian,
                  const std::size_t &dim_spherical)
                : angular_momentum(angular_momentum),
                  atomic_number(0),
                  dim_cartesian(dim_cartesian),
                  dim_spherical(dim_spherical),
                  pos(0),
                  xyz_coordinates({}),
                  contraction_coeffs({}),
                  contraction_exps({}),
                  normalization({}),
                  cartesian_exps({})
            {
            }

            const int angular_momentum;
            const int atomic_number;
            const std::size_t dim_cartesian;
            const std::size_t dim_spherical;
            const std::size_t pos;
            const std::array<double, 3> xyz_coordinates;
            const std::vector<double> contraction_coeffs;
            const std::vector<double> contraction_exps;
            const std::vector<double> normalization;
            const std::vector<std::array<int, 3>> cartesian_exps;
        };

        struct ShellPair
        {
            //TODO: remove this object
            ShellPair(const Shell &first, const Shell &second) : first(first), second(second) {}

            const Shell first;
            const Shell second;
        };
    }
}