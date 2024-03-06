#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace lible
{
    namespace ints
    {
        inline constexpr double doubleFactorial(const int &n)
        {
            /*
             * Definition from DOI:10.1002/9781119019572 in eq. (6.5.10) without the odd negative values.
             */

            double double_factorial = 1;
            if (n == 0)
                return double_factorial;
            else if (n == -1)
                return double_factorial;
            else if (n > 0 and (n % 2 == 0))
                for (int i = n; i >= 2; i -= 2)
                    double_factorial *= i;
            else if (n > 0 and (n % 2 == 1))
                for (int i = n; i >= 1; i -= 2)
                    double_factorial *= i;
            else
                throw std::runtime_error("False input value for the double factorial!");

            return double_factorial;
        }

        inline constexpr double calcPurePrimitiveNorm(const int angmom,
                                                      const double &exp)
        {
            double partial_norm = std::sqrt(std::pow(2 * exp / M_PI, 1.5) * std::pow(4 * exp, angmom) /
                                            doubleFactorial(2 * angmom - 1));

            return partial_norm;
        }

        inline constexpr int dimCartesians(const int &l)
        {
            return (l + 1) * (l + 2) / 2;
        }

        inline constexpr int dimSphericals(const int &l)
        {
            return 2 * l + 1;
        }

        inline std::vector<std::array<int, 3>> returnCartesianExps(const int &angmom)
        {
            std::size_t dim_cart = dimCartesians(angmom);

            std::vector<std::array<int, 3>> cartesian_exps(dim_cart);
            for (int x = angmom, pos = 0; x >= 0; x--)
                for (int y = angmom - x; y >= 0; y--)
                {
                    int z = angmom - x - y;
                    cartesian_exps[pos] = {x, y, z};
                    pos++;
                }

            return cartesian_exps;
        }

        std::string returnBasisPath(const std::string &basis_set);       
    }
}