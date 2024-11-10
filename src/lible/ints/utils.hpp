#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

#include <lible/types.hpp>

namespace lible
{
    namespace ints
    {
        /**
         *
         */
        double calcPurePrimitiveNorm(const int l, const double exp);

        /**
         *
         */
        double doubleFactorial(const int n);

        /** */
        consteval static int numSphericalsC(int l)
        {
            return 2 * l + 1;
        }

        /** */
        consteval static int numHermitesC(int l)
        {
            return (l + 1) * (l + 2) * (l + 3) / 6;
        }        

        /** Calculaces the index of a Cartesian Gaussian with given Cartesian exponents i, j, k. */
        constexpr int indexCart(const int i, const int j, const int k)
        {
            int jk = j + k;
            return jk * (jk + 1) / 2 + k;
        }

        /** */
        constexpr int numCartesians(const int l)
        {
            return (l + 1) * (l + 2) / 2;
        }

        /** */
        constexpr int numSphericals(const int l)
        {
            return 2 * l + 1;
        }

        /** */
        constexpr int numHermites(const int l)
        {
            return (l + 1) * (l + 2) * (l + 3) / 6;
        }

        /** */
        constexpr int numHermitesSum(const int l)
        {
            int sum = 0;
            for (int n = 0; n <= l; n++)
                sum += numHermites(n);

            return sum;
        }

        /**
         *
         */
        std::vector<std::array<int, 3>> cartExps(const int l);

        /**
         * Returns a list of angular momentum pairs {(0, 0), (1, 0), (1, 1), ..., (l_max, l_max)}.
         */
        std::vector<std::pair<int, int>> returnLPairs(const int l_max);

        /** */
        std::map<int, std::vector<std::pair<int, int>>> returnLPairsMap(const int l_max);

        /**
         *
         */
        vec3i returnHermiteGaussianPositions(const int l);

        /**
         *
         */
        std::vector<std::array<int, 3>> returnHermiteGaussianIdxs(const int l); 
        
    }
}