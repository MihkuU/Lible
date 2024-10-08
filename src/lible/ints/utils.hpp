#pragma once

#include <array>
#include <cassert>
#include <cmath>
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

        /**
         *
         */
        int dimCartesians(const int l);

        /**
         *
         */
        int dimSphericals(const int l);

        /**
         *
         */
        int dimHermiteGaussians(const int l);

        /**
         *
         */
        std::vector<std::array<int, 3>> cartExps(const int l);

        /**
         * Returns a list of angular momentum pairs {(0, 0), (1, 0), (1, 1), ..., (l_max, l_max)}.
         */
        std::vector<std::pair<int, int>> returnLPairs(const int l_max);

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