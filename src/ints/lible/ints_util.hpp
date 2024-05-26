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
        std::vector<std::array<int, 3>> returnCartesianExps(const int l);

        /**
         *
         */
        std::vector<std::pair<int, int>> returnLPairs(const int l_max);
    }
}