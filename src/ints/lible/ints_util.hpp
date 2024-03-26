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
        double calcPurePrimitiveNorm(int l, double exp);

        double doubleFactorial(int n);

        int dimCartesians(int l);

        int dimSphericals(int l);

        std::vector<std::array<int, 3>> returnCartesianExps(int l);

        std::string returnBasisPath(const std::string &basis_set);
    }
}