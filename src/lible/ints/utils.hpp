#pragma once

#include <array>
#include <map>
#include <vector>

#include <lible/types.hpp>

namespace lible::ints
{
    /// Implements eq. (6.5.10) from DOI:10.1002/9781119019572 without the odd negative values.
    double doubleFactorial(int n);

    /// Calculates the number of spherical Gaussians. Compile time only.
    consteval static int numSphericalsC(const int l)
    {
        return 2 * l + 1;
    }

    /// Calculates the number of Cartesian Gaussians. Compile time only.
    consteval static int numHermitesC(const int l)
    {
        return (l + 1) * (l + 2) * (l + 3) / 6;
    }

    /// Calculates the index (i, j, k) -> (ijk). Adopted from https://github.com/ValeevGroup/libintx.
    constexpr int indexCart(const int i, const int j, const int k) // TODO: remove i lol.
    {
        int jk = j + k;
        return jk * (jk + 1) / 2 + k;
    }

    /// Calculates the number of spherical Gaussians.
    constexpr int numSphericals(const int l)
    {
        return 2 * l + 1;
    }

    /// Calculates the number of Cartesian Gaussians.
    constexpr int numCartesians(const int l)
    {
        return (l + 1) * (l + 2) / 2;
    }

    /// Calculates the number of Hermite Gaussian triplets for 0,...,l.
    constexpr int numHermites(const int l)
    {
        return (l + 1) * (l + 2) * (l + 3) / 6;
    }

    /// Calculates the total number of Hermite Gaussian triplets for 0,...,l.
    constexpr int numHermitesSum(const int l)
    {
        int sum = 0;
        for (int n = 0; n <= l; n++)
            sum += numHermites(n);

        return sum;
    }

    /// Calculates a map {l, {{0, 0}, {1, 0}, ..., {l, l}}} such that la >= lb.
    std::map<int, std::vector<std::pair<int, int>>> getLPairsMap(int l_max);

    /// Calculates the positions for Hermite Gaussians (t, u, v) -> tuv.
    vec3i getHermiteGaussianPositions(int l);

    /// Calculates the positions for Hermite Gaussians (t, u, v) -> tuv.
    std::vector<std::array<int, 3>> getHermiteGaussianIdxs(int l);
}