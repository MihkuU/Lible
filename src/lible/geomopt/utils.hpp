#pragma once

#include <lible/geomopt/geomopt.hpp>

#include <array>

namespace lible::geomopt
{
    using std::size_t;

    /// Test for whether the vectors are parallel.
    bool areParallel(const std::array<double, 3> &u, const std::array<double, 3> &v,
                     double tol = tolerance);

    /// Kronecker delta-fun.
    double delta(size_t i, size_t j);

    /// The zeta defined in (18) in https://doi.org/10.1063/1.1515483.
    double zeta(size_t a, size_t m, size_t n);

    /// Calculates the L2 norm of a 3-vector.
    double norm(const std::array<double, 3> &u);

    /// Calculates the dot product between two 3-vectors.
    double dot(const std::array<double, 3> &u, const std::array<double, 3> &v);

    /// Calculates the angle between two given vectors. Throws if one of the vectors has
    /// approximately zero length. The angle is returned in radians.
    double angle(const std::array<double, 3> &u, const std::array<double, 3> &v,
                 double tol = tolerance);

    /// Calculates the cross-product between two vectors.
    std::array<double, 3> cross(const std::array<double, 3> &u, const std::array<double, 3> &v);

    /// Calculates the difference between two 3-vectors.
    std::array<double, 3> operator-(const std::array<double, 3> &u, const std::array<double, 3> &v);

    /// Divides the 3-vector elements by the given value.
    std::array<double, 3> operator/(const std::array<double, 3> &u, double x);
}