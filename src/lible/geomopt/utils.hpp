#pragma once

#include <lible/types.hpp>

#include <array>

namespace lible::geomopt
{
    /// Tolerance used for judging various things: vector zero length, vector parallelity, etc.
    constexpr double tolerance = 1e-10;

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

    /// Structure representing a hyper-dual number. Based on https://doi.org/10.3390/math13243909.
    struct HyperDual
    {
        /// Constructor for typical second derivative calculations.
        explicit HyperDual(const double x0) : x0_(x0)
        {
        }

        /// Generic constructor.
        HyperDual(const double x0, const double x1, const double x2, const double x12)
            : x0_(x0), x1_(x1), x2_(x2), x12_(x12)
        {
        }

        double x0_{};
        double x1_{};
        double x2_{};
        double x12_{};
    };

    /// Adds two hyper-dual numbers based on eq. (3) from https://doi.org/10.3390/math13243909.
    HyperDual operator+(const HyperDual &a, const HyperDual &b);

    /// Subtracts two hyper-dual numbers.
    HyperDual operator-(const HyperDual &a, const HyperDual &b);

    /// Multiplies two hyper-dual numbers based on eq. (4) from https://doi.org/10.3390/math13243909.
    HyperDual operator*(const HyperDual &a, const HyperDual &b);

    /// Scales a hyper-dual number.
    HyperDual operator*(double f, const HyperDual &a);

    /// Divides two hyper-dual numbers based on eq. (5) from https://doi.org/10.3390/math13243909.
    HyperDual operator/(const HyperDual &a, const HyperDual &b);

    /// Calculates a hyper-dual number inverse using eq. (5) from https://doi.org/10.3390/math13243909.
    HyperDual inv(const HyperDual &a);

    /// Calculates a square root of a hyper-dual number.
    HyperDual sqrt(const HyperDual &a);

    /// Calculates a sine function of a hyper-dual number.
    HyperDual sin(const HyperDual &a);

    /// Calculates the inverse cosine of a hyper-dual number.
    HyperDual acos(const HyperDual &a);

    /// Type alias for a 3-element array of hyper-dual numbers.
    using hd3_t = std::array<HyperDual, 3>;

    /// Calculates the norm of a 3-array of hyper-dual numbers.
    HyperDual norm(const hd3_t &u);

    /// Calculates the product of a 3-array of hyper-dual numbers.
    HyperDual dot(const hd3_t &u, const hd3_t &v);

    /// Calculates the cross-product of two 3-arrays of hyper-dual numbers.
    hd3_t cross(const hd3_t &u, const hd3_t &v);

    /// Subtracts two 3-arrays of hyper-dual numbers.
    hd3_t operator-(const hd3_t &u, const hd3_t &v);

    /// Divides the 3-array of hyper-dual numbers with another hyper-dual number, `x`.
    hd3_t operator/(const hd3_t &u, const HyperDual &x);

    /// Calculates the bond length with hyper-dual numbers. Based on FIG 1 from
    /// https://doi.org/10.1063/1.1515483.
    HyperDual bondLengthHD(const hd3_t &xyz_m, const hd3_t &xyz_n);

    /// Calculates the bond angle with hyper-dual numbers. Using eq. (23) from
    /// https://doi.org/10.1063/1.1515483.
    HyperDual bondAngleHD(const hd3_t &xyz_m, const hd3_t &xyz_o, const hd3_t &xyz_n);

    /// Calculates the dihedral angle with hyper-dual numbers. Using eq. (31) from
    /// https://doi.org/10.1063/1.1515483.
    HyperDual dihedralAngleHD(const hd3_t &xyz_m, const hd3_t &xyz_o, const hd3_t &xyz_p,
                              const hd3_t &xyz_n);

    /// Calculates the first derivatives of a bond length using hyper-dual numbers.
    std::array<double, 6> bondLengthGradientHD(const std::array<double, 3> &xyz_m,
                                               const std::array<double, 3> &xyz_n);

    /// Calculates the first derivatives of a bond angle using hyper-dual numbers.
    std::array<double, 9> bondAngleGradientHD(const std::array<double, 3> &xyz_m,
                                              const std::array<double, 3> &xyz_o,
                                              const std::array<double, 3> &xyz_n);

    /// Calculates the first derivatives of a dihedral angle using hyper-dual numbers.
    std::array<double, 12> dihedralAngleGradientHD(const std::array<double, 3> &xyz_m,
                                                   const std::array<double, 3> &xyz_o,
                                                   const std::array<double, 3> &xyz_p,
                                                   const std::array<double, 3> &xyz_n);

    /// Calculates the second derivatives of a bond length using hyper-dual numbers.
    arr2d<double, 6, 6> bondLengthHessianHD(const std::array<double, 3> &xyz_m,
                                            const std::array<double, 3> &xyz_n);

    /// Calculates the second derivatives of a bond angle using hyper-dual numbers.
    arr2d<double, 9, 9> bondAngleHessianHD(const std::array<double, 3> &xyz_m,
                                           const std::array<double, 3> &xyz_o,
                                           const std::array<double, 3> &xyz_n);

    /// Calculates the second derivatives of a dihedral angle using hyper-dual numbers.
    arr2d<double, 12, 12> dihedralAngleHessianHD(const std::array<double, 3> &xyz_m,
                                                 const std::array<double, 3> &xyz_o,
                                                 const std::array<double, 3> &xyz_p,
                                                 const std::array<double, 3> &xyz_n);
}
