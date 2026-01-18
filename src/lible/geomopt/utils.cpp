#include <lible/geomopt/utils.hpp>

#include <cmath>
#include <stdexcept>

namespace lgopt = lible::geomopt;

double lgopt::norm(const std::array<double, 3> &u)
{
    const auto &[x, y, z] = u;

    return std::sqrt(x * x + y * y + z * z);
}

double lgopt::dot(const std::array<double, 3> &u, const std::array<double, 3> &v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

double lgopt::angle(const std::array<double, 3> &u, const std::array<double, 3> &v,
                    const double tol)
{
    double norm_u = norm(u);
    double norm_v = norm(v);

    if (norm_u < tol || norm_v < tol)
        throw std::runtime_error("angle(): angle between vectors with zero length.");

    double cos_angle = dot(u, v) / (norm_u * norm_v);

    return std::acos(cos_angle);
}

std::array<double, 3> lgopt::cross(const std::array<double, 3> &u, const std::array<double, 3> &v)
{
    return {{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]}};
}

std::array<double, 3> lgopt::operator-(const std::array<double, 3> &u,
                                       const std::array<double, 3> &v)
{
    return {{u[0] - v[0], u[1] - v[1], u[2] - v[2]}};
}

std::array<double, 3> lgopt::operator/(const std::array<double, 3> &u, const double x)
{
    return {u[0] / x, u[1] / x, u[2] / x};
}

double lgopt::delta(const size_t i, const size_t j)
{
    return i == j ? 1 : 0;
}

double lgopt::zeta(const size_t a, const size_t m, const size_t n)
{
    return delta(a, m) - delta(a, n);
}

bool lgopt::areParallel(const std::array<double, 3> &u, const std::array<double, 3> &v,
                        const double tol)
{
    std::array<double, 3> u_normed = u / norm(u);
    std::array<double, 3> v_normed = v / norm(v);

    double dp = dot(u_normed, v_normed);
    if (std::fabs(std::fabs(dp) - 1.0) < tol)
        return true;

    return false;
}

lgopt::HyperDual lgopt::operator+(const HyperDual &a, const HyperDual &b)
{
    const auto &[a0, a1, a2, a12] = a;
    const auto &[b0, b1, b2, b12] = b;

    return {a0 + b0, a1 + b1, a2 + b2, a12 + b12};
}

lgopt::HyperDual lgopt::operator-(const HyperDual &a, const HyperDual &b)
{
    const auto &[a0, a1, a2, a12] = a;
    const auto &[b0, b1, b2, b12] = b;

    return {a0 - b0, a1 - b1, a2 - b2, a12 - b12};
}

lgopt::HyperDual lgopt::operator*(const HyperDual &a, const HyperDual &b)
{
    const auto &[a0, a1, a2, a12] = a;
    const auto &[b0, b1, b2, b12] = b;

    return {
        a0 * b0, a0 * b1 + a1 * b0, a0 * b2 + a2 * b0,
        a0 * b12 + a1 * b2 + a2 * b1 + a12 * b0
    };
}

lgopt::HyperDual lgopt::operator*(const double f, const HyperDual &a)
{
    return {f * a.x0_, f * a.x1_, f * a.x2_, f * a.x12_};
}

lgopt::HyperDual lgopt::operator/(const HyperDual &a, const HyperDual &b)
{
    return a * inv(b);
}

lgopt::HyperDual lgopt::inv(const HyperDual &a)
{
    const auto &[a0, a1, a2, a12] = a;

    return {
        1.0 / a0, -a1 / (a0 * a0), -a2 / (a0 * a0),
        2 * a1 * a2 / (a0 * a0 * a0) - a12 / (a0 * a0)
    };
}

lgopt::HyperDual lgopt::sqrt(const HyperDual &a)
{
    const auto &[a0, a1, a2, a12] = a;

    return {
        std::sqrt(a0), a1 / (2 * std::sqrt(a0)), a2 / (2 * std::sqrt(a0)),
        a12 / (2 * std::sqrt(a0)) - a1 * a2 / (4 * a0 * std::sqrt(a0))
    };
}

lgopt::HyperDual lgopt::sin(const HyperDual &a)
{
    const auto &[a0, a1, a2, a12] = a;

    return {
        std::sin(a0), a1 * std::cos(a0), a2 * std::cos(a0),
        (a12 * std::cos(a0) - a1 * a2 * std::sin(a0))
    };
}

lgopt::HyperDual lgopt::acos(const HyperDual &a)
{
    // TODO: error handling
    const auto &[a0, a1, a2, a12] = a;

    return {
        std::acos(a0), -a1 / std::sqrt(1 - a0 * a0), -a2 / std::sqrt(1 - a0 * a0),
        -(a12 / std::sqrt(1 - a0 * a0) + a0 * a1 * a2 / std::pow(1 - a0 * a0, 1.5)) // Gambled it right lul
    };
    // return {
    //     std::acos(a0), -a1 / std::sqrt(1 - a0 * a0), -a2 / std::sqrt(1 - a0 * a0),
    //     -a0 * a1 * a2 / std::pow(1 - a0 * a0, 1.5)
    // };
}

lgopt::HyperDual lgopt::norm(const hd3_t &u)
{
    const auto &[x, y, z] = u;

    return sqrt(x * x + y * y + z * z);
}

lgopt::hd3_t lgopt::operator-(const hd3_t &u, const hd3_t &v)
{
    return {{u[0] - v[0], u[1] - v[1], u[2] - v[2]}};
}

lgopt::hd3_t lgopt::operator/(const hd3_t &u, const HyperDual &x)
{
    return {u[0] / x, u[1] / x, u[2] / x};
}

lgopt::HyperDual lgopt::bondLengthHD(const hd3_t &xyz_m, const hd3_t &xyz_n)
{
    return norm(xyz_m - xyz_n);
}

lgopt::HyperDual lgopt::bondAngleHD(const hd3_t &xyz_m, const hd3_t &xyz_o, const hd3_t &xyz_n)
{
    hd3_t u = xyz_m - xyz_o;
    hd3_t v = xyz_n - xyz_o;

    return acos(dot(u, v) / (norm(u) * norm(v)));
}

lgopt::HyperDual lgopt::dihedralAngleHD(const hd3_t &xyz_m, const hd3_t &xyz_o, const hd3_t &xyz_p,
                                        const hd3_t &xyz_n)
{
    // TODO: error handling
    hd3_t u = xyz_m - xyz_o;
    hd3_t v = xyz_n - xyz_p;
    hd3_t w = xyz_p - xyz_o;

    u = u / norm(u);
    v = v / norm(v);
    w = w / norm(w);

    HyperDual dot_uw = dot(u, w);
    HyperDual dot_vw = dot(v, w);

    HyperDual sin_u = sqrt(HyperDual(1) - dot_uw * dot_uw);
    HyperDual sin_v = sqrt(HyperDual(1) - dot_vw * dot_vw);

    HyperDual arg = dot(cross(u, w), cross(v, w)) / (sin_u * sin_v);
    return acos(arg);


    // return acos(arg);
}

lgopt::HyperDual lgopt::dot(const hd3_t &u, const hd3_t &v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

lgopt::hd3_t lgopt::cross(const hd3_t &u, const hd3_t &v)
{
    return {{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]}};
}

std::array<double, 6> lgopt::bondLengthGradient(const std::array<double, 3> &xyz_m,
                                                const std::array<double, 3> &xyz_n)
{
    std::array<double, 6> result;
    for (int i = 0; i < 6; i++)
    {
        std::array<HyperDual, 6> coords_hd{
            HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
            HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2]),
        };

        coords_hd[i] = {coords_hd[i].x0_, 1, 0, 0};

        HyperDual bond_length = bondLengthHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                             {coords_hd[3], coords_hd[4], coords_hd[5]});

        result[i] = bond_length.x1_;
    }

    return result;
}

std::array<double, 9> lgopt::bondAngleGradient(const std::array<double, 3> &xyz_m,
                                               const std::array<double, 3> &xyz_o,
                                               const std::array<double, 3> &xyz_n)
{
    std::array<double, 9> result;
    for (int i = 0; i < 9; i++)
    {
        std::array<HyperDual, 9> coords_hd{
            HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
            HyperDual(xyz_o[0]), HyperDual(xyz_o[1]), HyperDual(xyz_o[2]),
            HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2]),
        };

        coords_hd[i] = {coords_hd[i].x0_, 1, 0, 0};

        HyperDual bond_angle = bondAngleHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                           {coords_hd[3], coords_hd[4], coords_hd[5]},
                                           {coords_hd[6], coords_hd[7], coords_hd[8]});

        result[i] = bond_angle.x1_;
    }

    return result;
}

std::array<double, 12> lgopt::dihedralAngleGradient(const std::array<double, 3> &xyz_m,
                                                    const std::array<double, 3> &xyz_o,
                                                    const std::array<double, 3> &xyz_p,
                                                    const std::array<double, 3> &xyz_n)
{
    std::array<double, 12> result;
    for (int i = 0; i < 12; i++)
    {
        std::array<HyperDual, 12> coords_hd{
            HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
            HyperDual(xyz_o[0]), HyperDual(xyz_o[1]), HyperDual(xyz_o[2]),
            HyperDual(xyz_p[0]), HyperDual(xyz_p[1]), HyperDual(xyz_p[2]),
            HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2]),
        };

        coords_hd[i] = {coords_hd[i].x0_, 1, 0, 0};

        HyperDual dihedral_angle = dihedralAngleHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                                   {coords_hd[3], coords_hd[4], coords_hd[5]},
                                                   {coords_hd[6], coords_hd[7], coords_hd[8]},
                                                   {coords_hd[9], coords_hd[10], coords_hd[11]});

        result[i] = dihedral_angle.x1_;
    }

    return result;
}

lible::arr2d<double, 6, 6> lgopt::bondLengthHessian(const std::array<double, 3> &xyz_m,
                                                    const std::array<double, 3> &xyz_n)
{
    arr2d<double, 6, 6> result{};

    // Diagonals
    for (int ideriv = 0; ideriv < 6; ideriv++)
    {
        std::array<HyperDual, 6> coords_hd{
            HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
            HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2])
        };

        coords_hd[ideriv] = {coords_hd[ideriv].x0_, 1, 1, 0};

        HyperDual bond_length = bondLengthHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                             {coords_hd[3], coords_hd[4], coords_hd[5]});

        result[ideriv][ideriv] = bond_length.x12_;
    }

    // Off-diagonals
    for (int ideriv = 0; ideriv < 6; ideriv++)
        for (int jderiv = ideriv + 1; jderiv < 6; jderiv++)
        {
            std::array<HyperDual, 6> coords_hd{
                HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
                HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2])
            };

            coords_hd[ideriv] = {coords_hd[ideriv].x0_, 1, 0, 0};
            coords_hd[jderiv] = {coords_hd[jderiv].x0_, 0, 1, 0};

            HyperDual bond_length = bondLengthHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                                 {coords_hd[3], coords_hd[4], coords_hd[5]});

            result[ideriv][jderiv] = bond_length.x12_;
            result[jderiv][ideriv] = result[ideriv][jderiv];
        }

    return result;
}

lible::arr2d<double, 9, 9> lgopt::bondAngleHessian(const std::array<double, 3> &xyz_m,
                                                   const std::array<double, 3> &xyz_o,
                                                   const std::array<double, 3> &xyz_n)
{
    arr2d<double, 9, 9> result{};

    // Diagonals
    for (int ideriv = 0; ideriv < 9; ideriv++)
    {
        std::array<HyperDual, 9> coords_hd{
            HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
            HyperDual(xyz_o[0]), HyperDual(xyz_o[1]), HyperDual(xyz_o[2]),
            HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2])
        };

        coords_hd[ideriv] = {coords_hd[ideriv].x0_, 1, 1, 0};

        HyperDual bond_angle = bondAngleHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                           {coords_hd[3], coords_hd[4], coords_hd[5]},
                                           {coords_hd[6], coords_hd[7], coords_hd[8]});

        result[ideriv][ideriv] = bond_angle.x12_;
    }

    // Off-diagonals
    for (int ideriv = 0; ideriv < 9; ideriv++)
        for (int jderiv = ideriv + 1; jderiv < 9; jderiv++)
        {
            std::array<HyperDual, 9> coords_hd{
                HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
                HyperDual(xyz_o[0]), HyperDual(xyz_o[1]), HyperDual(xyz_o[2]),
                HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2])
            };

            coords_hd[ideriv] = {coords_hd[ideriv].x0_, 1, 0, 0};
            coords_hd[jderiv] = {coords_hd[jderiv].x0_, 0, 1, 0};

            HyperDual bond_angle = bondAngleHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                               {coords_hd[3], coords_hd[4], coords_hd[5]},
                                               {coords_hd[6], coords_hd[7], coords_hd[8]});

            result[ideriv][jderiv] = bond_angle.x12_;
            result[jderiv][ideriv] = result[ideriv][jderiv];
        }

    return result;
}

lible::arr2d<double, 12, 12> lgopt::dihedralAngleHessian(const std::array<double, 3> &xyz_m,
                                                         const std::array<double, 3> &xyz_o,
                                                         const std::array<double, 3> &xyz_p,
                                                         const std::array<double, 3> &xyz_n)
{
    arr2d<double, 12, 12> result{};

    // Diagonals
    for (int i = 0; i < 12; i++)
    {
        std::array<HyperDual, 12> coords_hd{
            HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
            HyperDual(xyz_o[0]), HyperDual(xyz_o[1]), HyperDual(xyz_o[2]),
            HyperDual(xyz_p[0]), HyperDual(xyz_p[1]), HyperDual(xyz_p[2]),
            HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2]),
        };

        coords_hd[i] = {coords_hd[i].x0_, 1, 1, 0};

        HyperDual dihedral_angle = dihedralAngleHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                                   {coords_hd[3], coords_hd[4], coords_hd[5]},
                                                   {coords_hd[6], coords_hd[7], coords_hd[8]},
                                                   {coords_hd[9], coords_hd[10], coords_hd[11]});

        result[i][i] = dihedral_angle.x12_;
    }

    // Off-diagonals
    for (int i = 0; i < 12; i++)
        for (int j = i + 1; j < 12; j++)
        {
            std::array<HyperDual, 12> coords_hd{
                HyperDual(xyz_m[0]), HyperDual(xyz_m[1]), HyperDual(xyz_m[2]),
                HyperDual(xyz_o[0]), HyperDual(xyz_o[1]), HyperDual(xyz_o[2]),
                HyperDual(xyz_p[0]), HyperDual(xyz_p[1]), HyperDual(xyz_p[2]),
                HyperDual(xyz_n[0]), HyperDual(xyz_n[1]), HyperDual(xyz_n[2]),
            };

            coords_hd[i] = {coords_hd[i].x0_, 1, 0, 0};
            coords_hd[j] = {coords_hd[j].x0_, 0, 1, 0};

            HyperDual dihedral_angle = dihedralAngleHD({coords_hd[0], coords_hd[1], coords_hd[2]},
                                                       {coords_hd[3], coords_hd[4], coords_hd[5]},
                                                       {coords_hd[6], coords_hd[7], coords_hd[8]},
                                                       {coords_hd[9], coords_hd[10], coords_hd[11]});

            result[i][j] = dihedral_angle.x12_;
            result[j][i] = result[i][j];
        }

    return result;
}
