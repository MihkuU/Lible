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

lgopt::Dual lgopt::operator+(const Dual &a, const Dual &b)
{
    return {a.x0_ + b.x0_, a.x1_ + b.x1_};
}

lgopt::Dual lgopt::operator-(const Dual &a, const Dual &b)
{
    return {a.x0_ - b.x0_, a.x1_ - b.x1_};
}

lgopt::Dual lgopt::operator*(const Dual &a, const Dual &b)
{
    return {a.x0_ * b.x0_, a.x0_ * b.x1_ + a.x1_ * b.x0_};
}

lgopt::Dual lgopt::operator*(const double f, const Dual &a)
{
    return {a.x0_ * f, a.x1_ * f};
}

lgopt::Dual lgopt::inv(const Dual &a)
{
    return {1.0 / a.x0_, -a.x1_ / (a.x0_ * a.x0_)};
}

lgopt::Dual lgopt::operator/(const Dual &a, const Dual &b)
{
    return a * inv(b);
}

lgopt::Dual lgopt::sqrt(const Dual &a)
{
    return {std::sqrt(a.x0_), a.x1_ / (2 * std::sqrt(a.x0_))};
}

lgopt::Dual lgopt::sin(const Dual &a)
{
    return {std::sin(a.x0_), a.x1_ * std::cos(a.x0_)};
}

lgopt::Dual lgopt::acos(const Dual &a)
{
    return {std::acos(a.x0_), -a.x1_ / std::sqrt(1 - a.x0_ * a.x0_)};
}

lgopt::Dual lgopt::atan2(const Dual &y, const Dual &x)
{
    const auto &[y0, y1] = y;
    const auto &[x0, x1] = x;

    double denom = x0 * x0 + y0 * y0;

    return {std::atan2(y0, x0), (x0 * y1 - y0 * x1) / denom};
}

lgopt::Dual lgopt::norm(const dual3_t &u)
{
    const auto &[x, y, z] = u;

    return sqrt(x * x + y * y + z * z);
}

lgopt::Dual lgopt::dot(const dual3_t &u, const dual3_t &v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

lgopt::dual3_t lgopt::cross(const dual3_t &u, const dual3_t &v)
{
    return {{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]}};
}

lgopt::dual3_t lgopt::operator-(const dual3_t &u, const dual3_t &v)
{
    return {u[0] - v[0], u[1] - v[1], u[2] - v[2]};
}

lgopt::dual3_t lgopt::operator/(const dual3_t &u, const Dual &x)
{
    return {u[0] / x, u[1] / x, u[2] / x};
}

lgopt::Dual lgopt::dihedralAngleDual(const dual3_t &xyz_m, const dual3_t &xyz_o,
                                     const dual3_t &xyz_p, const dual3_t &xyz_n)
{
    // // TODO: error handling
    // dual3_t u = xyz_m - xyz_o;
    // dual3_t v = xyz_n - xyz_p;
    // dual3_t w = xyz_p - xyz_o;
    //
    // u = u / norm(u);
    // v = v / norm(v);
    // w = w / norm(w);
    //
    // Dual dot_uw = dot(u, w);
    // Dual dot_vw = dot(v, w);
    // Dual sin_u = sqrt(Dual(1) - dot_uw * dot_uw);
    // Dual sin_v = sqrt(Dual(1) - dot_vw * dot_vw);
    //
    // Dual arg = dot(cross(u, w), cross(v, w)) / (sin_u * sin_v);
    // return acos(arg);

    // New definition:
    dual3_t b1 = xyz_o - xyz_m;
    dual3_t b2 = xyz_p - xyz_o;
    dual3_t b3 = xyz_n - xyz_p;

    dual3_t n1 = cross(b1, b2);
    dual3_t n2 = cross(b2, b3);
    n1 = n1 / norm(n1);
    n2 = n2 / norm(n2);

    b2 = b2 / norm(b2);
    dual3_t m1 = cross(n1, b2);

    Dual x = dot(n1, n2);
    Dual y = dot(m1, n2);

    return atan2(y, x);
}

std::array<double, 12> lgopt::dihedralAngleGradientDual(const std::array<double, 3> &xyz_m,
                                                        const std::array<double, 3> &xyz_o,
                                                        const std::array<double, 3> &xyz_p,
                                                        const std::array<double, 3> &xyz_n)
{
    std::array<double, 12> result{};
    for (int i = 0; i < 12; i++)
    {
        std::array<Dual, 12> coords_dual{
            Dual(xyz_m[0]), Dual(xyz_m[1]), Dual(xyz_m[2]),
            Dual(xyz_o[0]), Dual(xyz_o[1]), Dual(xyz_o[2]),
            Dual(xyz_p[0]), Dual(xyz_p[1]), Dual(xyz_p[2]),
            Dual(xyz_n[0]), Dual(xyz_n[1]), Dual(xyz_n[2])
        };

        coords_dual[i] = {coords_dual[i].x0_, 1};

        Dual dihedral_angle = dihedralAngleDual(
            {coords_dual[0], coords_dual[1], coords_dual[2]},
            {coords_dual[3], coords_dual[4], coords_dual[5]},
            {coords_dual[6], coords_dual[7], coords_dual[8]},
            {coords_dual[9], coords_dual[10], coords_dual[11]});

        result[i] = dihedral_angle.x1_;
    }

    return result;
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
}

lgopt::HyperDual lgopt::atan2(const HyperDual &y, const HyperDual &x)
{
    // TODO: error handling
    const auto &[y0, y1, y2, y12] = y;
    const auto &[x0, x1, x2, x12] = x;

    double denom = x0 * x0 + y0 * y0;
    double term0 = std::atan2(y0, x0);
    double term1 = (x0 * y1 - y0 * x1) / denom;
    double term2 = (x0 * y2 - y0 * x2) / denom;
    double term12 = (x0 * y12 - y0 * x12 + x1 * y2 - y1 * x2) / denom
                    - 2 * (x0 * x1 + y0 * y1) * (x0 * y2 - y0 * x2) / (denom * denom);

    return {term0, term1, term2, term12};
}


lgopt::HyperDual lgopt::norm(const hd3_t &u)
{
    const auto &[x, y, z] = u;

    return sqrt(x * x + y * y + z * z);
}

lgopt::hd3_t lgopt::operator-(const hd3_t &u, const hd3_t &v)
{
    return {u[0] - v[0], u[1] - v[1], u[2] - v[2]};
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
    // // TODO: error handling
    // hd3_t u = xyz_m - xyz_o;
    // hd3_t v = xyz_n - xyz_p;
    // hd3_t w = xyz_p - xyz_o;
    //
    // u = u / norm(u);
    // v = v / norm(v);
    // w = w / norm(w);
    //
    // HyperDual dot_uw = dot(u, w);
    // HyperDual dot_vw = dot(v, w);
    //
    // HyperDual sin_u = sqrt(HyperDual(1) - dot_uw * dot_uw);
    // HyperDual sin_v = sqrt(HyperDual(1) - dot_vw * dot_vw);
    //
    // HyperDual arg = dot(cross(u, w), cross(v, w)) / (sin_u * sin_v);
    // return acos(arg);

    // New definition:
    hd3_t b1 = xyz_o - xyz_m;
    hd3_t b2 = xyz_p - xyz_o;
    hd3_t b3 = xyz_n - xyz_p;

    hd3_t n1 = cross(b1, b2);
    hd3_t n2 = cross(b2, b3);
    n1 = n1 / norm(n1);
    n2 = n2 / norm(n2);

    b2 = b2 / norm(b2);
    hd3_t m1 = cross(n1, b2);

    HyperDual x = dot(n1, n2);
    HyperDual y = dot(m1, n2);

    return atan2(y, x);
}

lgopt::HyperDual lgopt::dot(const hd3_t &u, const hd3_t &v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

lgopt::hd3_t lgopt::cross(const hd3_t &u, const hd3_t &v)
{
    return {{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]}};
}

std::array<double, 6> lgopt::bondLengthGradientHD(const std::array<double, 3> &xyz_m,
                                                  const std::array<double, 3> &xyz_n)
{
    std::array<double, 6> result{};
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

std::array<double, 9> lgopt::bondAngleGradientHD(const std::array<double, 3> &xyz_m,
                                                 const std::array<double, 3> &xyz_o,
                                                 const std::array<double, 3> &xyz_n)
{
    std::array<double, 9> result{};
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

std::array<double, 12> lgopt::dihedralAngleGradientHD(const std::array<double, 3> &xyz_m,
                                                      const std::array<double, 3> &xyz_o,
                                                      const std::array<double, 3> &xyz_p,
                                                      const std::array<double, 3> &xyz_n)
{
    std::array<double, 12> result{};
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

lible::arr2d<double, 6, 6> lgopt::bondLengthHessianHD(const std::array<double, 3> &xyz_m,
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

lible::arr2d<double, 9, 9> lgopt::bondAngleHessianHD(const std::array<double, 3> &xyz_m,
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

lible::arr2d<double, 12, 12> lgopt::dihedralAngleHessianHD(const std::array<double, 3> &xyz_m,
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
