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


size_t lgopt::delta(const size_t i, const size_t j)
{
    return i == j ? 1 : 0;
}

int lgopt::zeta(const size_t a, const size_t m, const size_t n)
{
    return static_cast<int>(delta(a, m)) - static_cast<int>(delta(a, n));
}


bool lgopt::areParallel(const std::array<double, 3> &u, const std::array<double, 3> &v,
                        const double tol)
{
    double ratio = u[0] / v[0];

    for (size_t i = 1; i < 3; ++i)
        if (std::fabs(u[i] - ratio * v[i]) > tol)
            return false;

    return true;
}
