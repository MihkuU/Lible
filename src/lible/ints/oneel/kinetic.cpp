#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <math.h>

namespace LI = lible::ints;
namespace LIO = lible::ints::one;

using std::array, std::vector;

namespace lible::ints
{
    // Forward declaration
    vec2d kineticEnergyKernel(const int ipair, const ShellPairData &sp_data);     

    // Forward declaration
    array<vec2d, 6> kineticEnergyD1Kernel(const int ipair, const ShellPairData &sp_data);

    // A helper function for 'kineticEnergyD1Kernel'.
    auto kinEDeriv1Helper = [](const double b, const double b2, const double fac,
                               const vec3d &ecoeffs_x, const vec3d &ecoeffs_y,
                               const vec3d &ecoeffs_z, const array<int, 3> &ijk,
                               const array<int, 3> &i_j_k_)
    {
        auto [i, j, k] = ijk;
        auto [i_, j_, k_] = i_j_k_;

        double Tx, Ty, Tz;
        if (i_ < 2)
            Tx = -2 * b2 * ecoeffs_x(i, i_ + 2, 0) +
                 b * (2 * i_ + 1) * ecoeffs_x(i, i_, 0);
        else
            Tx = -2 * b2 * ecoeffs_x(i, i_ + 2, 0) +
                 b * (2 * i_ + 1) * ecoeffs_x(i, i_, 0) -
                 0.5 * i_ * (i_ - 1) * ecoeffs_x(i, i_ - 2, 0);

        if (j_ < 2)
            Ty = -2 * b2 * ecoeffs_y(j, j_ + 2, 0) +
                 b * (2 * j_ + 1) * ecoeffs_y(j, j_, 0);
        else
            Ty = -2 * b2 * ecoeffs_y(j, j_ + 2, 0) +
                 b * (2 * j_ + 1) * ecoeffs_y(j, j_, 0) -
                 0.5 * j_ * (j_ - 1) * ecoeffs_y(j, j_ - 2, 0);

        if (k_ < 2)
            Tz = -2 * b2 * ecoeffs_z(k, k_ + 2, 0) +
                 b * (2 * k_ + 1) * ecoeffs_z(k, k_, 0);
        else
            Tz = -2 * b2 * ecoeffs_z(k, k_ + 2, 0) +
                 b * (2 * k_ + 1) * ecoeffs_z(k, k_, 0) -
                 0.5 * k_ * (k_ - 1) * ecoeffs_z(k, k_ - 2, 0);

        double integral = 0;
        integral += fac * Tx * ecoeffs_y(j, j_, 0) * ecoeffs_z(k, k_, 0);
        integral += fac * ecoeffs_x(i, i_, 0) * Ty * ecoeffs_z(k, k_, 0);
        integral += fac * ecoeffs_x(i, i_, 0) * ecoeffs_y(j, j_, 0) * Tz;

        return integral;
    };
}

lible::vec2d LI::kineticEnergyKernel(const int ipair, const ShellPairData &sp_data)
{
    // Formula taken from https://gqcg-res.github.io/knowdes/the-mcmurchie-davidson-integral-scheme.html.

    int la = sp_data.la;
    int lb = sp_data.lb;
    int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    int cofs_a = sp_data.coffsets[2 * ipair + 0];
    int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int n_cart_a = numCartesians(la);
    int n_cart_b = numCartesians(lb);

    vec2d ints_cart(Fill(0), n_cart_a, n_cart_b);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double b2 = b * b;
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double dadb = da * db;
            double fac = dadb * std::pow(M_PI / p, 1.5);

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb + 2, xyz_a, xyz_b);

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                {
                    double Tx, Ty, Tz;
                    if (i_ < 2)
                        Tx = -2 * b2 * Ex(i, i_ + 2, 0) +
                             b * (2 * i_ + 1) * Ex(i, i_, 0);
                    else
                        Tx = -2 * b2 * Ex(i, i_ + 2, 0) +
                             b * (2 * i_ + 1) * Ex(i, i_, 0) -
                             0.5 * i_ * (i_ - 1) * Ex(i, i_ - 2, 0);

                    if (j_ < 2)
                        Ty = -2 * b2 * Ey(j, j_ + 2, 0) +
                             b * (2 * j_ + 1) * Ey(j, j_, 0);
                    else
                        Ty = -2 * b2 * Ey(j, j_ + 2, 0) +
                             b * (2 * j_ + 1) * Ey(j, j_, 0) -
                             0.5 * j_ * (j_ - 1) * Ey(j, j_ - 2, 0);

                    if (k_ < 2)
                        Tz = -2 * b2 * Ez(k, k_ + 2, 0) +
                             b * (2 * k_ + 1) * Ez(k, k_, 0);
                    else
                        Tz = -2 * b2 * Ez(k, k_ + 2, 0) +
                             b * (2 * k_ + 1) * Ez(k, k_, 0) -
                             0.5 * k_ * (k_ - 1) * Ez(k, k_ - 2, 0);

                    ints_cart(mu, nu) += fac * Tx * Ey(j, j_, 0) * Ez(k, k_, 0);
                    ints_cart(mu, nu) += fac * Ex(i, i_, 0) * Ty * Ez(k, k_, 0);
                    ints_cart(mu, nu) += fac * Ex(i, i_, 0) * Ey(j, j_, 0) * Tz;
                }
        }

    vec2d ints_sph = trafo2Spherical(la, lb, ints_cart);

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (size_t mu = 0; mu < ints_sph.dim<0>(); mu++)
        for (size_t nu = 0; nu < ints_sph.dim<1>(); nu++)
        {
            double norm_a = sp_data.norms[ofs_norm_a + mu];
            double norm_b = sp_data.norms[ofs_norm_b + nu];
            ints_sph(mu, nu) *= norm_a * norm_b;
        }

    return ints_sph;
}

array<lible::vec2d, 6> LI::kineticEnergyD1Kernel(const int ipair, const ShellPairData &sp_data)
{
    int la = sp_data.la;
    int lb = sp_data.lb;
    int cdepth_a = sp_data.cdepths[2 * ipair + 0];
    int cdepth_b = sp_data.cdepths[2 * ipair + 1];
    int cofs_a = sp_data.coffsets[2 * ipair + 0];
    int cofs_b = sp_data.coffsets[2 * ipair + 1];

    const double *cexps_a = &sp_data.exps[cofs_a];
    const double *cexps_b = &sp_data.exps[cofs_b];
    const double *ccoeffs_a = &sp_data.coeffs[cofs_a];
    const double *ccoeffs_b = &sp_data.coeffs[cofs_b];
    const double *xyz_a = &sp_data.coords[6 * ipair + 0];
    const double *xyz_b = &sp_data.coords[6 * ipair + 3];

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    int n_cart_a = numCartesians(la);
    int n_cart_b = numCartesians(lb);

    array<vec2d, 6> ints_cart;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_cart[ideriv] = vec2d(Fill(0), n_cart_a, n_cart_b);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double b2 = b * b;
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double dadb = da * db;
            double fac = dadb * std::pow(M_PI / p, 1.5);

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb + 2, xyz_a, xyz_b);

            auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb + 2, xyz_a, xyz_b,
                                                           {Ex, Ey, Ez});

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                {                    
                    double kin_x = kinEDeriv1Helper(b, b2, fac, E1x, Ey, Ez, {i, j, k}, {i_, j_, k_});
                    double kin_y = kinEDeriv1Helper(b, b2, fac, Ex, E1y, Ez, {i, j, k}, {i_, j_, k_});
                    double kin_z = kinEDeriv1Helper(b, b2, fac, Ex, Ey, E1z, {i, j, k}, {i_, j_, k_});

                    // d/dA
                    ints_cart[0](mu, nu) += kin_x;
                    ints_cart[1](mu, nu) += kin_y;
                    ints_cart[2](mu, nu) += kin_z;

                    // d/dB
                    ints_cart[3](mu, nu) -= kin_x;
                    ints_cart[4](mu, nu) -= kin_y;
                    ints_cart[5](mu, nu) -= kin_z;
                }
        }

    array<vec2d, 6> ints_sph;
    for (int ideriv = 0; ideriv < 6; ideriv++)
        ints_sph[ideriv] = trafo2Spherical(la, lb, ints_cart[ideriv]);

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (int ideriv = 0; ideriv < 6; ideriv++)
        for (size_t mu = 0; mu < ints_sph[ideriv].dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_sph[ideriv].dim<1>(); nu++)
            {
                double norm_a = sp_data.norms[ofs_norm_a + mu];
                double norm_b = sp_data.norms[ofs_norm_b + nu];
                ints_sph[ideriv](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

template <>
void LIO::kernel<LIO::Option::kinetic_energy>(const int la, const int lb,
                                              const ShellPairData &sp_data,
                                              vec2d &ints_out)
{
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        vec2d ints_ipair = kineticEnergyKernel(ipair, sp_data);

        int ofs_a = sp_data.offsets_sph[2 * ipair + 0];
        int ofs_b = sp_data.offsets_sph[2 * ipair + 1];
        for (size_t mu = 0; mu < ints_ipair.dim<0>(); mu++)
            for (size_t nu = 0; nu < ints_ipair.dim<1>(); nu++)
            {
                ints_out(ofs_a + mu, ofs_b + nu) = ints_ipair(mu, nu);
                ints_out(ofs_b + nu, ofs_a + mu) = ints_ipair(mu, nu);
            }
    }
}