#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

namespace LI = lible::ints;
namespace LIO = lible::ints::one;

using std::array, std::vector;

namespace lible::ints
{
    // Forward declaration
    array<vec2d, 3> dipoleMomentKernel(const int ipair, const array<double, 3> &origin,
                                       const ShellPairData &sp_data);
}

array<lible::vec2d, 3> LI::dipoleMomentKernel(const int ipair, const array<double, 3> &origin,
                                              const ShellPairData &sp_data)
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

    array<vec2d, 3> ints_cart;
    for (int icart = 0; icart < 3; icart++)
        ints_cart[icart] = vec2d(n_cart_a, n_cart_b, 0);

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double dadb = da * db;
            double fac = dadb * std::pow(M_PI / p, 1.5);

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                   (a * xyz_a[1] + b * xyz_b[1]) / p,
                                   (a * xyz_a[2] + b * xyz_b[2]) / p};

            array<double, 3> xyz_po{xyz_p[0] - origin[0],
                                    xyz_p[1] - origin[1],
                                    xyz_p[2] - origin[2]};

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                {
                    double valx = fac *
                                  (Ex(i, i_, 1) + xyz_po[0] * Ex(i, i_, 0)) *
                                  Ey(j, j_, 0) * Ez(k, k_, 0);

                    double valy = fac *
                                  (Ex(i, i_, 0)) *
                                  (Ey(j, j_, 1) + xyz_po[1] * Ey(j, j_, 0)) *
                                  Ez(k, k_, 0);

                    double valz = fac *
                                  (Ex(i, i_, 0)) * Ey(j, j_, 0) *
                                  (Ez(k, k_, 1) + xyz_po[2] * Ez(k, k_, 0));                    

                    ints_cart[0](mu, nu) += valx;
                    ints_cart[1](mu, nu) += valy;
                    ints_cart[2](mu, nu) += valz;
                }
        }

    array<vec2d, 3> ints_sph;
    for (int icart = 0; icart < 3; icart++)
        ints_sph[icart] = trafo2Spherical(la, lb, ints_cart[icart]);

    int ofs_norm_a = sp_data.offsets_norms[2 * ipair + 0];
    int ofs_norm_b = sp_data.offsets_norms[2 * ipair + 1];
    for (int icart = 0; icart < 3; icart++)
        for (size_t mu = 0; mu < ints_sph[icart].getDim(0); mu++)
            for (size_t nu = 0; nu < ints_sph[icart].getDim(1); nu++)
            {
                double norm_a = sp_data.norms[ofs_norm_a + mu];
                double norm_b = sp_data.norms[ofs_norm_b + nu];
                ints_sph[icart](mu, nu) *= norm_a * norm_b;
            }

    return ints_sph;
}

template <>
void LIO::kernel<LIO::Option::dipole_moment, array<double, 3>>(const int la, const int lb,
                                                               const ShellPairData &sp_data,
                                                               const array<double, 3> &origin,
                                                               array<lible::vec2d, 3> &ints_out)
{
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        array<vec2d, 3> ints_ipair = dipoleMomentKernel(ipair, origin, sp_data);

        int ofs_a = sp_data.offsets_sph[2 * ipair + 0];
        int ofs_b = sp_data.offsets_sph[2 * ipair + 1];
        for (int icart = 0; icart < 3; icart++)
            for (size_t mu = 0; mu < ints_ipair[icart].getDim(0); mu++)
                for (size_t nu = 0; nu < ints_ipair[icart].getDim(1); nu++)
                {
                    ints_out[icart](ofs_a + mu, ofs_b + nu) = ints_ipair[icart](mu, nu);
                    ints_out[icart](ofs_b + nu, ofs_a + mu) = ints_ipair[icart](mu, nu);
                }
    }
}