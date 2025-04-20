#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <math.h>

namespace LI = lible::ints;
namespace LIO = lible::ints::one;

using std::array, std::size_t, std::vector;

namespace lible::ints
{
    // Forward declarations
    void overlapKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                       const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                       const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                       double *ints_contracted);
    
    void overlapDeriv1Kernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                             const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                             const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                             double *intderivs_contracted);
}

void LI::overlapKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                       const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                       const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                       double *ints_contracted)
{    
    int n_a_cart = numCartesians(la);
    int n_b_cart = numCartesians(lb);
    int n_ab_cart = n_a_cart * n_b_cart;

    std::fill(ints_contracted, ints_contracted + n_ab_cart, 0);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double mu = a * b / p;
            double dadb = da * db;
            double fac = dadb * std::pow(M_PI / p, 1.5);

            auto [ecoeffs_x, ecoeffs_y, ecoeffs_z] = ecoeffsPrimitivePair(a, b, la, lb, xyz_a, xyz_b);

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    ints_contracted[mu * n_b_cart + nu] += fac *
                                                           ecoeffs_x(i, i_, 0) *
                                                           ecoeffs_y(j, j_, 0) *
                                                           ecoeffs_z(k, k_, 0);
        }
}

void LI::overlapDeriv1Kernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                             const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                             const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                             double *intderivs_contracted)
{
    int n_a_cart = numCartesians(la);
    int n_b_cart = numCartesians(lb);
    int n_ab_cart = n_a_cart * n_b_cart;

    std::fill(intderivs_contracted, intderivs_contracted + 6 * n_ab_cart, 0);

    const auto &cart_exps_a = cart_exps[la];
    const auto &cart_exps_b = cart_exps[lb];

    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = cexps_a[ia];
            double b = cexps_b[ib];
            double da = ccoeffs_a[ia];
            double db = ccoeffs_b[ib];

            double p = a + b;
            double mu = a * b / p;
            double dadb = da * db;
            double fac = dadb * std::pow(M_PI / p, 1.5);

            auto [Ex, Ey, Ez] = ecoeffsPrimitivePair(a, b, la, lb + 2, xyz_a, xyz_b);

            auto [E1x, E1y, E1z] = ecoeffsPrimitivePair_n1(a, b, la, lb + 2, xyz_a, xyz_b, {Ex, Ey, Ez});

            for (const auto &[i, j, k, mu] : cart_exps_a)
                for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                {
                    int munu = mu * n_b_cart + nu;

                    // d/dA
                    intderivs_contracted[0 * n_ab_cart + munu] += fac * E1x(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
                    intderivs_contracted[1 * n_ab_cart + munu] += fac * Ex(i, i_, 0) * E1y(j, j_, 0) * Ez(k, k_, 0);
                    intderivs_contracted[2 * n_ab_cart + munu] += fac * Ex(i, i_, 0) * Ey(j, j_, 0) * E1z(k, k_, 0);

                    // d/dB
                    intderivs_contracted[3 * n_ab_cart + munu] -= fac * E1x(i, i_, 0) * Ey(j, j_, 0) * Ez(k, k_, 0);
                    intderivs_contracted[4 * n_ab_cart + munu] -= fac * Ex(i, i_, 0) * E1y(j, j_, 0) * Ez(k, k_, 0);
                    intderivs_contracted[5 * n_ab_cart + munu] -= fac * Ex(i, i_, 0) * Ey(j, j_, 0) * E1z(k, k_, 0);
                }
        }
}

template <>
void LIO::kernel<LIO::Option::overlap>(const int la, const int lb,
                                       const ShellPairData &sp_data, vec2d &ints_out)
{    
    int n_a_cart = numCartesians(sp_data.la);
    int n_b_cart = numCartesians(sp_data.lb);

    vector<CartExps> cart_exps_a = cart_exps[la];
    vector<CartExps> cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(sp_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(sp_data.lb).t();

    arma::dmat ints_contracted(n_b_cart, n_a_cart);
    arma::dmat ints_sph;
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        ints_contracted.zeros();

        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int cofs_a = sp_data.coffsets[2 * ipair + 0];
        int cofs_b = sp_data.coffsets[2 * ipair + 1];

        overlapKernel(la, lb, cdepth_a, cdepth_b, &sp_data.exps[cofs_a], &sp_data.exps[cofs_b],
                      &sp_data.coeffs[cofs_a], &sp_data.coeffs[cofs_b],
                      &sp_data.coords[6 * ipair + 0], &sp_data.coords[6 * ipair + 3],
                      ints_contracted.memptr());

        ints_sph = sph_trafo_bra * ints_contracted.t() * sph_trafo_ket;

        transferInts1El(ipair, sp_data, ints_sph, ints_out);
    }
}
