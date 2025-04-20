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
    void dipoleMomentKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                            const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                            const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                            const double *origin, double *ints_contracted);
}

void LI::dipoleMomentKernel(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                            const double *cexps_a, const double *cexps_b, const double *ccoeffs_a,
                            const double *ccoeffs_b, const double *xyz_a, const double *xyz_b,
                            const double *origin, double *ints_contracted)
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

                    int munu = mu * n_b_cart + nu;

                    ints_contracted[0 * n_ab_cart + munu] += valx;
                    ints_contracted[1 * n_ab_cart + munu] += valy;
                    ints_contracted[2 * n_ab_cart + munu] += valz;
                }
        }
}

template <>
void LIO::kernel<LIO::Option::dipole_moment, array<double, 3>>(const int la, const int lb,
                                                               const ShellPairData &sp_data,
                                                               const array<double, 3> &origin,
                                                               array<lible::vec2d, 3> &ints_out)
{
    int dim_a_cart = numCartesians(sp_data.la);
    int dim_b_cart = numCartesians(sp_data.lb);
    int dim_ab_cart = dim_a_cart * dim_b_cart;

    vector<CartExps> cart_exps_a = cart_exps[la];
    vector<CartExps> cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(sp_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(sp_data.lb).t();

    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        int cdepth_a = sp_data.cdepths[2 * ipair + 0];
        int cdepth_b = sp_data.cdepths[2 * ipair + 1];
        int cofs_a = sp_data.coffsets[2 * ipair + 0];
        int cofs_b = sp_data.coffsets[2 * ipair + 1];

        vector<double> ints_contracted(dim_ab_cart * 3);
        dipoleMomentKernel(la, lb, cdepth_a, cdepth_b, &sp_data.exps[cofs_a], &sp_data.exps[cofs_b],
                           &sp_data.coeffs[cofs_a], &sp_data.coeffs[cofs_b],
                           &sp_data.coords[6 * ipair + 0], &sp_data.coords[6 * ipair + 3],
                           origin.data(), ints_contracted.data());

        std::array<arma::dmat, 3> ints_contracted_arma;
        for (int i = 0; i < 3; i++)
            ints_contracted_arma[i] = arma::dmat(ints_contracted.data() + i * dim_ab_cart, dim_b_cart, dim_a_cart);

        std::array<arma::dmat, 3> ints_sph;
        ints_sph[0] = sph_trafo_bra * ints_contracted_arma[0].t() * sph_trafo_ket;
        ints_sph[1] = sph_trafo_bra * ints_contracted_arma[1].t() * sph_trafo_ket;
        ints_sph[2] = sph_trafo_bra * ints_contracted_arma[2].t() * sph_trafo_ket;

        transferInts1El(ipair, sp_data, ints_sph[0], ints_out[0]);
        transferInts1El(ipair, sp_data, ints_sph[1], ints_out[1]);
        transferInts1El(ipair, sp_data, ints_sph[2], ints_out[2]);
    }
}