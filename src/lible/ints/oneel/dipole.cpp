#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

namespace LIO = lible::ints::one;

using std::array, std::vector;

template <>
void LIO::kernel<LIO::Option::dipole_moment, array<double, 3>>(const int la, const int lb,
                                                               const ShellPairData &sp_data,
                                                               array<lible::vec2d, 3> &ints_out,
                                                               const array<double, 3> &origin)
{
    vector<vector<vec4d>> ecoeffs;
    ecoeffsShellPairs4D(la, lb, 1, sp_data, ecoeffs);

    int dim_a_cart = numCartesians(sp_data.la);
    int dim_b_cart = numCartesians(sp_data.lb);

    vector<CartExps> cart_exps_a = cart_exps[la];
    vector<CartExps> cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(sp_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(sp_data.lb).t();

    std::array<arma::dmat, 3> ints_contracted;
    for (int i = 0; i < 3; i++)
        ints_contracted[i].resize(dim_a_cart, dim_b_cart);

    std::array<arma::dmat, 3> ints_sph;
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        for (int i = 0; i < 3; i++)
            ints_contracted[i].zeros();

        array<double, 3> xyz_a{sp_data.coords[6 * ipair + 0],
                               sp_data.coords[6 * ipair + 1],
                               sp_data.coords[6 * ipair + 2]};

        array<double, 3> xyz_b{sp_data.coords[6 * ipair + 3],
                               sp_data.coords[6 * ipair + 4],
                               sp_data.coords[6 * ipair + 5]};

        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        const vector<vec4d> &Exyz = ecoeffs[ipair];

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double p = a + b;

                array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                       (a * xyz_a[1] + b * xyz_b[1]) / p,
                                       (a * xyz_a[2] + b * xyz_b[2]) / p};

                array<double, 3> xyz_po{xyz_p[0] - origin[0],
                                        xyz_p[1] - origin[1],
                                        xyz_p[2] - origin[2]};

                double da = sp_data.coeffs[pos_a + ia];
                double db = sp_data.coeffs[pos_b + ib];

                double dadb = da * db;
                double fac = dadb * std::pow(M_PI / p, 1.5);

                for (const auto &[i, j, k, mu] : cart_exps_a)
                    for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                    {
                        double valx = fac *
                                      (Exyz[iab](0, i, i_, 1) + xyz_po[0] * Exyz[iab](0, i, i_, 0)) *
                                      Exyz[iab](1, j, j_, 0) * Exyz[iab](2, k, k_, 0);

                        double valy = fac *
                                      (Exyz[iab](0, i, i_, 0)) *
                                      (Exyz[iab](1, j, j_, 1) + xyz_po[1] * Exyz[iab](1, j, j_, 0)) *
                                      Exyz[iab](2, k, k_, 0);

                        double valz = fac *
                                      (Exyz[iab](0, i, i_, 0)) * Exyz[iab](1, j, j_, 1) *
                                      (Exyz[iab](2, k, k_, 1) + xyz_po[2] * Exyz[iab](2, k, k_, 0));

                        ints_contracted[0](mu, nu) += valx;
                        ints_contracted[1](mu, nu) += valy;
                        ints_contracted[2](mu, nu) += valz;
                    }
            }

        ints_sph[0] = sph_trafo_bra * ints_contracted[0] * sph_trafo_ket;
        ints_sph[1] = sph_trafo_bra * ints_contracted[1] * sph_trafo_ket;
        ints_sph[2] = sph_trafo_bra * ints_contracted[2] * sph_trafo_ket;

        transferInts1El(ipair, sp_data, ints_sph[0], ints_out[0]);
        transferInts1El(ipair, sp_data, ints_sph[1], ints_out[1]);
        transferInts1El(ipair, sp_data, ints_sph[2], ints_out[2]);
    }
}