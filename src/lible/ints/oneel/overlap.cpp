#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <math.h>

namespace LIO = lible::ints::one;

using std::array, std::size_t, std::vector;

template <>
void LIO::kernel<LIO::Option::overlap>(const int la, const int lb,
                                       const ShellPairData &sp_data, vec2d &ints_out)
{
    vector<vector<vec3d>> ecoeffs = ecoeffsSPData_Eij0(la, lb, sp_data);

    int dim_a_cart = numCartesians(sp_data.la);
    int dim_b_cart = numCartesians(sp_data.lb);

    vector<CartExps> cart_exps_a = cart_exps[la];
    vector<CartExps> cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(sp_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(sp_data.lb).t();

    arma::dmat ints_contracted(dim_a_cart, dim_b_cart);
    arma::dmat ints_sph;
    for (int ipair = 0; ipair < sp_data.n_pairs; ipair++)
    {
        ints_contracted.zeros();

        int dim_a = sp_data.cdepths[2 * ipair + 0];
        int dim_b = sp_data.cdepths[2 * ipair + 1];
        int pos_a = sp_data.coffsets[2 * ipair + 0];
        int pos_b = sp_data.coffsets[2 * ipair + 1];

        const vector<vec3d> &ecoeffs_ipair = ecoeffs[ipair];

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                double a = sp_data.exps[pos_a + ia];
                double b = sp_data.exps[pos_b + ib];
                double da = sp_data.coeffs[pos_a + ia];
                double db = sp_data.coeffs[pos_b + ib];

                double p = a + b;
                double dadb = da * db;
                double fac = dadb * std::pow(M_PI / p, 1.5);

                for (const auto &[i, j, k, mu] : cart_exps_a)
                    for (const auto &[i_, j_, k_, nu] : cart_exps_b)
                        ints_contracted(mu, nu) += fac *
                                                   ecoeffs_ipair[iab](0, i, i_) *
                                                   ecoeffs_ipair[iab](1, j, j_) *
                                                   ecoeffs_ipair[iab](2, k, k_);
            }

        ints_sph = sph_trafo_bra * ints_contracted * sph_trafo_ket;

        transferInts1El(ipair, sp_data, ints_sph, ints_out);
    }
}