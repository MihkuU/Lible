#include <lible/oneel_detail.hpp>
#include <lible/cart_exps.hpp>
#include <lible/mcmurchie_davidson.hpp>
#include <lible/spherical_trafo.hpp>

#include <math.h>

namespace LIO = lible::ints::one;

using std::array, std::size_t, std::vector;

template <>
void LIO::kernel<LIO::Option::overlap>(const int la, const int lb,
                                       const ShellPairData &shell_pair_data, vec2d &ints_out)
{
    vector<vector<vec3d>> ecoeffs;
    MD::calcECoeffs(la, lb, shell_pair_data, ecoeffs);

    int dim_a_cart = dimCartesians(shell_pair_data.la);
    int dim_b_cart = dimCartesians(shell_pair_data.lb);

    auto cart_exps_a = cart_exps[la];
    auto cart_exps_b = cart_exps[lb];

    arma::dmat sph_trafo_bra = returnSphericalTrafo(shell_pair_data.la);
    arma::dmat sph_trafo_ket = returnSphericalTrafo(shell_pair_data.lb).t();

    arma::dmat ints_contracted(dim_a_cart, dim_b_cart);
    arma::dmat ints_sph;
    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        ints_contracted.zeros();

        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
        const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
        const auto &[xyz_a, xyz_b] = shell_pair_data.coords[ipair];

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();

        const vector<vec3d> &Exyz = ecoeffs[ipair];

        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double p = exps_a[ia] + exps_b[ib];
                double dadb = ccoeffs_a[ia] * ccoeffs_b[ib];

                double fac = dadb * std::pow(M_PI / p, 1.5);

                for (const auto [i, j, k, mu] : cart_exps_a)
                    for (const auto [i_, j_, k_, nu] : cart_exps_b)
                        ints_contracted(mu, nu) += fac *
                                                   Exyz[iab](0, i, i_) *
                                                   Exyz[iab](1, j, j_) *
                                                   Exyz[iab](2, k, k_);
            }

        ints_sph = sph_trafo_bra * ints_contracted * sph_trafo_ket;

        transferIntegrals(ipair, shell_pair_data, ints_sph, ints_out);
    }
}