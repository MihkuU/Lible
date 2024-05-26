#include <lible/oneel_detail.hpp>
#include <lible/cart_exps.hpp>
#include <lible/spherical_trafo.hpp>

#include <math.h>

namespace LIO = lible::ints::one;

using std::vector;

template <>
void LIO::kernel<LIO::Option::kinetic_energy>(const int la, const int lb,
                                              const ShellPairData &shell_pair_data,
                                              vec2d &ints_out)
{
    // Formula taken from https://gqcg-res.github.io/knowdes/the-mcmurchie-davidson-integral-scheme.html.

    vector<vector<vec3d>> ecoeffs;
    MD::calcECoeffs(la, lb + 2, shell_pair_data, ecoeffs);
    
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

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();

        const vector<vec3d> &Exyz = ecoeffs[ipair];

        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {
                double p = exps_a[ia] + exps_b[ib];
                double dadb = ccoeffs_a[ia] * ccoeffs_b[ib];

                double b = exps_b[ib];
                double b2 = std::pow(b, 2);

                double fac = dadb * std::pow(M_PI / p, 1.5);

                for (const auto [i, j, k, mu] : cart_exps_a)
                    for (const auto [i_, j_, k_, nu] : cart_exps_b)
                    {
                        double Tx, Ty, Tz;
                        if (i_ < 2)
                            Tx = -2 * b2 * Exyz[iab](0, i, i_ + 2) +
                                 b * (2 * i_ + 1) * Exyz[iab](0, i, i_);
                        else
                            Tx = -2 * b2 * Exyz[iab](0, i, i_ + 2) +
                                 b * (2 * i_ + 1) * Exyz[iab](0, i, i_) -
                                 0.5 * i_ * (i_ - 1) * Exyz[iab](0, i, i_ - 2);

                        if (j_ < 2)
                            Ty = -2 * b2 * Exyz[iab](1, j, j_ + 2) +
                                 b * (2 * j_ + 1) * Exyz[iab](1, j, j_);
                        else
                            Ty = -2 * b2 * Exyz[iab](1, j, j_ + 2) +
                                 b * (2 * j_ + 1) * Exyz[iab](1, j, j_) -
                                 0.5 * j_ * (j_ - 1) * Exyz[iab](1, j, j_ - 2);

                        if (k_ < 2)
                            Tz = -2 * b2 * Exyz[iab](2, k, k_ + 2) +
                                 b * (2 * k_ + 1) * Exyz[iab](2, k, k_);
                        else
                            Tz = -2 * b2 * Exyz[iab](2, k, k_ + 2) +
                                 b * (2 * k_ + 1) * Exyz[iab](2, k, k_) -
                                 0.5 * k_ * (k_ - 1) * Exyz[iab](2, k, k_ - 2);

                        ints_contracted(mu, nu) += fac * Tx * Exyz[iab](1, j, j_) * Exyz[iab](2, k, k_);
                        ints_contracted(mu, nu) += fac * Exyz[iab](0, i, i_) * Ty * Exyz[iab](2, k, k_);
                        ints_contracted(mu, nu) += fac * Exyz[iab](0, i, i_) * Exyz[iab](1, j, j_) * Tz;
                    }
            }

        ints_sph = sph_trafo_bra * ints_contracted * sph_trafo_ket;

        transferIntegrals(ipair, shell_pair_data, ints_sph, ints_out);
    }
}