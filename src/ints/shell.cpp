#include <lible/shell.hpp>
#include <lible/ints_util.hpp>
#include <lible/mcmurchie_davidson.hpp>
#include <lible/spherical_trafo.hpp>

namespace LI = lible::ints;
namespace LIMD = lible::ints::MD;

using std::array, std::size_t, std::vector;

vector<double> LI::calcShellNormalization(const int &angmom,
                                          const vector<double> &coeffs,
                                          const vector<double> &exps)
{
    int dim_cart = dimCartesians(angmom);
    int dim_sph = dimSphericals(angmom);

    vector<arma::dmat> h_coeffs;
    LIMD::calcHCoeffs(angmom, exps, h_coeffs);

    arma::dmat ints_contracted(dim_cart, dim_cart, arma::fill::zeros);

    size_t k = exps.size();

    size_t iab = 0;
    for (size_t ia = 0; ia < k; ia++)
        for (size_t ib = 0; ib < k; ib++)
        {
            double a = exps[ia];
            double b = exps[ib];
            double da = coeffs[ia];
            double db = coeffs[ib];

            double p = a + b;
            double da_x_db = da * db;

            double fac = da_x_db * std::pow(M_PI / p, 1.5);

            for (int mu = 0; mu < dim_cart; mu++)
                for (int nu = 0; nu < dim_cart; nu++)
                    ints_contracted(mu, nu) += fac * h_coeffs[iab](mu, nu);

            iab++;
        }

    arma::dmat sph_trafo1 = returnSphericalTrafo(angmom);
    arma::dmat sph_trafo2 = sph_trafo1.t();
    arma::dmat ints_sph = sph_trafo1 * ints_contracted * sph_trafo2;

    vector<double> normalization(dim_sph);
    for (int i = 0; i < dim_sph; i++)
        normalization[i] = 1.0 / std::sqrt(ints_sph(i, i));

    return normalization;
}