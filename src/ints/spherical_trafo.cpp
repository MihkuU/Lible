#include <lible/spherical_trafo.hpp>
#include <lible/ints_defs.hpp>
#include <lible/ints_util.hpp>

#include <cassert>
#include <stdexcept>

namespace LI = lible::ints;

arma::dmat LI::returnSphericalTrafo(const int angmom)
{
    int dim_cart = dimCartesians(angmom);
    int dim_sph = dimSphericals(angmom);

    arma::dmat trafo(dim_sph, dim_cart, arma::fill::zeros);

    switch (angmom)
    {
    case (0):
    {
        trafo(0, 0) = 1.0000000000000000;

        return trafo;
    }
    case (1):
    {
        trafo(0, 2) = 1.0000000000000000;
        trafo(1, 0) = 1.0000000000000000;
        trafo(2, 1) = 1.0000000000000000;

        return trafo;
    }
    case (2):
    {
        trafo(0, 5) = 1.00000000000000;
        trafo(0, 0) = -0.50000000000000;
        trafo(0, 3) = -0.50000000000000;
        trafo(1, 2) = 1.73205080756888;
        trafo(2, 4) = 1.73205080756888;
        trafo(3, 0) = 0.86602540378444;
        trafo(3, 3) = -0.86602540378444;
        trafo(4, 1) = 1.73205080756888;

        return trafo;
    }
    case (3):
    {
        trafo(0, 2) = -1.50000000000000;
        trafo(0, 7) = -1.50000000000000;
        trafo(0, 9) = 1.00000000000000;
        trafo(1, 0) = -0.61237243569579;
        trafo(1, 3) = -0.61237243569579;
        trafo(1, 5) = 2.44948974278318;
        trafo(2, 1) = -0.61237243569579;
        trafo(2, 6) = -0.61237243569579;
        trafo(2, 8) = 2.44948974278318;
        trafo(3, 2) = 1.93649167310371;
        trafo(3, 7) = -1.93649167310371;
        trafo(4, 4) = 3.87298334620742;
        trafo(5, 0) = 0.79056941504209;
        trafo(5, 3) = -2.37170824512628;
        trafo(6, 1) = 2.37170824512628;
        trafo(6, 6) = -0.79056941504209;

        return trafo;
    }
    case (4):
    {
        trafo(0, 0) = 0.37500000000000;
        trafo(0, 3) = 0.75000000000000;
        trafo(0, 5) = -3.00000000000000;
        trafo(0, 10) = 0.37500000000000;
        trafo(0, 12) = -3.00000000000000;
        trafo(0, 14) = 1.00000000000000;
        trafo(1, 2) = -2.37170824512628;
        trafo(1, 7) = -2.37170824512628;
        trafo(1, 9) = 3.16227766016838;
        trafo(2, 4) = -2.37170824512628;
        trafo(2, 11) = -2.37170824512628;
        trafo(2, 13) = 3.16227766016838;
        trafo(3, 0) = -0.55901699437495;
        trafo(3, 5) = 3.35410196624968;
        trafo(3, 10) = 0.55901699437495;
        trafo(3, 12) = -3.35410196624968;
        trafo(4, 1) = -1.11803398874989;
        trafo(4, 6) = -1.11803398874989;
        trafo(4, 8) = 6.70820393249937;
        trafo(5, 2) = 2.09165006633519;
        trafo(5, 7) = -6.27495019900557;
        trafo(6, 4) = 6.27495019900557;
        trafo(6, 11) = -2.09165006633519;
        trafo(7, 0) = 0.73950997288745;
        trafo(7, 3) = -4.43705983732471;
        trafo(7, 10) = 0.73950997288745;
        trafo(8, 1) = 2.95803989154981;
        trafo(8, 6) = -2.95803989154981;

        return trafo;
    }
    case (5):
    {
        trafo(0, 2) = 1.87500000000000;
        trafo(0, 7) = 3.75000000000000;
        trafo(0, 9) = -5.00000000000000;
        trafo(0, 16) = 1.87500000000000;
        trafo(0, 18) = -5.00000000000000;
        trafo(0, 20) = 1.00000000000000;
        trafo(1, 0) = 0.48412291827593;
        trafo(1, 3) = 0.96824583655185;
        trafo(1, 5) = -5.80947501931113;
        trafo(1, 10) = 0.48412291827593;
        trafo(1, 12) = -5.80947501931113;
        trafo(1, 14) = 3.87298334620742;
        trafo(2, 1) = 0.48412291827593;
        trafo(2, 6) = 0.96824583655185;
        trafo(2, 8) = -5.80947501931113;
        trafo(2, 15) = 0.48412291827593;
        trafo(2, 17) = -5.80947501931113;
        trafo(2, 19) = 3.87298334620742;
        trafo(3, 2) = -2.56173769148990;
        trafo(3, 9) = 5.12347538297980;
        trafo(3, 16) = 2.56173769148990;
        trafo(3, 18) = -5.12347538297980;
        trafo(4, 4) = -5.12347538297980;
        trafo(4, 11) = -5.12347538297980;
        trafo(4, 13) = 10.24695076595960;
        trafo(5, 0) = -0.52291251658380;
        trafo(5, 3) = 1.04582503316759;
        trafo(5, 5) = 4.18330013267038;
        trafo(5, 10) = 1.56873754975139;
        trafo(5, 12) = -12.54990039801113;
        trafo(6, 1) = -1.56873754975139;
        trafo(6, 6) = -1.04582503316759;
        trafo(6, 8) = 12.54990039801113;
        trafo(6, 15) = 0.52291251658380;
        trafo(6, 17) = -4.18330013267038;
        trafo(7, 2) = 2.21852991866236;
        trafo(7, 7) = -13.31117951197414;
        trafo(7, 16) = 2.21852991866236;
        trafo(8, 4) = 8.87411967464942;
        trafo(8, 11) = -8.87411967464942;
        trafo(9, 0) = 0.70156076002011;
        trafo(9, 3) = -7.01560760020114;
        trafo(9, 10) = 3.50780380010057;
        trafo(10, 1) = 3.50780380010057;
        trafo(10, 6) = -7.01560760020114;
        trafo(10, 15) = 0.70156076002011;

        return trafo;
    }
    case (6):
    {
        trafo(0, 0) = -0.31250000000000;
        trafo(0, 3) = -0.93750000000000;
        trafo(0, 5) = 5.62500000000000;
        trafo(0, 10) = -0.93750000000000;
        trafo(0, 12) = 11.25000000000000;
        trafo(0, 14) = -7.50000000000000;
        trafo(0, 21) = -0.31250000000000;
        trafo(0, 23) = 5.62500000000000;
        trafo(0, 25) = -7.50000000000000;
        trafo(0, 27) = 1.00000000000000;
        trafo(1, 2) = 2.86410980934740;
        trafo(1, 7) = 5.72821961869480;
        trafo(1, 9) = -11.45643923738960;
        trafo(1, 16) = 2.86410980934740;
        trafo(1, 18) = -11.45643923738960;
        trafo(1, 20) = 4.58257569495584;
        trafo(2, 4) = 2.86410980934740;
        trafo(2, 11) = 5.72821961869480;
        trafo(2, 13) = -11.45643923738960;
        trafo(2, 22) = 2.86410980934740;
        trafo(2, 24) = -11.45643923738960;
        trafo(2, 26) = 4.58257569495584;
        trafo(3, 0) = 0.45285552331842;
        trafo(3, 3) = 0.45285552331842;
        trafo(3, 5) = -7.24568837309472;
        trafo(3, 10) = -0.45285552331842;
        trafo(3, 14) = 7.24568837309472;
        trafo(3, 21) = -0.45285552331842;
        trafo(3, 23) = 7.24568837309472;
        trafo(3, 25) = -7.24568837309472;
        trafo(4, 1) = 0.90571104663684;
        trafo(4, 6) = 1.81142209327368;
        trafo(4, 8) = -14.49137674618944;
        trafo(4, 15) = 0.90571104663684;
        trafo(4, 17) = -14.49137674618944;
        trafo(4, 19) = 14.49137674618944;
        trafo(5, 2) = -2.71713313991052;
        trafo(5, 7) = 5.43426627982104;
        trafo(5, 9) = 7.24568837309472;
        trafo(5, 16) = 8.15139941973156;
        trafo(5, 18) = -21.73706511928416;
        trafo(6, 4) = -8.15139941973156;
        trafo(6, 11) = -5.43426627982104;
        trafo(6, 13) = 21.73706511928416;
        trafo(6, 22) = 2.71713313991052;
        trafo(6, 24) = -7.24568837309472;
        trafo(7, 0) = -0.49607837082461;
        trafo(7, 3) = 2.48039185412305;
        trafo(7, 5) = 4.96078370824611;
        trafo(7, 10) = 2.48039185412305;
        trafo(7, 12) = -29.76470224947665;
        trafo(7, 21) = -0.49607837082461;
        trafo(7, 23) = 4.96078370824611;
        trafo(8, 1) = -1.98431348329844;
        trafo(8, 8) = 19.84313483298443;
        trafo(8, 15) = 1.98431348329844;
        trafo(8, 17) = -19.84313483298443;
        trafo(9, 2) = 2.32681380862329;
        trafo(9, 7) = -23.26813808623286;
        trafo(9, 16) = 11.63406904311643;
        trafo(10, 4) = 11.63406904311643;
        trafo(10, 11) = -23.26813808623286;
        trafo(10, 22) = 2.32681380862329;
        trafo(11, 0) = 0.67169328938140;
        trafo(11, 3) = -10.07539934072094;
        trafo(11, 10) = 10.07539934072094;
        trafo(11, 21) = -0.67169328938140;
        trafo(12, 1) = 4.03015973628838;
        trafo(12, 6) = -13.43386578762792;
        trafo(12, 15) = 4.03015973628838;

        return trafo;
    }
    default:
    {
        throw std::runtime_error("Inappropriate angular momentum given!\n");
    }
    }
}

void LI::sphericalTrafo(const arma::dmat &trafo_a, const arma::dmat &trafo_b,
                        const arma::dmat &trafo_c, const arma::dmat &trafo_d,
                        const vec4d &eri4_shells_cart, vec4d &eri4_shells_sph)
{
    vec4d eri4_(trafo_a.n_cols, trafo_b.n_cols, trafo_c.n_cols, trafo_d.n_rows, 0);
    for (size_t mu_ = 0; mu_ < trafo_a.n_cols; mu_++)
        for (size_t nu_ = 0; nu_ < trafo_b.n_cols; nu_++)
            for (size_t kappa_ = 0; kappa_ < trafo_c.n_cols; kappa_++)
                for (size_t tau_ = 0; tau_ < trafo_d.n_cols; tau_++)
                    for (size_t tau = 0; tau < trafo_d.n_rows; tau++)
                        eri4_(mu_, nu_, kappa_, tau) += eri4_shells_cart(mu_, nu_, kappa_, tau_) *
                                                        trafo_d(tau, tau_);

    vec4d eri4__(trafo_a.n_cols, trafo_b.n_cols, trafo_c.n_rows, trafo_d.n_rows, 0);
    for (size_t mu_ = 0; mu_ < trafo_a.n_cols; mu_++)
        for (size_t nu_ = 0; nu_ < trafo_b.n_cols; nu_++)
            for (size_t kappa_ = 0; kappa_ < trafo_c.n_cols; kappa_++)
                for (size_t kappa = 0; kappa < trafo_c.n_rows; kappa++)
                    for (size_t tau = 0; tau < trafo_d.n_rows; tau++)
                        eri4__(mu_, nu_, kappa, tau) += eri4_(mu_, nu_, kappa_, tau) *
                                                        trafo_c(kappa, kappa_);

    vec4d eri4___(trafo_a.n_cols, trafo_b.n_rows, trafo_c.n_rows, trafo_d.n_rows, 0);
    for (size_t mu_ = 0; mu_ < trafo_a.n_cols; mu_++)
        for (size_t nu_ = 0; nu_ < trafo_b.n_cols; nu_++)
            for (size_t nu = 0; nu < trafo_b.n_rows; nu++)
                for (size_t kappa = 0; kappa < trafo_c.n_rows; kappa++)
                    for (size_t tau = 0; tau < trafo_d.n_rows; tau++)
                        eri4___(mu_, nu, kappa, tau) += eri4__(mu_, nu_, kappa, tau) *
                                                        trafo_b(nu, nu_);

    for (size_t mu_ = 0; mu_ < trafo_a.n_cols; mu_++)
        for (size_t mu = 0; mu < trafo_a.n_rows; mu++)
            for (size_t nu = 0; nu < trafo_b.n_rows; nu++)
                for (size_t kappa = 0; kappa < trafo_c.n_rows; kappa++)
                    for (size_t tau = 0; tau < trafo_d.n_rows; tau++)
                        eri4_shells_sph(mu, nu, kappa, tau) += eri4___(mu_, nu, kappa, tau) *
                                                               trafo_a(mu, mu_);
}

void LI::transferIntegrals(const size_t ipair,
                           const ShellPairData &shell_pair_data,
                           const arma::dmat &ints_sph,
                           vec2d &ints)
{
    auto [pos_a, pos_b] = shell_pair_data.offsets[ipair];
    auto [norms_a, norms_b] = shell_pair_data.norms[ipair];
    for (size_t mu = 0; mu < ints_sph.n_rows; mu++)
        for (size_t nu = 0; nu < ints_sph.n_cols; nu++)
        {
            double normalized_int = ints_sph(mu, nu) * norms_a[mu] * norms_b[nu];
            ints(pos_a + mu, pos_b + nu) = normalized_int;
            ints(pos_b + nu, pos_a + mu) = normalized_int;
        }
}

void LI::transferIntegrals(const size_t ipair_ab, const size_t ipair_cd,
                           const ShellPairData &shell_pair_data_ab,
                           const ShellPairData &shell_pair_data_cd,
                           const vec4d &eri4_shells_sph, vec4d &eri4)
{
    auto [pos_a, pos_b] = shell_pair_data_ab.offsets[ipair_ab];
    auto [pos_c, pos_d] = shell_pair_data_cd.offsets[ipair_cd];

    auto [norms_a, norms_b] = shell_pair_data_ab.norms[ipair_ab];
    auto [norms_c, norms_d] = shell_pair_data_cd.norms[ipair_cd];

    size_t dim_a = eri4_shells_sph.getDim(0);
    size_t dim_b = eri4_shells_sph.getDim(1);
    size_t dim_c = eri4_shells_sph.getDim(2);
    size_t dim_d = eri4_shells_sph.getDim(3);

    for (size_t mu = 0; mu < dim_a; mu++)
        for (size_t nu = 0; nu < dim_b; nu++)
            for (size_t kappa = 0; kappa < dim_c; kappa++)
                for (size_t tau = 0; tau < dim_d; tau++)
                {
                    double normalized_int = eri4_shells_sph(mu, nu, kappa, tau) *
                                            norms_a[mu] * norms_b[nu] *
                                            norms_c[kappa] * norms_d[tau];

                    size_t a = pos_a + mu;
                    size_t b = pos_b + nu;
                    size_t c = pos_c + kappa;
                    size_t d = pos_d + tau;

                    eri4(a, b, c, d) = normalized_int;
                    eri4(a, b, d, c) = normalized_int;
                    eri4(b, a, c, d) = normalized_int;
                    eri4(b, a, d, c) = normalized_int;
                    eri4(c, d, a, b) = normalized_int;
                    eri4(c, d, b, a) = normalized_int;
                    eri4(d, c, a, b) = normalized_int;
                    eri4(d, c, b, a) = normalized_int;
                }
}

void LI::transferIntegrals(const size_t ipair_ab, const size_t ipair_cd,
                           const ShellPairData &shell_pair_data_ab,
                           const ShellPairData &shell_pair_data_cd,
                           const arma::dmat &eri4_shells_sph, vec4d &eri4)
{
    auto [pos_a, pos_b] = shell_pair_data_ab.offsets[ipair_ab];
    auto [pos_c, pos_d] = shell_pair_data_cd.offsets[ipair_cd];

    const auto &[norms_a, norms_b] = shell_pair_data_ab.norms[ipair_ab];
    const auto &[norms_c, norms_d] = shell_pair_data_cd.norms[ipair_cd];

    size_t dim_a = norms_a.size();
    size_t dim_b = norms_b.size();
    size_t dim_c = norms_c.size();
    size_t dim_d = norms_d.size();

    for (size_t mu = 0, munu = 0; mu < dim_a; mu++)
        for (size_t nu = 0; nu < dim_b; nu++, munu++)
            for (size_t ka = 0, kata = 0; ka < dim_c; ka++)
                for (size_t ta = 0; ta < dim_d; ta++, kata++)
                {
                    double normalized_int = eri4_shells_sph(munu, kata) *
                                            norms_a[mu] * norms_b[nu] *
                                            norms_c[ka] * norms_d[ta];

                    size_t a = pos_a + mu;
                    size_t b = pos_b + nu;
                    size_t c = pos_c + ka;
                    size_t d = pos_d + ta;

                    eri4(a, b, c, d) = normalized_int;
                    eri4(a, b, d, c) = normalized_int;
                    eri4(b, a, c, d) = normalized_int;
                    eri4(b, a, d, c) = normalized_int;
                    eri4(c, d, a, b) = normalized_int;
                    eri4(c, d, b, a) = normalized_int;
                    eri4(d, c, a, b) = normalized_int;
                    eri4(d, c, b, a) = normalized_int;
                }

}