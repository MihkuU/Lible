#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/util.hpp>
#include <lible/ints/ints_util.hpp>
#include <lible/ints/spherical_trafo.hpp>

#include <fmt/core.h>

namespace LIT = lible::ints::two;

using std::pair, std::vector;

lible::vec4d LIT::calcERI4(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}", "4-center 2-electron integrals..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> shell_pair_datas(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        shell_pair_datas[ipair] = ShellPairData(la, lb, structure);
    }

    vector<vector<vector<vec4d>>> ecoeffs(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        vector<vector<vec4d>> ecoeffs_ipair;
        MD::calcECoeffs(la, lb, shell_pair_datas[ipair], ecoeffs_ipair);
        ecoeffs[ipair] = std::move(ecoeffs_ipair);
    }

    size_t dim_ao = structure.getDimAO();
    vec4d eri4(dim_ao, 0);
    for (int lalb = l_pairs.size() - 1; lalb >= 0; lalb--)
        for (int lcld = lalb; lcld >= 0; lcld--)
        {
            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            int lab = la + lb;
            int lcd = lc + ld;
            int labcd = lab + lcd;

            int dim_a_cart = dimCartesians(la);
            int dim_b_cart = dimCartesians(lb);
            int dim_c_cart = dimCartesians(lc);
            int dim_d_cart = dimCartesians(ld);

            int dim_a_sph = dimSphericals(la);
            int dim_b_sph = dimSphericals(lb);
            int dim_c_sph = dimSphericals(lc);
            int dim_d_sph = dimSphericals(ld);

            BoysF boys_f(labcd);

            auto cart_exps_a = cart_exps[la];
            auto cart_exps_b = cart_exps[lb];
            auto cart_exps_c = cart_exps[lc];
            auto cart_exps_d = cart_exps[ld];

            const auto &ecoeffs_lalb = ecoeffs[lalb];
            const auto &ecoeffs_lcld = ecoeffs[lcld];

            const auto &shell_pair_data_ab = shell_pair_datas[lalb];
            const auto &shell_pair_data_cd = shell_pair_datas[lcld];

            arma::dmat sph_trafo_a = returnSphericalTrafo(la);
            arma::dmat sph_trafo_b = returnSphericalTrafo(lb);
            arma::dmat sph_trafo_c = returnSphericalTrafo(lc);
            arma::dmat sph_trafo_d = returnSphericalTrafo(ld);

            size_t n_pairs_ab = shell_pair_data_ab.n_pairs;
            size_t n_pairs_cd = shell_pair_data_cd.n_pairs;

            // size_t n_shells_abcd = 0;
            if (lalb == lcld)
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        vec4d eri4_shells_cart(dim_a_cart, dim_b_cart, dim_c_cart, dim_d_cart, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_lalb, ecoeffs_lcld,
                                   cart_exps_a, cart_exps_b, cart_exps_c, cart_exps_d,
                                   shell_pair_data_ab, shell_pair_data_cd, boys_f,
                                   eri4_shells_cart);

                        vec4d eri4_shells_sph(dim_a_sph, dim_b_sph, dim_c_sph, dim_d_sph, 0);
                        sphericalTrafo(sph_trafo_a, sph_trafo_b, sph_trafo_c, sph_trafo_d,
                                       eri4_shells_cart, eri4_shells_sph);

                        transferIntegrals(ipair_ab, ipair_cd, shell_pair_data_ab,
                                          shell_pair_data_cd, eri4_shells_sph, eri4);
                    }
            else
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vec4d eri4_shells_cart(dim_a_cart, dim_b_cart, dim_c_cart, dim_d_cart, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_lalb, ecoeffs_lcld,
                                   cart_exps_a, cart_exps_b, cart_exps_c, cart_exps_d,
                                   shell_pair_data_ab, shell_pair_data_cd, boys_f,
                                   eri4_shells_cart);

                        vec4d eri4_shells_sph(dim_a_sph, dim_b_sph, dim_c_sph, dim_d_sph, 0);
                        sphericalTrafo(sph_trafo_a, sph_trafo_b, sph_trafo_c, sph_trafo_d,
                                       eri4_shells_cart, eri4_shells_sph);

                        transferIntegrals(ipair_ab, ipair_cd, shell_pair_data_ab,
                                          shell_pair_data_cd, eri4_shells_sph, eri4);
                    }
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));

    return eri4;
}

void LIT::calcERI4Benchmark(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}\n", "ERI4 benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> shell_pair_datas(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        shell_pair_datas[ipair] = ShellPairData(la, lb, structure);
    }

    auto t0_ecoeffs{std::chrono::steady_clock::now()};
    vector<vector<vector<vec4d>>> ecoeffs(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        vector<vector<vec4d>> ecoeffs_ipair;
        MD::calcECoeffs(la, lb, shell_pair_datas[ipair], ecoeffs_ipair);
        ecoeffs[ipair] = std::move(ecoeffs_ipair);
    }  

    size_t dim_ao = structure.getDimAO();
    for (int lalb = l_pairs.size() - 1; lalb >= 0; lalb--)
        for (int lcld = lalb; lcld >= 0; lcld--)
        {
            auto start{std::chrono::steady_clock::now()};

            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            int lab = la + lb;
            int lcd = lc + ld;
            int labcd = lab + lcd;

            int dim_a_cart = dimCartesians(la);
            int dim_b_cart = dimCartesians(lb);
            int dim_c_cart = dimCartesians(lc);
            int dim_d_cart = dimCartesians(ld);

            int dim_a_sph = dimSphericals(la);
            int dim_b_sph = dimSphericals(lb);
            int dim_c_sph = dimSphericals(lc);
            int dim_d_sph = dimSphericals(ld);

            BoysF boys_f(labcd);

            auto cart_exps_a = cart_exps[la];
            auto cart_exps_b = cart_exps[lb];
            auto cart_exps_c = cart_exps[lc];
            auto cart_exps_d = cart_exps[ld];

            const auto &ecoeffs_lalb = ecoeffs[lalb];
            const auto &ecoeffs_lcld = ecoeffs[lcld];

            const auto &shell_pair_data_ab = shell_pair_datas[lalb];
            const auto &shell_pair_data_cd = shell_pair_datas[lcld];

            arma::dmat sph_trafo_a = returnSphericalTrafo(la);
            arma::dmat sph_trafo_b = returnSphericalTrafo(lb);
            arma::dmat sph_trafo_c = returnSphericalTrafo(lc);
            arma::dmat sph_trafo_d = returnSphericalTrafo(ld);

            size_t n_pairs_ab = shell_pair_data_ab.n_pairs;
            size_t n_pairs_cd = shell_pair_data_cd.n_pairs;

            size_t n_shells_abcd = 0;
            if (lalb == lcld)
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        vec4d eri4_shells_cart(dim_a_cart, dim_b_cart, dim_c_cart, dim_d_cart, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_lalb, ecoeffs_lcld,
                                   cart_exps_a, cart_exps_b, cart_exps_c, cart_exps_d,
                                   shell_pair_data_ab, shell_pair_data_cd, boys_f,
                                   eri4_shells_cart);

                        vec4d eri4_shells_sph(dim_a_sph, dim_b_sph, dim_c_sph, dim_d_sph, 0);
                        sphericalTrafo(sph_trafo_a, sph_trafo_b, sph_trafo_c, sph_trafo_d,
                                       eri4_shells_cart, eri4_shells_sph);

                        n_shells_abcd += 4;
                    }
            else
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vec4d eri4_shells_cart(dim_a_cart, dim_b_cart, dim_c_cart, dim_d_cart, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_lalb, ecoeffs_lcld,
                                   cart_exps_a, cart_exps_b, cart_exps_c, cart_exps_d,
                                   shell_pair_data_ab, shell_pair_data_cd, boys_f,
                                   eri4_shells_cart);

                        vec4d eri4_shells_sph(dim_a_sph, dim_b_sph, dim_c_sph, dim_d_sph, 0);
                        sphericalTrafo(sph_trafo_a, sph_trafo_b, sph_trafo_c, sph_trafo_d,
                                       eri4_shells_cart, eri4_shells_sph);

                        n_shells_abcd += 4;
                    }

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            palPrint(fmt::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
                                 n_shells_abcd, duration.count()));
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format("done {:.2e} s\n", duration.count()));
}

void LIT::kernelERI4(const int lab, const int lcd,
                     const size_t ipair_ab, const size_t ipair_cd,
                     const vector<vector<vec4d>> &ecoeffs_lalb,
                     const vector<vector<vec4d>> &ecoeffs_lcld,
                     const vector<CartExps> &cart_exps_a,
                     const vector<CartExps> &cart_exps_b,
                     const vector<CartExps> &cart_exps_c,
                     const vector<CartExps> &cart_exps_d,
                     const ShellPairData &shell_pair_data_ab,
                     const ShellPairData &shell_pair_data_cd,
                     const BoysF &boys_f, vec4d &eri4_shells_cart)
{
    int labcd = lab + lcd;

    const auto &[exps_a, exps_b] = shell_pair_data_ab.exps[ipair_ab];
    const auto &[exps_c, exps_d] = shell_pair_data_cd.exps[ipair_cd];
    const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data_ab.ccoeffs[ipair_ab];
    const auto &[ccoeffs_c, ccoeffs_d] = shell_pair_data_cd.ccoeffs[ipair_cd];
    const auto &[xyz_a, xyz_b] = shell_pair_data_ab.coords[ipair_ab];
    const auto &[xyz_c, xyz_d] = shell_pair_data_cd.coords[ipair_cd];

    const vector<vec4d> &Exyz_ab = ecoeffs_lalb[ipair_ab];
    const vector<vec4d> &Exyz_cd = ecoeffs_lcld[ipair_cd];

    size_t ka = exps_a.size();
    size_t kb = exps_b.size();
    size_t kc = exps_c.size();
    size_t kd = exps_d.size();

    arma::vec::fixed<3> A{xyz_a[0], xyz_a[1], xyz_a[2]};
    arma::vec::fixed<3> B{xyz_b[0], xyz_b[1], xyz_b[2]};
    arma::vec::fixed<3> C{xyz_c[0], xyz_c[1], xyz_c[2]};
    arma::vec::fixed<3> D{xyz_d[0], xyz_d[1], xyz_d[2]};

    vector<double> fnx(labcd + 1, 0);
    vec3d rints(labcd + 1, 0);
    vec4d rints_tmp(labcd + 1, 0);

    for (size_t ia = 0, iab = 0; ia < ka; ia++)
        for (size_t ib = 0; ib < kb; ib++, iab++)
            for (size_t ic = 0, icd = 0; ic < kc; ic++)
                for (size_t id = 0; id < kd; id++, icd++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double c = exps_c[ic];
                    double d = exps_d[id];

                    double p = a + b;
                    double q = c + d;
                    double alpha = p * q / (p + q);

                    arma::vec::fixed<3> P = (a * A + b * B) / p;
                    arma::vec::fixed<3> Q = (c * C + d * D) / q;
                    arma::vec::fixed<3> RPQ = P - Q;

                    double x = alpha * arma::dot(RPQ, RPQ);

                    boys_f.calcFnx(labcd, x, fnx);

                    MD::calcRInts(lab, lcd, alpha, RPQ, fnx, rints_tmp, rints);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q))) *
                                 ccoeffs_a[ia] * ccoeffs_b[ib] * ccoeffs_c[ic] * ccoeffs_d[id];

                    for (const auto [i, j, k, mu] : cart_exps_a)
                        for (const auto [i_, j_, k_, nu] : cart_exps_b)
                            for (const auto [l, m, n, kappa] : cart_exps_c)
                                for (const auto [l_, m_, n_, tau] : cart_exps_d)
                                    for (int t_ = 0; t_ <= l + l_; t_++)
                                        for (int u_ = 0; u_ <= m + m_; u_++)
                                            for (int v_ = 0; v_ <= n + n_; v_++)
                                            {
                                                double sign = std::pow(-1.0, t_ + u_ + v_);
                                                for (int t = 0; t <= i + i_; t++)
                                                    for (int u = 0; u <= j + j_; u++)
                                                        for (int v = 0; v <= k + k_; v++)
                                                        {
                                                            eri4_shells_cart(mu, nu, kappa, tau) +=
                                                                fac * sign *
                                                                Exyz_ab[iab](0, i, i_, t) *
                                                                Exyz_ab[iab](1, j, j_, u) *
                                                                Exyz_ab[iab](2, k, k_, v) *
                                                                Exyz_cd[icd](0, l, l_, t_) *
                                                                Exyz_cd[icd](1, m, m_, u_) *
                                                                Exyz_cd[icd](2, n, n_, v_) *
                                                                rints(t + t_, u + u_, v + v_);
                                                        }
                                            }
                }
}