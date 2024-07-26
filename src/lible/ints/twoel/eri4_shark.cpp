#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/util.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/util.hpp>

#include <fmt/core.h>

namespace LIT = lible::ints::two;

using std::pair, std::vector;

lible::vec4d LIT::calcERI4Shark(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}", "SHARK ERI4..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    auto start_spdata{std::chrono::steady_clock::now()};
    vector<ShellPairData> shell_pair_datas(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        shell_pair_datas[ipair] = ShellPairData(la, lb, structure);
    }
    auto stop_spdata{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration_spdata{stop_spdata - start_spdata};
    palPrint(fmt::format("\nduration_spdata: {:.2e} s\n", duration_spdata.count()));

    auto start_ecoeffs{std::chrono::steady_clock::now()};
    vector<vector<vector<arma::dmat>>> ecoeffs(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        vector<vector<arma::dmat>> ecoeffs_ipair;
        calcECoeffsSpherical(la, lb, shell_pair_datas[ipair], ecoeffs_ipair);
        ecoeffs[ipair] = std::move(ecoeffs_ipair);
    }
    auto stop_ecoeffs{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration_ecoeffs{stop_ecoeffs - start_ecoeffs};
    palPrint(fmt::format("duration_ecoeffs: {:.2e} s\n", duration_ecoeffs.count()));

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

            int dim_a_sph = dimSphericals(la);
            int dim_b_sph = dimSphericals(lb);
            int dim_c_sph = dimSphericals(lc);
            int dim_d_sph = dimSphericals(ld);

            int dim_ab_sph = dim_a_sph * dim_b_sph;
            int dim_cd_sph = dim_c_sph * dim_d_sph;

            BoysF boys_f(labcd);

            const auto &ecoeffs_ab = ecoeffs[lalb];
            const auto &ecoeffs_cd = ecoeffs[lcld];

            const auto &shell_pair_data_ab = shell_pair_datas[lalb];
            const auto &shell_pair_data_cd = shell_pair_datas[lcld];

            size_t n_pairs_ab = shell_pair_data_ab.n_pairs;
            size_t n_pairs_cd = shell_pair_data_cd.n_pairs;

            vector<IdxsTUV> idxs_tuv_ab = returnIdxsTUV(lab);
            vector<IdxsTUV> idxs_tuv_cd = returnIdxsTUV(lcd);

            if (lalb == lcld)
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        arma::dmat eri4_shells_sph(dim_ab_sph, dim_cd_sph, arma::fill::zeros);
                        kernelERI4Shark(lab, lcd, ipair_ab, ipair_cd,
                                        ecoeffs_ab, ecoeffs_cd,
                                        idxs_tuv_ab, idxs_tuv_cd,
                                        shell_pair_data_ab, shell_pair_data_cd,
                                        boys_f, eri4_shells_sph);

                        transferIntegrals(ipair_ab, ipair_cd, shell_pair_data_ab,
                                          shell_pair_data_cd, eri4_shells_sph, eri4);
                    }
            else
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        arma::dmat eri4_shells_sph(dim_ab_sph, dim_cd_sph, arma::fill::zeros);
                        kernelERI4Shark(lab, lcd, ipair_ab, ipair_cd,
                                        ecoeffs_ab, ecoeffs_cd,
                                        idxs_tuv_ab, idxs_tuv_cd,
                                        shell_pair_data_ab, shell_pair_data_cd,
                                        boys_f, eri4_shells_sph);

                        transferIntegrals(ipair_ab, ipair_cd, shell_pair_data_ab,
                                          shell_pair_data_cd, eri4_shells_sph, eri4);
                    }
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));

    return eri4;
}

void LIT::calcERI4BenchmarkShark(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}\n", "ERI4 (Shark) benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> shell_pair_datas(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        shell_pair_datas[ipair] = ShellPairData(la, lb, structure);
    }

    vector<vector<vector<arma::dmat>>> ecoeffs(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        vector<vector<arma::dmat>> ecoeffs_ipair;
        calcECoeffsSpherical(la, lb, shell_pair_datas[ipair], ecoeffs_ipair);
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

            int dim_a_sph = dimSphericals(la);
            int dim_b_sph = dimSphericals(lb);
            int dim_c_sph = dimSphericals(lc);
            int dim_d_sph = dimSphericals(ld);

            int dim_ab_sph = dim_a_sph * dim_b_sph;
            int dim_cd_sph = dim_c_sph * dim_d_sph;

            BoysF boys_f(labcd);

            const auto &ecoeffs_ab = ecoeffs[lalb];
            const auto &ecoeffs_cd = ecoeffs[lcld];

            const auto &shell_pair_data_ab = shell_pair_datas[lalb];
            const auto &shell_pair_data_cd = shell_pair_datas[lcld];

            size_t n_pairs_ab = shell_pair_data_ab.n_pairs;
            size_t n_pairs_cd = shell_pair_data_cd.n_pairs;

            vector<IdxsTUV> idxs_tuv_ab = returnIdxsTUV(lab);
            vector<IdxsTUV> idxs_tuv_cd = returnIdxsTUV(lcd);

            size_t n_shells_abcd = 0;
            if (lalb == lcld)
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        arma::dmat eri4_shells_sph(dim_ab_sph, dim_cd_sph, arma::fill::zeros);
                        kernelERI4Shark(lab, lcd, ipair_ab, ipair_cd,
                                        ecoeffs_ab, ecoeffs_cd,
                                        idxs_tuv_ab, idxs_tuv_cd,
                                        shell_pair_data_ab, shell_pair_data_cd,
                                        boys_f, eri4_shells_sph);

                        n_shells_abcd += 4;
                    }
            else
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        arma::dmat eri4_shells_sph(dim_ab_sph, dim_cd_sph, arma::fill::zeros);
                        kernelERI4Shark(lab, lcd, ipair_ab, ipair_cd,
                                        ecoeffs_ab, ecoeffs_cd,
                                        idxs_tuv_ab, idxs_tuv_cd,
                                        shell_pair_data_ab, shell_pair_data_cd,
                                        boys_f, eri4_shells_sph);

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

void LIT::kernelERI4Shark(const int lab, const int lcd,
                          const size_t ipair_ab, const size_t ipair_cd,
                          const vector<vector<arma::dmat>> &ecoeffs_lalb,
                          const vector<vector<arma::dmat>> &ecoeffs_lcld,
                          const vector<IdxsTUV> &idxs_tuv_ab,
                          const vector<IdxsTUV> &idxs_tuv_cd,
                          const ShellPairData &shell_pair_data_ab,
                          const ShellPairData &shell_pair_data_cd,
                          const BoysF &boys_f, arma::dmat &eri4_shells_sph)
{
    int labcd = lab + lcd;

    const auto &[exps_a, exps_b] = shell_pair_data_ab.exps[ipair_ab];
    const auto &[exps_c, exps_d] = shell_pair_data_cd.exps[ipair_cd];
    const auto &[xyz_a, xyz_b] = shell_pair_data_ab.coords[ipair_ab];
    const auto &[xyz_c, xyz_d] = shell_pair_data_cd.coords[ipair_cd];

    const vector<arma::dmat> &Exyz_ab = ecoeffs_lalb[ipair_ab];
    const vector<arma::dmat> &Exyz_cd = ecoeffs_lcld[ipair_cd];

    size_t ka = exps_a.size();
    size_t kb = exps_b.size();
    size_t kc = exps_c.size();
    size_t kd = exps_d.size();

    arma::vec::fixed<3> A{xyz_a[0], xyz_a[1], xyz_a[2]};
    arma::vec::fixed<3> B{xyz_b[0], xyz_b[1], xyz_b[2]};
    arma::vec::fixed<3> C{xyz_c[0], xyz_c[1], xyz_c[2]};
    arma::vec::fixed<3> D{xyz_d[0], xyz_d[1], xyz_d[2]};

    vector<double> fnx(labcd + 1, 0);
    vec4d rints_tmp(labcd + 1, 0);
    vec3d rints_3d_tmp(labcd + 1, 0);

    vec3i tuv_poss_ab = returnTUVPoss(lab);
    vec3i tuv_poss_cd = returnTUVPoss(lcd);

    int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
    int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;
    arma::dmat rints(dim_tuv_ab, dim_tuv_cd, arma::fill::zeros);

    int dim_sph_cd = Exyz_cd[0].n_rows;
    vector<arma::dmat> X(ka * kb, arma::zeros(dim_tuv_ab, dim_sph_cd));

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

                    calcRInts(lab, lcd, alpha, RPQ, fnx, idxs_tuv_ab, idxs_tuv_cd,
                              rints_tmp, rints);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));

                    X[iab] += fac * rints * Exyz_cd[icd].t();
                }

    for (size_t ia = 0, iab = 0; ia < ka; ia++)
        for (size_t ib = 0; ib < kb; ib++, iab++)
            eri4_shells_sph += Exyz_ab[iab] * X[iab];
}