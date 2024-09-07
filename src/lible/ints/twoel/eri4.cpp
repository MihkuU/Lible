#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/util.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <fmt/core.h>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LIT = lible::ints::two;

using std::pair, std::vector;

namespace lible::ints::two
{
    void kernelERI4(const int lab, const int lcd, const int ipair_ab, const int ipair_cd,
                    const vector<double> &ecoeffs_lalb, const vector<double> &ecoeffs_lcld_tsp,
                    const vector<IdxsTUV> &idxs_tuv_ab, const vector<IdxsTUV> &idxs_tuv_cd_tsp,
                    const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                    const BoysF &boys_f, vector<double> &eri4_shells_sph, vector<double> &rints,
                    vector<double> &fnx, vec4d &rints_tmp)
    {
        int labcd = lab + lcd;

        int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
        int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
        int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
        int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
        int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
        int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
        int pos_c = sp_data_cd.coffsets[2 * ipair_cd];
        int pos_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

        arma::vec::fixed<3> xyz_a{sp_data_ab.coords[6 * ipair_ab],
                                  sp_data_ab.coords[6 * ipair_ab + 1],
                                  sp_data_ab.coords[6 * ipair_ab + 2]};

        arma::vec::fixed<3> xyz_b{sp_data_ab.coords[6 * ipair_ab + 3],
                                  sp_data_ab.coords[6 * ipair_ab + 4],
                                  sp_data_ab.coords[6 * ipair_ab + 5]};

        arma::vec::fixed<3> xyz_c{sp_data_cd.coords[6 * ipair_cd],
                                  sp_data_cd.coords[6 * ipair_cd + 1],
                                  sp_data_cd.coords[6 * ipair_cd + 2]};

        arma::vec::fixed<3> xyz_d{sp_data_cd.coords[6 * ipair_cd + 3],
                                  sp_data_cd.coords[6 * ipair_cd + 4],
                                  sp_data_cd.coords[6 * ipair_cd + 5]};

        int dim_a = dimSphericals(sp_data_ab.la);
        int dim_b = dimSphericals(sp_data_ab.lb);
        int dim_c = dimSphericals(sp_data_cd.la);
        int dim_d = dimSphericals(sp_data_cd.lb);
        int dim_tuv_ab = dimHermiteGaussians(lab);
        int dim_tuv_cd = dimHermiteGaussians(lcd);
        int dim_sph_ab = dim_a * dim_b;
        int dim_sph_cd = dim_c * dim_d;

        int dim_E_ab = dim_sph_ab * dim_tuv_ab;
        int dim_E_cd = dim_sph_cd * dim_tuv_cd;
        int dim_R_x_E = dim_sph_cd * dim_tuv_ab;
        vector<double> X(cdepth_a * cdepth_b * dim_R_x_E, 0);

        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++)
            {
                int pos_X = iab * dim_R_x_E;

                for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                    for (int id = 0; id < cdepth_d; id++)
                    {
                        double a = sp_data_ab.exps[pos_a + ia];
                        double b = sp_data_ab.exps[pos_b + ib];
                        double c = sp_data_cd.exps[pos_c + ic];
                        double d = sp_data_cd.exps[pos_d + id];

                        double p = a + b;
                        double q = c + d;
                        double alpha = p * q / (p + q);

                        arma::vec::fixed<3> xyz_p = (a * xyz_a + b * xyz_b) / p;
                        arma::vec::fixed<3> xyz_q = (c * xyz_c + d * xyz_d) / q;
                        arma::vec::fixed<3> xyz_pq = xyz_p - xyz_q;

                        double x = alpha * arma::dot(xyz_pq, xyz_pq);
                        boys_f.calcFnx(labcd, x, fnx);

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                        calcRInts(lab, lcd, alpha, fac, xyz_pq, fnx, idxs_tuv_ab,
                                  idxs_tuv_cd_tsp, rints_tmp, rints);

                        int pos_ecoeffs_cd = sp_data_cd.offsets_ecoeffs[ipair_cd] + icd * dim_E_cd;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_tuv_ab,
                                    dim_sph_cd, dim_tuv_cd, 1.0, &rints[0], dim_tuv_cd,
                                    &ecoeffs_lcld_tsp[pos_ecoeffs_cd], dim_sph_cd, 1.0, &X[pos_X],
                                    dim_sph_cd);
                        icd++;
                    }
                iab++;
            }

        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
            {
                int pos_X = iab * dim_R_x_E;
                int pos_ecoeffs_ab = sp_data_ab.offsets_ecoeffs[ipair_ab] + iab * dim_E_ab;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_ab, dim_sph_cd,
                            dim_tuv_ab, 1.0, &ecoeffs_lalb[pos_ecoeffs_ab], dim_tuv_ab, &X[pos_X],
                            dim_sph_cd, 1.0, &eri4_shells_sph[0], dim_sph_cd);
            }
    }
}

lible::vec4d LIT::calcERI4(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}", "SHARK ERI4 flat..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> sp_datas;
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        sp_datas.emplace_back(constructShellPairData(la, lb, structure));
    }

    vector<vector<double>> ecoeffs(l_pairs.size());
    vector<vector<double>> ecoeffs_tsp(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        int lab = la + lb;
        int n_ecoeffs_sph = dimSphericals(la) * dimSphericals(lb) * dimHermiteGaussians(lab) *
                            sp_datas[ipair].n_prim_pairs;

        vector<double> ecoeffs_ipair(n_ecoeffs_sph);
        vector<double> ecoeffs_tsp_ipair(n_ecoeffs_sph);

        calcECoeffsSpherical(la, lb, sp_datas[ipair], ecoeffs_ipair, ecoeffs_tsp_ipair);

        ecoeffs[ipair] = std::move(ecoeffs_ipair);
        ecoeffs_tsp[ipair] = std::move(ecoeffs_tsp_ipair);
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

            int dim_a_sph = dimSphericals(la);
            int dim_b_sph = dimSphericals(lb);
            int dim_c_sph = dimSphericals(lc);
            int dim_d_sph = dimSphericals(ld);

            int dim_ab_sph = dim_a_sph * dim_b_sph;
            int dim_cd_sph = dim_c_sph * dim_d_sph;

            BoysF boys_f(labcd);

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sp_data_cd = sp_datas[lcld];

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            vector<IdxsTUV> idxs_tuv_ab = returnIdxsTUV(lab);
            vector<IdxsTUV> idxs_tuv_cd = returnIdxsTUV(lcd);

            int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
            int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

            vector<double> rints(dim_tuv_ab * dim_tuv_cd, 0);
            vector<double> fnx(labcd + 1, 0);
            vec4d rints_tmp(labcd + 1, 0);

            if (lalb == lcld)
                for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (int ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
                                   idxs_tuv_ab, idxs_tuv_cd, sp_data_ab, sp_data_cd, boys_f,
                                   eri4_shells_sph, rints, fnx, rints_tmp);

                        transferIntegrals(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                          eri4_shells_sph, eri4);
                    }
            else
                for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (int ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
                                   idxs_tuv_ab, idxs_tuv_cd, sp_data_ab, sp_data_cd, boys_f,
                                   eri4_shells_sph, rints, fnx, rints_tmp);

                        transferIntegrals(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                          eri4_shells_sph, eri4);
                    }
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));

    return eri4;
}

void LIT::calcERI4Benchmark(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}\n", "ERI4 (Shark flat) benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> sp_datas;
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        sp_datas.emplace_back(constructShellPairData(la, lb, structure));
    }

    vector<vector<double>> ecoeffs(l_pairs.size());
    vector<vector<double>> ecoeffs_tsp(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        int lab = la + lb;
        int n_ecoeffs_sph = dimSphericals(la) * dimSphericals(lb) * dimHermiteGaussians(lab) *
                            sp_datas[ipair].n_prim_pairs;

        vector<double> ecoeffs_ipair(n_ecoeffs_sph);
        vector<double> ecoeffs_tsp_ipair(n_ecoeffs_sph);

        calcECoeffsSpherical(la, lb, sp_datas[ipair], ecoeffs_ipair, ecoeffs_tsp_ipair);

        ecoeffs[ipair] = std::move(ecoeffs_ipair);
        ecoeffs_tsp[ipair] = std::move(ecoeffs_tsp_ipair);
    }

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

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sp_data_cd = sp_datas[lcld];

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            vector<IdxsTUV> idxs_tuv_ab = returnIdxsTUV(lab);
            vector<IdxsTUV> idxs_tuv_cd = returnIdxsTUV(lcd);

            int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
            int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

            vector<double> rints(dim_tuv_ab * dim_tuv_cd, 0);
            vector<double> fnx(labcd + 1, 0);
            vec4d rints_tmp(labcd + 1, 0);

            size_t n_shells_abcd = 0;
            if (lalb == lcld)
                for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (int ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
                                   idxs_tuv_ab, idxs_tuv_cd, sp_data_ab, sp_data_cd, boys_f,
                                   eri4_shells_sph, rints, fnx, rints_tmp);
                    }
            else
                for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (int ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
                                   idxs_tuv_ab, idxs_tuv_cd, sp_data_ab, sp_data_cd, boys_f,
                                   eri4_shells_sph, rints, fnx, rints_tmp);
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