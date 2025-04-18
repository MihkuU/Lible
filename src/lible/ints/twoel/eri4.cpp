#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/utils.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <map>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LIT = lible::ints::two;

using std::array, std::map, std::pair, std::vector;

namespace lible::ints::two
{
    void kernelERI4(const int lab, const int lcd, const int ipair_ab, const int ipair_cd,
                    const vector<double> &ecoeffs_ab, const vector<double> &ecoeffs_cd_tsp,
                    const vector<array<int, 3>> &idxs_tuv_ab,
                    const vector<array<int, 3>> &idxs_tuv_cd,
                    const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                    const BoysF &boys_f, vector<double> &eri4_shells_sph, vector<double> &rints,
                    vector<double> &fnx, vec4d &rints_tmp)
    {
        int labcd = lab + lcd;

        int dim_a = sp_data_ab.cdepths[2 * ipair_ab];
        int dim_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
        int dim_c = sp_data_cd.cdepths[2 * ipair_cd];
        int dim_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
        int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
        int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
        int pos_c = sp_data_cd.coffsets[2 * ipair_cd];
        int pos_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

        array<double, 3> xyz_a{sp_data_ab.coords[6 * ipair_ab],
                               sp_data_ab.coords[6 * ipair_ab + 1],
                               sp_data_ab.coords[6 * ipair_ab + 2]};

        array<double, 3> xyz_b{sp_data_ab.coords[6 * ipair_ab + 3],
                               sp_data_ab.coords[6 * ipair_ab + 4],
                               sp_data_ab.coords[6 * ipair_ab + 5]};

        array<double, 3> xyz_c{sp_data_cd.coords[6 * ipair_cd],
                               sp_data_cd.coords[6 * ipair_cd + 1],
                               sp_data_cd.coords[6 * ipair_cd + 2]};

        array<double, 3> xyz_d{sp_data_cd.coords[6 * ipair_cd + 3],
                               sp_data_cd.coords[6 * ipair_cd + 4],
                               sp_data_cd.coords[6 * ipair_cd + 5]};

        int dim_a_sph = numSphericals(sp_data_ab.la);
        int dim_b_sph = numSphericals(sp_data_ab.lb);
        int dim_c_sph = numSphericals(sp_data_cd.la);
        int dim_d_sph = numSphericals(sp_data_cd.lb);
        int dim_tuv_ab = numHermites(lab);
        int dim_tuv_cd = numHermites(lcd);
        int dim_sph_ab = dim_a_sph * dim_b_sph;
        int dim_sph_cd = dim_c_sph * dim_d_sph;

        int dim_ecoeffs_ab = dim_sph_ab * dim_tuv_ab;
        int dim_ecoeffs_cd = dim_sph_cd * dim_tuv_cd;
        int dim_rints_x_ecoeffs = dim_sph_cd * dim_tuv_ab;
        vector<double> rints_x_ecoeffs(dim_a * dim_b * dim_rints_x_ecoeffs, 0);

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                int pos_rints_x_ecoeffs = iab * dim_rints_x_ecoeffs;
                for (int ic = 0, icd = 0; ic < dim_c; ic++)
                    for (int id = 0; id < dim_d; id++, icd++)
                    {
                        double a = sp_data_ab.exps[pos_a + ia];
                        double b = sp_data_ab.exps[pos_b + ib];
                        double c = sp_data_cd.exps[pos_c + ic];
                        double d = sp_data_cd.exps[pos_d + id];

                        double p = a + b;
                        double q = c + d;
                        double alpha = p * q / (p + q);

                        array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                               (a * xyz_a[1] + b * xyz_b[1]) / p,
                                               (a * xyz_a[2] + b * xyz_b[2]) / p};

                        array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                               (c * xyz_c[1] + d * xyz_d[1]) / q,
                                               (c * xyz_c[2] + d * xyz_d[2]) / q};

                        array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],
                                                xyz_p[2] - xyz_q[2]};

                        double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                        double x = alpha * (xx * xx + xy * xy + xz * xz);
                        boys_f.calcFnx(labcd, x, fnx);

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                        calcRInts_(lab, lcd, fac, alpha, xyz_pq, fnx, idxs_tuv_ab, idxs_tuv_cd,
                                   rints_tmp, rints);

                        int pos_ecoeffs_cd = sp_data_cd.offsets_ecoeffs[ipair_cd] +
                                             icd * dim_ecoeffs_cd;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_tuv_ab,
                                    dim_sph_cd, dim_tuv_cd, 1.0, &rints[0], dim_tuv_cd,
                                    &ecoeffs_cd_tsp[pos_ecoeffs_cd], dim_sph_cd, 1.0,
                                    &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_cd);
                    }
            }

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                int pos_rints_x_ecoeffs = iab * dim_rints_x_ecoeffs;
                int pos_ecoeffs_ab = sp_data_ab.offsets_ecoeffs[ipair_ab] + iab * dim_ecoeffs_ab;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_ab, dim_sph_cd,
                            dim_tuv_ab, 1.0, &ecoeffs_ab[pos_ecoeffs_ab], dim_tuv_ab,
                            &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_cd, 1.0,
                            &eri4_shells_sph[0], dim_sph_cd);
            }
    }

    void kernelERI4Diagonal(const int lab, const int ipair_ab, const vector<double> &ecoeffs_ab,
                            const vector<double> &ecoeffs_ab_tsp,
                            const vector<array<int, 3>> &idxs_tuv,
                            const BoysF &boys_f, const ShellPairData &sp_data_ab,
                            vec4d &rints_tmp, vector<double> &eri4_shells_sph, vector<double> &fnx,
                            vector<double> &rints)
    {
        int labab = lab + lab;

        int dim_a = sp_data_ab.cdepths[2 * ipair_ab];
        int dim_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
        int dim_c = dim_a;
        int dim_d = dim_b;
        int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
        int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
        int pos_c = pos_a;
        int pos_d = pos_b;

        array<double, 3> xyz_a{sp_data_ab.coords[6 * ipair_ab],
                               sp_data_ab.coords[6 * ipair_ab + 1],
                               sp_data_ab.coords[6 * ipair_ab + 2]};

        array<double, 3> xyz_b{sp_data_ab.coords[6 * ipair_ab + 3],
                               sp_data_ab.coords[6 * ipair_ab + 4],
                               sp_data_ab.coords[6 * ipair_ab + 5]};

        int dim_a_sph = numSphericals(sp_data_ab.la);
        int dim_b_sph = numSphericals(sp_data_ab.lb);
        int dim_tuv_ab = numHermites(lab);
        int dim_sph_ab = dim_a_sph * dim_b_sph;

        int dim_ecoeffs_ab = dim_sph_ab * dim_tuv_ab;
        int dim_rints_x_ecoeffs = dim_sph_ab * dim_tuv_ab;
        vector<double> rints_x_ecoeffs(dim_a * dim_b * dim_rints_x_ecoeffs, 0);

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++)
            {
                int pos_rints_x_ecoeffs = iab * dim_rints_x_ecoeffs;

                for (int ic = 0, icd = 0; ic < dim_c; ic++)
                    for (int id = 0; id < dim_d; id++)
                    {
                        double a = sp_data_ab.exps[pos_a + ia];
                        double b = sp_data_ab.exps[pos_b + ib];
                        double c = sp_data_ab.exps[pos_c + ic];
                        double d = sp_data_ab.exps[pos_d + id];

                        double p = a + b;
                        double q = c + d;
                        array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                               (a * xyz_a[1] + b * xyz_b[1]) / p,
                                               (a * xyz_a[2] + b * xyz_b[2]) / p};

                        array<double, 3> xyz_q{(c * xyz_a[0] + d * xyz_b[0]) / q,
                                               (c * xyz_a[1] + d * xyz_b[1]) / q,
                                               (c * xyz_a[2] + d * xyz_b[2]) / q};

                        array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],
                                                xyz_p[2] - xyz_q[2]};

                        double alpha = p * q / (p + q);
                        double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                        double x = alpha * (xx * xx + xy * xy + xz * xz);
                        boys_f.calcFnx(labab, x, fnx);

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                        calcRInts_(lab, lab, fac, alpha, xyz_pq, fnx, idxs_tuv, idxs_tuv, rints_tmp, rints);

                        int pos_ecoeffs_cd = sp_data_ab.offsets_ecoeffs[ipair_ab] +
                                             icd * dim_ecoeffs_ab;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_tuv_ab,
                                    dim_sph_ab, dim_tuv_ab, 1.0, &rints[0], dim_tuv_ab,
                                    &ecoeffs_ab_tsp[pos_ecoeffs_cd], dim_sph_ab, 1.0,
                                    &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_ab);

                        icd++;
                    }
                iab++;
            }

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                int pos_rints_x_ecoeffs = iab * dim_rints_x_ecoeffs;
                int pos_ecoeffs_ab = sp_data_ab.offsets_ecoeffs[ipair_ab] + iab * dim_ecoeffs_ab;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_ab, dim_sph_ab,
                            dim_tuv_ab, 1.0, &ecoeffs_ab[pos_ecoeffs_ab], dim_tuv_ab,
                            &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_ab, 1.0,
                            &eri4_shells_sph[0], dim_sph_ab);
            }
    }
}

lible::vec4d LIT::calcERI4(const Structure &structure)
{
    vector<pair<int, int>> l_pairs = returnLPairs(structure.getMaxL());

    vector<ShellPairData> sp_datas = constructShellPairDatas(l_pairs, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(l_pairs, sp_datas);

    size_t dim_ao = structure.getDimAO();
    vec4d eri4(dim_ao, 0);
    for (size_t lalb = 0; lalb < l_pairs.size(); lalb++)
        for (size_t lcld = 0; lcld <= lalb; lcld++)
        {
            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sp_data_cd = sp_datas[lcld];

            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            int lab = la + lb;
            int lcd = lc + ld;

            int dim_a_sph = numSphericals(la);
            int dim_b_sph = numSphericals(lb);
            int dim_c_sph = numSphericals(lc);
            int dim_d_sph = numSphericals(ld);
            int dim_ab_sph = dim_a_sph * dim_b_sph;
            int dim_cd_sph = dim_c_sph * dim_d_sph;
            int dim_tuv_ab = numHermites(lab);
            int dim_tuv_cd = numHermites(lcd);

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            int labcd = lab + lcd;
            BoysF boys_f(labcd);

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

            vector<array<int, 3>> idxs_tuv_ab = returnHermiteGaussianIdxs(lab);
            vector<array<int, 3>> idxs_tuv_cd = returnHermiteGaussianIdxs(lcd);

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

                        transferIntsERI4(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
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

                        transferIntsERI4(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                         eri4_shells_sph, eri4);
                    }
        }

    return eri4;
}

lible::vec4d LIT::calcERI4New(const Structure &structure)
{
    vector<pair<int, int>> l_pairs = returnLPairs(structure.getMaxL());

    vector<ShellPairData> sp_datas = constructShellPairDatas(l_pairs, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(l_pairs, sp_datas);

    size_t dim_ao = structure.getDimAO();
    vec4d eri4(dim_ao, 0);
    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
        for (int lcld = 0; lcld <= lalb; lcld++)
        {
            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);
            int n_sph_c = numSphericals(lc);
            int n_sph_d = numSphericals(ld);
            int n_sph_ab = n_sph_a * n_sph_b;
            int n_sph_cd = n_sph_c * n_sph_d;

            const ShellPairData &sp_data_ab = sp_datas[lalb];
            const ShellPairData &sp_data_cd = sp_datas[lcld];

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

            kernel_eri4_t kernel_eri4 = deployERI4Kernel(la, lb, lc, ld);

            vector<double> eri4_batch(n_sph_ab * n_sph_cd);
            for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
            {
                int bound_cd = (lalb == lcld) ? ipair_ab + 1 : n_pairs_cd;
                for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
                    int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
                    int pos_c = sp_data_cd.coffsets[2 * ipair_cd];
                    int pos_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

                    kernel_eri4(sp_data_ab.cdepths[2 * ipair_ab],
                                sp_data_ab.cdepths[2 * ipair_ab + 1],
                                sp_data_cd.cdepths[2 * ipair_cd],
                                sp_data_cd.cdepths[2 * ipair_cd + 1],
                                &sp_data_ab.exps[pos_a], &sp_data_ab.exps[pos_b],
                                &sp_data_cd.exps[pos_c], &sp_data_cd.exps[pos_d],
                                &sp_data_ab.coords[6 * ipair_ab],
                                &sp_data_ab.coords[6 * ipair_ab + 3],
                                &sp_data_cd.coords[6 * ipair_cd],
                                &sp_data_cd.coords[6 * ipair_cd + 3],
                                &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]],
                                &ecoeffs_cd_tsp[sp_data_cd.offsets_ecoeffs[ipair_cd]],
                                &eri4_batch[0]);

                    transferIntsERI4(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                     eri4_batch, eri4);
                }
            }
        }

    return eri4;
}

void LIT::calcERI4Benchmark(const Structure &structure)
{
    palPrint(std::format("Lible::{:<40}\n", "ERI4 (Shark flat) benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    vector<pair<int, int>> l_pairs = returnLPairs(structure.getMaxL());

    vector<ShellPairData> sp_datas = constructShellPairDatas(l_pairs, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(l_pairs, sp_datas);

    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
        for (int lcld = 0; lcld <= lalb; lcld++)
        {
            auto start{std::chrono::steady_clock::now()};
            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sp_data_cd = sp_datas[lcld];

            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            int lab = la + lb;
            int lcd = lc + ld;

            int dim_a_sph = numSphericals(la);
            int dim_b_sph = numSphericals(lb);
            int dim_c_sph = numSphericals(lc);
            int dim_d_sph = numSphericals(ld);
            int dim_ab_sph = dim_a_sph * dim_b_sph;
            int dim_cd_sph = dim_c_sph * dim_d_sph;

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            int labcd = lab + lcd;
            BoysF boys_f(labcd);

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

            vector<array<int, 3>> idxs_tuv_ab = returnHermiteGaussianIdxs(lab);
            vector<array<int, 3>> idxs_tuv_cd = returnHermiteGaussianIdxs(lcd);

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

                        n_shells_abcd++;
                    }
            else
                for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (int ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
                                   idxs_tuv_ab, idxs_tuv_cd, sp_data_ab, sp_data_cd, boys_f,
                                   eri4_shells_sph, rints, fnx, rints_tmp);

                        n_shells_abcd++;
                    }

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            palPrint(std::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
                                 n_shells_abcd, duration.count()));
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(std::format("done {:.2e} s\n", duration.count()));
}

void LIT::calcERI4BenchmarkNew(const Structure &structure)
{
    palPrint(std::format("Lible::{:<40}\n", "ERI4 (Shark flat new) benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    vector<pair<int, int>> l_pairs = returnLPairs(structure.getMaxL());

    vector<ShellPairData> sp_datas = constructShellPairDatas(l_pairs, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(l_pairs, sp_datas);

    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
        for (int lcld = 0; lcld <= lalb; lcld++)
        {
            auto start{std::chrono::steady_clock::now()};

            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);
            int n_sph_c = numSphericals(lc);
            int n_sph_d = numSphericals(ld);
            int n_sph_ab = n_sph_a * n_sph_b;
            int n_sph_cd = n_sph_c * n_sph_d;

            const ShellPairData &sp_data_ab = sp_datas[lalb];
            const ShellPairData &sp_data_cd = sp_datas[lcld];

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

            kernel_eri4_t kernel_eri4 = deployERI4Kernel(la, lb, lc, ld);

            size_t n_shells_abcd = 0;
            vector<double> eri4_batch(n_sph_ab * n_sph_cd);
            for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
            {
                int bound_cd = (lalb == lcld) ? ipair_ab + 1 : n_pairs_cd;
                for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
                    int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
                    int pos_c = sp_data_cd.coffsets[2 * ipair_cd];
                    int pos_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

                    kernel_eri4(sp_data_ab.cdepths[2 * ipair_ab],
                                sp_data_ab.cdepths[2 * ipair_ab + 1],
                                sp_data_cd.cdepths[2 * ipair_cd],
                                sp_data_cd.cdepths[2 * ipair_cd + 1],
                                &sp_data_ab.exps[pos_a], &sp_data_ab.exps[pos_b],
                                &sp_data_cd.exps[pos_c], &sp_data_cd.exps[pos_d],
                                &sp_data_ab.coords[6 * ipair_ab],
                                &sp_data_ab.coords[6 * ipair_ab + 3],
                                &sp_data_cd.coords[6 * ipair_cd],
                                &sp_data_cd.coords[6 * ipair_cd + 3],
                                &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]],
                                &ecoeffs_cd_tsp[sp_data_cd.offsets_ecoeffs[ipair_cd]],
                                &eri4_batch[0]);

                    n_shells_abcd++;
                }
            }

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            palPrint(std::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
                                 n_shells_abcd, duration.count()));
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(std::format("done {:.2e} s\n", duration.count()));
}

lible::vec2d LIT::calcERI4Diagonal(const Structure &structure)
{
    vector<pair<int, int>> l_pairs = returnLPairs(structure.getMaxL());

    vector<ShellPairData> sp_datas = constructShellPairDatas(l_pairs, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(l_pairs, sp_datas);

    size_t dim_ao = structure.getDimAO();
    vec2d eri4_diagonal(dim_ao, dim_ao, 0);
    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
    {
        const auto &sp_data_ab = sp_datas[lalb];

        auto [la, lb] = l_pairs[lalb];

        int lab = la + lb;

        int dim_a_sph = numSphericals(la);
        int dim_b_sph = numSphericals(lb);
        int dim_ab_sph = dim_a_sph * dim_b_sph;
        int dim_tuv_ab = numHermites(lab);

        int n_pairs_ab = sp_data_ab.n_pairs;

        int labab = 2 * lab;
        BoysF boys_f(labab);

        const vector<double> &ecoeffs_ab = ecoeffs[lalb];
        const vector<double> &ecoeffs_tsp_ab = ecoeffs_tsp[lalb];

        vector<array<int, 3>> idxs_tuv_ab = returnHermiteGaussianIdxs(lab);

        vector<double> rints(dim_tuv_ab * dim_tuv_ab, 0);
        vector<double> fnx(labab + 1, 0);
        vec4d rints_tmp(labab + 1, 0);

        for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
        {
            vector<double> eri4_shells_sph(dim_ab_sph * dim_ab_sph, 0);
            kernelERI4Diagonal(lab, ipair_ab, ecoeffs_ab, ecoeffs_tsp_ab, idxs_tuv_ab, boys_f,
                               sp_data_ab, rints_tmp, eri4_shells_sph, fnx, rints);

            transferIntsERI4Diag(ipair_ab, sp_data_ab, eri4_shells_sph, eri4_diagonal);
        }
    }

    return eri4_diagonal;
}