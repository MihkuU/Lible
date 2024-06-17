#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/util.hpp>
#include <lible/ints/ints_util.hpp>
#include <lible/ints/spherical_trafo.hpp>

#include <fmt/core.h>
#include <libxsmm_source.h>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LIT = lible::ints::two;

using std::pair, std::vector;

namespace lible::ints::two
{
    void kernelERI4SharkFlat_libxsmm(const int lab, const int lcd,
                                     const size_t ipair_ab, const size_t ipair_cd,
                                     const vector<double> &ecoeffs_lalb,
                                     const vector<double> &ecoeffs_lcld,
                                     const vector<MD::IdxsTUV> &idxs_tuv_ab,
                                     const vector<MD::IdxsTUV> &idxs_tuv_cd,
                                     const ShellPairData &shell_pair_data_ab,
                                     const ShellPairData &shell_pair_data_cd,
                                     const BoysF &boys_f, vector<double> &eri4_shells_sph,
                                     const libxsmm_mmfunction<double> &kernel_E_X_R,
                                     const libxsmm_mmfunction<double> &kernel_T_X_E,
                                     vector<double> &rints, vector<double> &fnx,
                                     vec4d &rints_tmp)
    {
        int labcd = lab + lcd;

        const auto &[exps_a, exps_b] = shell_pair_data_ab.exps[ipair_ab];
        const auto &[exps_c, exps_d] = shell_pair_data_cd.exps[ipair_cd];
        const auto &[xyz_a, xyz_b] = shell_pair_data_ab.coords[ipair_ab];
        const auto &[xyz_c, xyz_d] = shell_pair_data_cd.coords[ipair_cd];

        const auto &offsets_ecoeffs_ab = shell_pair_data_ab.offsets_ecoeffs;
        const auto &offsets_ecoeffs_cd = shell_pair_data_cd.offsets_ecoeffs;

        size_t offset_prim_ab = shell_pair_data_ab.offsets_prims[ipair_ab];
        size_t offset_prim_cd = shell_pair_data_cd.offsets_prims[ipair_cd];

        size_t ka = exps_a.size();
        size_t kb = exps_b.size();
        size_t kc = exps_c.size();
        size_t kd = exps_d.size();

        arma::vec::fixed<3> A{xyz_a[0], xyz_a[1], xyz_a[2]};
        arma::vec::fixed<3> B{xyz_b[0], xyz_b[1], xyz_b[2]};
        arma::vec::fixed<3> C{xyz_c[0], xyz_c[1], xyz_c[2]};
        arma::vec::fixed<3> D{xyz_d[0], xyz_d[1], xyz_d[2]};

        int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
        int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

        int la = shell_pair_data_ab.la, lb = shell_pair_data_ab.lb;
        int lc = shell_pair_data_cd.la, ld = shell_pair_data_cd.lb;
        int dim_sph_a = dimSphericals(la);
        int dim_sph_b = dimSphericals(lb);
        int dim_sph_c = dimSphericals(lc);
        int dim_sph_d = dimSphericals(ld);
        int dim_sph_ab = dim_sph_a * dim_sph_b;
        int dim_sph_cd = dim_sph_c * dim_sph_d;
        int dim_E_x_R = dim_tuv_ab * dim_sph_cd;
        vector<double> X(ka * kb * dim_E_x_R, 0);

        // libxsmm_mmfunction<double> kernel_E_X_R_test(LIBXSMM_GEMM_FLAG_NONE, dim_tuv_cd,
        //                                              dim_sph_ab, dim_tuv_ab, 1.0, 1.0);

        // libxsmm_mmfunction<double> kernel_T_X_E_test(LIBXSMM_GEMM_FLAG_NONE, dim_sph_ab,
        //                                              dim_sph_cd, dim_tuv_cd, 1.0, 1.0);

        size_t iab = 0;
        for (size_t ia = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++)
            {
                size_t pos_X = iab * dim_E_x_R;

                size_t icd = 0;
                for (size_t ic = 0; ic < kc; ic++)
                    for (size_t id = 0; id < kd; id++)
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

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                        MD::calcRInts(lab, lcd, alpha, fac, RPQ, fnx, idxs_tuv_ab, idxs_tuv_cd,
                                      rints_tmp, rints);                                                

                        size_t iprim_cd = offset_prim_cd + icd;
                        size_t pos_ecoeffs_cd = offsets_ecoeffs_cd[iprim_cd];
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim_tuv_ab,
                                    dim_sph_cd, dim_tuv_cd, 1.0, &rints[0], dim_tuv_cd,
                                    &ecoeffs_lcld[pos_ecoeffs_cd], dim_tuv_cd, 1.0, &X[pos_X],
                                    dim_sph_cd);

                        // kernel_E_X_R(&ecoeffs_lalb[pos_ecoeffs_ab], &rints[0], &X[pos_X]);
                        // kernel_E_X_R_test(&rints[0], &ecoeffs_lalb[pos_ecoeffs_ab], &X[pos_X]);
                        icd++;
                    }
                iab++;
            }

        for (size_t ia = 0, iab = 0; ia < ka; ia++)
            for (size_t ib = 0; ib < kb; ib++, iab++)
            {            
                size_t iprim_ab = offset_prim_ab + iab;
                size_t pos_ecoeffs_ab = offsets_ecoeffs_ab[iprim_ab];
                size_t pos_X = iab * dim_E_x_R;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_ab, dim_sph_cd,
                            dim_tuv_ab, 1.0, &ecoeffs_lalb[pos_ecoeffs_ab], dim_tuv_ab,
                            &X[pos_X], dim_sph_cd, 1.0, &eri4_shells_sph[0], dim_sph_cd);

                // kernel_T_X_E(&X[pos_X], &ecoeffs_lcld[pos_ecoeffs_cd], &eri4_shells_sph[0]);
                // kernel_T_X_E_test(&X[pos_X], &ecoeffs_lcld[pos_ecoeffs_cd], &eri4_shells_sph[0]);
            }
    }
}

lible::vec4d LIT::calcERI4SharkFlat(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}", "SHARK ERI4 flat..."));

    auto start{std::chrono::steady_clock::now()};

    libxsmm_init();

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> shell_pair_datas(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        shell_pair_datas[ipair] = ShellPairData(la, lb, structure);
    }

    vector<vector<double>> ecoeffs(l_pairs.size()), ecoeffs_t(l_pairs.size());    
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        int lab = la + lb;
        int n_ecoeffs_sph = (2 * la + 1) * (2 * lb + 1) *
                            ((lab + 1) * (lab + 2) * (lab + 3) / 6) *
                            shell_pair_datas[ipair].n_prim_pairs;

        vector<double> ecoeffs_ipair(n_ecoeffs_sph), ecoeffs_t_ipair(n_ecoeffs_sph);
        MD::calcECoeffsSpherical(la, lb, shell_pair_datas[ipair], ecoeffs_ipair, ecoeffs_t_ipair);
        ecoeffs[ipair] = std::move(ecoeffs_ipair);
        ecoeffs_t[ipair] = std::move(ecoeffs_t_ipair);
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
            const vector<double> &ecoeffs_cd = ecoeffs[lcld];
            // const vector<double> &ecoeffs_cd = ecoeffs_t[lcld];

            const auto &shell_pair_data_ab = shell_pair_datas[lalb];
            const auto &shell_pair_data_cd = shell_pair_datas[lcld];

            size_t n_pairs_ab = shell_pair_data_ab.n_pairs;
            size_t n_pairs_cd = shell_pair_data_cd.n_pairs;

            vector<MD::IdxsTUV> idxs_tuv_ab = MD::returnIdxsTUV(lab);
            vector<MD::IdxsTUV> idxs_tuv_cd = MD::returnIdxsTUV(lcd);
            
            int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
            int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

            libxsmm_mmfunction<double> kernel_E_X_R(LIBXSMM_GEMM_FLAG_NONE, dim_ab_sph,
                                                    dim_tuv_cd, dim_tuv_ab, 1.0, 1.0);
            libxsmm_mmfunction<double> kernel_T_X_E(LIBXSMM_GEMM_FLAG_NONE, dim_cd_sph,
                                                    dim_ab_sph, dim_tuv_cd, 1.0, 1.0);

            vector<double> rints(dim_tuv_ab * dim_tuv_cd, 0);
            vector<double> fnx(labcd + 1, 0);
            vec4d rints_tmp(labcd + 1, 0);

            if (lalb == lcld)
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4SharkFlat_libxsmm(lab, lcd, ipair_ab, ipair_cd,
                                                    ecoeffs_ab, ecoeffs_cd,
                                                    idxs_tuv_ab, idxs_tuv_cd,
                                                    shell_pair_data_ab, shell_pair_data_cd,
                                                    boys_f, eri4_shells_sph,
                                                    kernel_E_X_R, kernel_T_X_E,
                                                    rints, fnx, rints_tmp);

                        // kernelERI4SharkFlat(lab, lcd, ipair_ab, ipair_cd,
                        //                     ecoeffs_ab, ecoeffs_cd,
                        //                     idxs_tuv_ab, idxs_tuv_cd,
                        //                     shell_pair_data_ab, shell_pair_data_cd,
                        //                     boys_f, eri4_shells_sph, rints,
                        //                     fnx, rints_tmp);

                        transferIntegrals(ipair_ab, ipair_cd, shell_pair_data_ab,
                                          shell_pair_data_cd, eri4_shells_sph, eri4);
                    }
            else
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4SharkFlat_libxsmm(lab, lcd, ipair_ab, ipair_cd,
                                                    ecoeffs_ab, ecoeffs_cd,
                                                    idxs_tuv_ab, idxs_tuv_cd,
                                                    shell_pair_data_ab, shell_pair_data_cd,
                                                    boys_f, eri4_shells_sph,
                                                    kernel_E_X_R, kernel_T_X_E,
                                                    rints, fnx, rints_tmp);

                        // kernelERI4SharkFlat(lab, lcd, ipair_ab, ipair_cd,
                        //                     ecoeffs_ab, ecoeffs_cd,
                        //                     idxs_tuv_ab, idxs_tuv_cd,
                        //                     shell_pair_data_ab, shell_pair_data_cd,
                        //                     boys_f, eri4_shells_sph, rints,
                        //                     fnx, rints_tmp);

                        transferIntegrals(ipair_ab, ipair_cd, shell_pair_data_ab,
                                          shell_pair_data_cd, eri4_shells_sph, eri4);
                    }
        }

    libxsmm_finalize();

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));    

    return eri4;
}

void LIT::calcERI4BenchmarkSharkFlat(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}\n", "ERI4 (Shark flat) benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> shell_pair_datas(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        shell_pair_datas[ipair] = ShellPairData(la, lb, structure);
    }

    vector<vector<double>> ecoeffs(l_pairs.size()), ecoeffs_t(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        int lab = la + lb;
        int n_ecoeffs_sph = (2 * la + 1) * (2 * lb + 1) *
                            ((lab + 1) * (lab + 2) * (lab + 3) / 6) *
                            shell_pair_datas[ipair].n_prim_pairs;

        vector<double> ecoeffs_ipair(n_ecoeffs_sph), ecoeffs_t_ipair(n_ecoeffs_sph);
        MD::calcECoeffsSpherical(la, lb, shell_pair_datas[ipair], ecoeffs_ipair, ecoeffs_t_ipair);
        ecoeffs[ipair] = std::move(ecoeffs_ipair);
        ecoeffs_t[ipair] = std::move(ecoeffs_t_ipair);        
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

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            // const vector<double> &ecoeffs_cd = ecoeffs[lcld];
            const vector<double> &ecoeffs_cd = ecoeffs_t[lcld];

            const auto &shell_pair_data_ab = shell_pair_datas[lalb];
            const auto &shell_pair_data_cd = shell_pair_datas[lcld];

            size_t n_pairs_ab = shell_pair_data_ab.n_pairs;
            size_t n_pairs_cd = shell_pair_data_cd.n_pairs;

            vector<MD::IdxsTUV> idxs_tuv_ab = MD::returnIdxsTUV(lab);
            vector<MD::IdxsTUV> idxs_tuv_cd = MD::returnIdxsTUV(lcd);

            int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
            int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

            libxsmm_mmfunction<double> kernel_E_X_R(LIBXSMM_GEMM_FLAG_NONE, dim_ab_sph,
                                                    dim_tuv_cd, dim_tuv_ab, 1.0, 1.0);
            libxsmm_mmfunction<double> kernel_T_X_E(LIBXSMM_GEMM_FLAG_NONE, dim_ab_sph,
                                                    dim_cd_sph, dim_tuv_cd, 1.0, 1.0);

            vector<double> rints(dim_tuv_ab * dim_tuv_cd, 0);
            vector<double> fnx(labcd + 1, 0);
            vec4d rints_tmp(labcd + 1, 0);

            size_t n_shells_abcd = 0;
            if (lalb == lcld)
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd <= ipair_ab; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4SharkFlat_libxsmm(lab, lcd, ipair_ab, ipair_cd,
                                                    ecoeffs_ab, ecoeffs_cd,
                                                    idxs_tuv_ab, idxs_tuv_cd,
                                                    shell_pair_data_ab, shell_pair_data_cd,
                                                    boys_f, eri4_shells_sph,
                                                    kernel_E_X_R, kernel_T_X_E,
                                                    rints, fnx, rints_tmp);

                        // kernelERI4SharkFlat(lab, lcd, ipair_ab, ipair_cd,
                        //                     ecoeffs_ab, ecoeffs_cd,
                        //                     idxs_tuv_ab, idxs_tuv_cd,
                        //                     shell_pair_data_ab, shell_pair_data_cd,
                        //                     boys_f, eri4_shells_sph, rints,
                        //                     fnx, rints_tmp);

                        n_shells_abcd += 4;
                    }
            else
                for (size_t ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (size_t ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4SharkFlat_libxsmm(lab, lcd, ipair_ab, ipair_cd,
                                                    ecoeffs_ab, ecoeffs_cd,
                                                    idxs_tuv_ab, idxs_tuv_cd,
                                                    shell_pair_data_ab, shell_pair_data_cd,
                                                    boys_f, eri4_shells_sph,
                                                    kernel_E_X_R, kernel_T_X_E,
                                                    rints, fnx, rints_tmp);

                        // kernelERI4SharkFlat(lab, lcd, ipair_ab, ipair_cd,
                        //                     ecoeffs_ab, ecoeffs_cd,
                        //                     idxs_tuv_ab, idxs_tuv_cd,
                        //                     shell_pair_data_ab, shell_pair_data_cd,
                        //                     boys_f, eri4_shells_sph, rints,
                        //                     fnx, rints_tmp);

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

void LIT::kernelERI4SharkFlat(const int lab, const int lcd,
                              const size_t ipair_ab, const size_t ipair_cd,
                              const vector<double> &ecoeffs_lalb,
                              const vector<double> &ecoeffs_lcld,                            
                              const vector<MD::IdxsTUV> &idxs_tuv_ab,
                              const vector<MD::IdxsTUV> &idxs_tuv_cd,
                              const ShellPairData &shell_pair_data_ab,
                              const ShellPairData &shell_pair_data_cd,
                              const BoysF &boys_f, vector<double> &eri4_shells_sph,
                              vector<double> &rints, vector<double> &fnx,
                              vec4d &rints_tmp)
{
    int labcd = lab + lcd;

    const auto &[exps_a, exps_b] = shell_pair_data_ab.exps[ipair_ab];
    const auto &[exps_c, exps_d] = shell_pair_data_cd.exps[ipair_cd];
    const auto &[xyz_a, xyz_b] = shell_pair_data_ab.coords[ipair_ab];
    const auto &[xyz_c, xyz_d] = shell_pair_data_cd.coords[ipair_cd];

    const auto &offsets_ecoeffs_ab = shell_pair_data_ab.offsets_ecoeffs;
    const auto &offsets_ecoeffs_cd = shell_pair_data_cd.offsets_ecoeffs;

    size_t offset_prim_ab = shell_pair_data_ab.offsets_prims[ipair_ab];
    size_t offset_prim_cd = shell_pair_data_cd.offsets_prims[ipair_cd];

    size_t ka = exps_a.size();
    size_t kb = exps_b.size();
    size_t kc = exps_c.size();
    size_t kd = exps_d.size();

    arma::vec::fixed<3> A{xyz_a[0], xyz_a[1], xyz_a[2]};
    arma::vec::fixed<3> B{xyz_b[0], xyz_b[1], xyz_b[2]};
    arma::vec::fixed<3> C{xyz_c[0], xyz_c[1], xyz_c[2]};
    arma::vec::fixed<3> D{xyz_d[0], xyz_d[1], xyz_d[2]};

    int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
    int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;    

    int la = shell_pair_data_ab.la, lb = shell_pair_data_ab.lb;
    int lc = shell_pair_data_cd.la, ld = shell_pair_data_cd.lb;
    int dim_sph_a = dimSphericals(la);
    int dim_sph_b = dimSphericals(lb);
    int dim_sph_c = dimSphericals(lc);
    int dim_sph_d = dimSphericals(ld);
    int dim_sph_ab = dim_sph_a * dim_sph_b;
    int dim_sph_cd = dim_sph_c * dim_sph_d;
    int dim_E_x_R = dim_sph_ab * dim_tuv_cd;
    vector<double> X(kc * kd * dim_E_x_R, 0);

    size_t iab = 0;
    for (size_t ia = 0, iab = 0; ia < ka; ia++)
        for (size_t ib = 0; ib < kb; ib++)
        {
            size_t iprim_ab = offset_prim_ab + iab;
            size_t pos_ecoeffs_ab = offsets_ecoeffs_ab[iprim_ab];

            size_t icd = 0;
            for (size_t ic = 0, icd = 0; ic < kc; ic++)
                for (size_t id = 0; id < kd; id++)
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

                    MD::calcRInts(lab, lcd, alpha, RPQ, fnx, idxs_tuv_ab, idxs_tuv_cd,
                                  rints_tmp, rints);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));

                    size_t pos_X = icd * dim_E_x_R;
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_ab, dim_tuv_cd,
                                dim_tuv_ab, fac, &ecoeffs_lalb[pos_ecoeffs_ab], dim_tuv_ab,
                                &rints[0], dim_tuv_cd, 1.0, &X[pos_X], dim_tuv_cd);
                    icd++;
                }
            iab++;
        }

    for (size_t ic = 0, icd = 0; ic < kc; ic++)
        for (size_t id = 0; id < kd; id++, icd++)
        {
            size_t iprim_cd = offset_prim_cd + icd;
            size_t pos_ecoeffs_cd = offsets_ecoeffs_cd[iprim_cd];
            size_t pos_X = icd * dim_E_x_R;            

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim_sph_ab,
                        dim_sph_cd, dim_tuv_cd, 1.0, &X[pos_X], dim_tuv_cd,
                        &ecoeffs_lcld[pos_ecoeffs_cd], dim_tuv_cd, 1.0,
                        &eri4_shells_sph[0], dim_sph_cd);
        }
}