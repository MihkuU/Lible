#include <lible/utils.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <chrono>
#include <cstring>
#include <format>
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
                    const BoysF &boys_f, vector<double> &eri4_shells_sph, vector<double> &fnx)
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

                        vector<double> rints = calcRIntsMatrix(labcd, fac, alpha, xyz_pq.data(),
                                                               fnx.data(), idxs_tuv_ab,
                                                               idxs_tuv_cd);

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

        int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
        int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
        int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];
        int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
        for (int mu = 0, idx = 0; mu < dim_a_sph; mu++)
            for (int nu = 0; nu < dim_b_sph; nu++)
                for (int ka = 0; ka < dim_c_sph; ka++)
                    for (int ta = 0; ta < dim_d_sph; ta++, idx++)
                    {
                        double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                        double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                        double norm_c = sp_data_cd.norms[ofs_norm_c + ka];
                        double norm_d = sp_data_cd.norms[ofs_norm_d + ta];

                        eri4_shells_sph[idx] *= norm_a * norm_b * norm_c * norm_d;
                    }
    }

    // vec4d kernelERI4Test(const int ipair_ab, const int ipair_cd, const vector<double> &ecoeffs_bra,
    //                      const vector<double> &ecoeffs_ket, const ShellPairData &sp_data_ab,
    //                      const ShellPairData &sp_data_cd, const BoysGrid &boys_grid)
    vec4d kernelERI4Test(const int ipair_ab, const int ipair_cd, const vector<double> &ecoeffs_bra,
                         const vector<double> &ecoeffs_ket, const ShellPairData &sp_data_ab,
                         const ShellPairData &sp_data_cd, const BoysGrid &boys_grid)
    {
        const int la = sp_data_ab.la;
        const int lb = sp_data_ab.lb;
        const int lc = sp_data_cd.la;
        const int ld = sp_data_cd.lb;
        const int lab = la + lb;
        const int lcd = lc + ld;
        const int labcd = lab + lcd;

        const int n_sph_a = numSphericals(la);
        const int n_sph_b = numSphericals(lb);
        const int n_sph_c = numSphericals(lc);
        const int n_sph_d = numSphericals(ld);
        const int n_sph_ab = n_sph_a * n_sph_b;
        const int n_sph_cd = n_sph_c * n_sph_d;
        const int n_hermite_ab = numHermites(lab);
        const int n_hermite_cd = numHermites(lcd);

        const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
        const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
        const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
        const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
        const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
        const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
        const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
        const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

        const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab];
        const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
        const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd];
        const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];

        vector<array<int, 3>> idxs_tuv_ab = getHermiteGaussianIdxs(lab);
        vector<array<int, 3>> idxs_tuv_cd = getHermiteGaussianIdxs(lcd);

        int n_r_rows = (cdepth_a * cdepth_b * n_hermite_ab);
        int n_r_cols = (cdepth_c * cdepth_d * n_hermite_cd);
        int n_r_ints = n_r_rows * n_r_cols;
        vector<double> rints(n_r_ints, 0);
        for (int ia = 0, iab = 0, iabcd = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
                for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                    for (int id = 0; id < cdepth_d; id++, icd++)
                    {
                        double a = sp_data_ab.exps[cofs_a + ia];
                        double b = sp_data_ab.exps[cofs_b + ib];
                        double c = sp_data_cd.exps[cofs_c + ic];
                        double d = sp_data_cd.exps[cofs_d + id];

                        double p = a + b;
                        double q = c + d;
                        double alpha = p * q / (p + q);

                        array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                               (a * xyz_a[1] + b * xyz_b[1]) / p,
                                               (a * xyz_a[2] + b * xyz_b[2]) / p};

                        array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                               (c * xyz_c[1] + d * xyz_d[1]) / q,
                                               (c * xyz_c[2] + d * xyz_d[2]) / q};

                        array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                                xyz_p[1] - xyz_q[1],
                                                xyz_p[2] - xyz_q[2]};

                        double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                        double x = alpha * (xx * xx + xy * xy + xz * xz);
                        vector<double> fnx = calcBoysF(labcd, x, boys_grid);

                        // std::fill(fnx.begin(), fnx.end(), 1); // tmp

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                        int ofs_row = iab * n_hermite_ab;
                        int ofs_col = icd * n_hermite_cd;

                        calcRIntsMatrixTest(labcd, n_r_cols, ofs_row, ofs_col, fac, alpha,
                                            xyz_pq.data(), fnx.data(), idxs_tuv_ab, idxs_tuv_cd,
                                            &rints[0]);

                        iabcd++;
                    }

        int m = cdepth_a * cdepth_b * n_hermite_ab;
        int n = n_sph_cd;
        int k = cdepth_c * cdepth_d * n_hermite_cd;
        int n_R_x_E = cdepth_a * cdepth_b * (n_hermite_ab * n_sph_cd);
        int ofs_E_bra = sp_data_ab.offsets_ecoeffs[ipair_ab];
        int ofs_E_ket = sp_data_cd.offsets_ecoeffs[ipair_cd];

        vector<double> R_x_E(n_R_x_E, 0);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0, &rints[0], k,
                    &ecoeffs_ket[ofs_E_ket], k, 1.0, &R_x_E[0], n);

        m = n_sph_ab;
        n = n_sph_cd;
        k = cdepth_a * cdepth_b * n_hermite_ab;

        vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs_bra[ofs_E_bra],
                    k, &R_x_E[0], n, 1.0, &eri4_batch[0], n);

        return eri4_batch;
    }

    void kernelERI4Diagonal(const int lab, const int ipair_ab, const vector<double> &ecoeffs_ab,
                            const vector<double> &ecoeffs_ab_tsp,
                            const vector<array<int, 3>> &idxs_tuv,
                            const BoysF &boys_f, const ShellPairData &sp_data_ab,
                            vector<double> &eri4_shells_sph, vector<double> &fnx)
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
                        boys_f.calcFnx(labab, x, fnx); // TODO: change this shit to *double fnx

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));

                        vector<double> rints = calcRIntsMatrix(labab, fac, alpha, xyz_pq.data(),
                                                               fnx.data(), idxs_tuv, idxs_tuv);

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

    void transferIntsERI4Diag(const int ipair_ab, const ShellPairData &sp_data_ab,
                              const vector<double> &eri4_shells_sph, vec2d &eri4_diagonal)
    {
        int dim_a = numSphericals(sp_data_ab.la);
        int dim_b = numSphericals(sp_data_ab.lb);
        int dim_ab = dim_a * dim_b;

        int pos_a = sp_data_ab.offsets_sph[2 * ipair_ab];
        int pos_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
        int pos_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
        int pos_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];

        for (int mu = 0, munu = 0; mu < dim_a; mu++)
            for (int nu = 0; nu < dim_b; nu++, munu++)
            {
                int munumunu = munu * dim_ab + munu;

                double norm_a = sp_data_ab.norms[pos_norm_a + mu];
                double norm_b = sp_data_ab.norms[pos_norm_b + nu];

                double normalized_int = eri4_shells_sph[munumunu] * norm_a * norm_b * norm_a * norm_b;

                int a = pos_a + mu;
                int b = pos_b + nu;

                eri4_diagonal(a, b) = normalized_int;
                eri4_diagonal(b, a) = normalized_int;
            }
    }
}

lible::vec4d LIT::calcERI4(const Structure &structure)
{
    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

    // auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(sp_datas);

    size_t dim_ao = structure.getDimAO();
    vec4d eri4(Fill(0), dim_ao);
    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
        for (int lcld = 0; lcld <= lalb; lcld++)
        {
            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);
            int n_sph_c = numSphericals(lc);
            int n_sph_d = numSphericals(ld);

            const ShellPairData &sp_data_ab = sp_datas[lalb];
            const ShellPairData &sp_data_cd = sp_datas[lcld];

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            // const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            // const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

            // kernel_eri4_t kernel_eri4 = deployERI4Kernel(la, lb, lc, ld);
            ERI4Kernel eri4_kernel = deployERI4Kernel(sp_data_ab, sp_data_cd);

            for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
            {
                int bound_cd = (lalb == lcld) ? ipair_ab + 1 : n_pairs_cd;
                for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd);

                    int ofs_a = sp_data_ab.offsets_sph[2 * ipair_ab];
                    int ofs_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
                    int ofs_c = sp_data_cd.offsets_sph[2 * ipair_cd];
                    int ofs_d = sp_data_cd.offsets_sph[2 * ipair_cd + 1];

                    for (int ia = 0; ia < n_sph_a; ia++)
                        for (int ib = 0; ib < n_sph_b; ib++)
                            for (int ic = 0; ic < n_sph_c; ic++)
                                for (int id = 0; id < n_sph_d; id++)
                                {
                                    int mu = ofs_a + ia;
                                    int nu = ofs_b + ib;
                                    int ka = ofs_c + ic;
                                    int ta = ofs_d + id;

                                    double integral = eri4_batch(ia, ib, ic, id);
                                    eri4(mu, nu, ka, ta) = integral;
                                    eri4(mu, nu, ta, ka) = integral;
                                    eri4(nu, mu, ka, ta) = integral;
                                    eri4(nu, mu, ta, ka) = integral;
                                    eri4(ka, ta, mu, nu) = integral;
                                    eri4(ka, ta, nu, mu) = integral;
                                    eri4(ta, ka, mu, nu) = integral;
                                    eri4(ta, ka, nu, mu) = integral;
                                }
                }
            }
        }

    return eri4;
}

void LIT::calcERI4Benchmark(const Structure &structure)
{
    palPrint(std::format("Lible::{:<40}\n", "ERI4 (Shark flat) benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(sp_datas);

    double sum_eri4 = 0;

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
            // vector<double> ecoeffs_ab = ecoeffs[lalb];                  // tmp
            // vector<double> ecoeffs_cd_tsp = ecoeffs_tsp[lcld];          // tmp
            // std::fill(ecoeffs_ab.begin(), ecoeffs_ab.end(), 1);         // tmp
            // std::fill(ecoeffs_cd_tsp.begin(), ecoeffs_cd_tsp.end(), 1); // tmp

            vector<array<int, 3>> idxs_tuv_ab = getHermiteGaussianIdxs(lab);
            vector<array<int, 3>> idxs_tuv_cd = getHermiteGaussianIdxs(lcd);

            int dim_tuv_ab = (lab + 1) * (lab + 2) * (lab + 3) / 6;
            int dim_tuv_cd = (lcd + 1) * (lcd + 2) * (lcd + 3) / 6;

            vector<double> rints(dim_tuv_ab * dim_tuv_cd, 0);
            vector<double> fnx(labcd + 1, 0);
            vec4d rints_tmp(Fill(0), labcd + 1);

            size_t n_shells_abcd = 0;
            for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
            {
                int bound_cd = (lalb == lcld) ? ipair_ab + 1 : n_pairs_cd;
                for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                    kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
                               idxs_tuv_ab, idxs_tuv_cd, sp_data_ab, sp_data_cd, boys_f,
                               eri4_shells_sph, fnx);

                    for (double x : eri4_shells_sph)
                        sum_eri4 += std::fabs(x);

                    n_shells_abcd++;
                }
            }

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            palPrint(std::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
                                 n_shells_abcd, duration.count()));
        }

    palPrint(std::format("   sum_eri4 = {:16.12f}\n", sum_eri4));

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(std::format("done {:.2e} s\n", duration.count()));
}

void LIT::calcERI4BenchmarkTest(const Structure &structure)
{
    palPrint(std::format("Lible::{:<40}\n", "ERI4 (Shark test) benchmark..."));

    auto start{std::chrono::steady_clock::now()};

    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

    double sum_eri4 = 0;

    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
        for (int lcld = 0; lcld <= lalb; lcld++)
        {
            auto start{std::chrono::steady_clock::now()};

            auto [la, lb] = l_pairs[lalb];
            auto [lc, ld] = l_pairs[lcld];

            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sp_data_cd = sp_datas[lcld];

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_pairs_cd = sp_data_cd.n_pairs;

            // int labcd = la + lb + lc + ld;
            // BoysGrid boys_grid(labcd);

            // vector<double> ecoeffs_bra = ecoeffsSHARKSPData(sp_data_ab);
            // vector<double> ecoeffs_ket = ecoeffsSHARKSPData(sp_data_cd);

            ERI4Kernel eri4_kernel = deployERI4Kernel(sp_data_ab, sp_data_cd);

            size_t n_shells_abcd = 0;
            for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
            {
                int bound_cd = (lalb == lcld) ? ipair_ab + 1 : n_pairs_cd;
                for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
                {
                    // vec4d eri4_batch = kernelERI4Test(ipair_ab, ipair_cd, ecoeffs_bra, ecoeffs_ket,
                    //                                   sp_data_ab, sp_data_cd, boys_grid);

                    vec4d eri4_batch = eri4_kernel(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd);

                    for (double x : eri4_batch)
                        sum_eri4 += std::fabs(x);

                    n_shells_abcd++;
                }
            }

            auto end{std::chrono::steady_clock::now()};
            std::chrono::duration<double> duration{end - start};

            palPrint(std::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
                                 n_shells_abcd, duration.count()));
        }

    palPrint(std::format("   sum_eri4 = {:16.12f}\n", sum_eri4));

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(std::format("done {:.2e} s\n", duration.count()));
}

// void LIT::calcERI4BenchmarkNew(const Structure &structure)
// {
//     palPrint(std::format("Lible::{:<40}\n", "ERI4 (Shark flat new) benchmark..."));

//     auto start{std::chrono::steady_clock::now()};

//     vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

//     vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

//     // auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(sp_datas);

//     for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
//         for (int lcld = 0; lcld <= lalb; lcld++)
//         {
//             auto start{std::chrono::steady_clock::now()};

//             auto [la, lb] = l_pairs[lalb];
//             auto [lc, ld] = l_pairs[lcld];

//             const ShellPairData &sp_data_ab = sp_datas[lalb];
//             const ShellPairData &sp_data_cd = sp_datas[lcld];

//             int n_pairs_ab = sp_data_ab.n_pairs;
//             int n_pairs_cd = sp_data_cd.n_pairs;

//             // const vector<double> &ecoeffs_ab = ecoeffs[lalb];
//             // const vector<double> &ecoeffs_cd_tsp = ecoeffs_tsp[lcld];

//             // kernel_eri4_t kernel_eri4 = deployERI4Kernel(la, lb, lc, ld);
//             ERI4Kernel eri4_kernel = deployERI4Kernel(sp_data_ab, sp_data_cd);

//             size_t n_shells_abcd = 0;
//             for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
//             {
//                 int bound_cd = (lalb == lcld) ? ipair_ab + 1 : n_pairs_cd;
//                 for (int ipair_cd = 0; ipair_cd < bound_cd; ipair_cd++)
//                 {
//                     vec4d eri4_batch = kernel_eri4(ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
//                                                    sp_data_ab, sp_data_cd);

//                     n_shells_abcd++;
//                 }
//             }

//             auto end{std::chrono::steady_clock::now()};
//             std::chrono::duration<double> duration{end - start};

//             palPrint(std::format("   {} {} {} {} ; {:10} ; {:.2e} s\n", la, lb, lc, ld,
//                                  n_shells_abcd, duration.count()));
//         }

//     auto end{std::chrono::steady_clock::now()};
//     std::chrono::duration<double> duration{end - start};
//     palPrint(std::format("done {:.2e} s\n", duration.count()));
// }

lible::vec2d LIT::calcERI4Diagonal(const Structure &structure)
{
    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalSPDatas_BraKet(sp_datas);

    size_t dim_ao = structure.getDimAO();
    vec2d eri4_diagonal(Fill(0), dim_ao, dim_ao);
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

        vector<array<int, 3>> idxs_tuv_ab = getHermiteGaussianIdxs(lab);

        vector<double> rints(dim_tuv_ab * dim_tuv_ab, 0);
        vector<double> fnx(labab + 1, 0);
        vec4d rints_tmp(Fill(0), labab + 1);

        for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
        {
            vector<double> eri4_shells_sph(dim_ab_sph * dim_ab_sph, 0);
            kernelERI4Diagonal(lab, ipair_ab, ecoeffs_ab, ecoeffs_tsp_ab, idxs_tuv_ab, boys_f,
                               sp_data_ab, eri4_shells_sph, fnx);

            transferIntsERI4Diag(ipair_ab, sp_data_ab, eri4_shells_sph, eri4_diagonal); // TODO remove this function, do it here
        }
    }

    return eri4_diagonal;
}

array<lible::vec4d, 12> LIT::kernelERI4Deriv1(const int ipair_ab, const int ipair_cd,
                                              const vector<double> &ecoeffs_ab,
                                              const vector<double> &ecoeffs1_ab,
                                              const vector<double> &ecoeffs_cd_tsp,
                                              const vector<double> &ecoeffs1_cd_tsp,
                                              const BoysGrid &boys_grid,
                                              const ShellPairData &sp_data_ab,
                                              const ShellPairData &sp_data_cd)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sp_data_cd.la;
    int ld = sp_data_cd.lb;
    int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
    int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
    int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
    int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
    int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
    int eofs0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    int eofs1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    int eofs0_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];
    int eofs1_cd = sp_data_cd.offsets_ecoeffs_deriv1[ipair_cd];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
    const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
    const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
    const double *pecoeffs_ab = &ecoeffs_ab[eofs0_ab];
    const double *pecoeffs1_ab = &ecoeffs1_ab[eofs1_ab];
    const double *pecoeffs_cd_tsp = &ecoeffs_cd_tsp[eofs0_cd];
    const double *pecoeffs_deriv1_cd_tsp = &ecoeffs1_cd_tsp[eofs1_cd];

    int lab = la + lb;
    int lcd = lc + ld;
    int labcd = lab + lcd;

    vector<array<int, 3>> idxs_tuv_ab = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> idxs_tuv_cd = getHermiteGaussianIdxs(lcd);

    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_c = numSphericals(lc);
    int n_sph_d = numSphericals(ld);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_sph_cd = n_sph_c * n_sph_d;
    int n_sph_abcd = n_sph_ab * n_sph_cd;
    int n_hermite_ab = numHermites(lab);
    int n_hermite_cd = numHermites(lcd);
    int n_hermite_abcd = n_hermite_ab * n_hermite_cd;
    int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;

    int n_R_x_E = n_hermite_ab * n_sph_cd;

    vector<double> R_x_E(12 * cdepth_a * cdepth_b * n_R_x_E, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            int ofs_R_x_E = 12 * iab * n_R_x_E;
            vector<double> R_x_E_P_R_(6 * n_R_x_E, 0);
            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double c = exps_c[ic];
                    double d = exps_d[id];

                    double p = a + b;
                    double q = c + d;
                    double alpha = p * q / (p + q);

                    array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                           (a * xyz_a[1] + b * xyz_b[1]) / p,
                                           (a * xyz_a[2] + b * xyz_b[2]) / p};

                    array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                           (c * xyz_c[1] + d * xyz_d[1]) / q,
                                           (c * xyz_c[2] + d * xyz_d[2]) / q};

                    array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                            xyz_p[1] - xyz_q[1],
                                            xyz_p[2] - xyz_q[2]};

                    double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                    double x = alpha * (xx * xx + xy * xy + xz * xz);
                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));

                    vector<double> fns = calcBoysF(labcd + 1, x, boys_grid);

                    vector<double> rints = calcRInts_ERI4_Deriv1(labcd, fac, alpha, xyz_pq.data(),
                                                                 fns.data(), idxs_tuv_ab,
                                                                 idxs_tuv_cd);

                    int ofs_E0_cd = icd * n_ecoeffs_cd;

                    // P
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[0 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 0 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[1 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 1 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[2 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 2 * n_R_x_E], n_sph_cd);

                    // R
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 3 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 4 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 5 * n_R_x_E], n_sph_cd);

                    std::fill(R_x_E_P_R_.begin(), R_x_E_P_R_.end(), 0);

                    // P'
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[4 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[0 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[5 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[1 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[6 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_cd_tsp[ofs_E0_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[2 * n_R_x_E], n_sph_cd);

                    // R'
                    int ofs_E1_cd = 3 * icd * n_ecoeffs_cd;
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_deriv1_cd_tsp[ofs_E1_cd + 0 * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[3 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_deriv1_cd_tsp[ofs_E1_cd + 1 * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[4 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &pecoeffs_deriv1_cd_tsp[ofs_E1_cd + 2 * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[5 * n_R_x_E], n_sph_cd);

                    // P'R' -> CD
                    for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                        for (int ka = 0, kata = 0; ka < n_sph_c; ka++)
                            for (int ta = 0; ta < n_sph_d; ta++, kata++)
                            {
                                int idx = tuv * n_sph_cd + kata;

                                int idx0 = 0 * n_R_x_E + idx;
                                int idx1 = 1 * n_R_x_E + idx;
                                int idx2 = 2 * n_R_x_E + idx;

                                int idx3 = 3 * n_R_x_E + idx;
                                int idx4 = 4 * n_R_x_E + idx;
                                int idx5 = 5 * n_R_x_E + idx;

                                int idx0_ = ofs_R_x_E + 6 * n_R_x_E + idx;
                                int idx1_ = ofs_R_x_E + 7 * n_R_x_E + idx;
                                int idx2_ = ofs_R_x_E + 8 * n_R_x_E + idx;

                                int idx3_ = ofs_R_x_E + 9 * n_R_x_E + idx;
                                int idx4_ = ofs_R_x_E + 10 * n_R_x_E + idx;
                                int idx5_ = ofs_R_x_E + 11 * n_R_x_E + idx;

                                // TODO: try BLAS here? lol actually not, CBLAS was slower here

                                // C
                                R_x_E[idx0_] += (c / q) * R_x_E_P_R_[idx0] + R_x_E_P_R_[idx3];
                                R_x_E[idx1_] += (c / q) * R_x_E_P_R_[idx1] + R_x_E_P_R_[idx4];
                                R_x_E[idx2_] += (c / q) * R_x_E_P_R_[idx2] + R_x_E_P_R_[idx5];

                                // D
                                R_x_E[idx3_] += (d / q) * R_x_E_P_R_[idx0] - R_x_E_P_R_[idx3];
                                R_x_E[idx4_] += (d / q) * R_x_E_P_R_[idx1] - R_x_E_P_R_[idx4];
                                R_x_E[idx5_] += (d / q) * R_x_E_P_R_[idx2] - R_x_E_P_R_[idx5];
                            }
                }
        }

    array<vec4d, 12> eri4_batch;
    for (int ideriv = 0; ideriv < 12; ideriv++)
        eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

    vector<double> E_x_R_x_E_PR(6 * n_sph_abcd, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            int ofs_R_x_E = iab * 12 * n_R_x_E;

            int ofs_ecoeffs0_ab = iab * n_ecoeffs_ab;
            int ofs_ecoeffs1_ab = iab * 3 * n_ecoeffs_ab;

            std::fill(E_x_R_x_E_PR.begin(), E_x_R_x_E_PR.end(), 0);

            // P
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 0 * n_R_x_E],
                        n_sph_cd, 1.0, &E_x_R_x_E_PR[0 * n_sph_abcd], n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 1 * n_R_x_E],
                        n_sph_cd, 1.0, &E_x_R_x_E_PR[1 * n_sph_abcd], n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 2 * n_R_x_E],
                        n_sph_cd, 1.0, &E_x_R_x_E_PR[2 * n_sph_abcd], n_sph_cd);

            // R
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs1_ab[ofs_ecoeffs1_ab + 0 * n_ecoeffs_ab], n_hermite_ab,
                        &R_x_E[ofs_R_x_E + 3 * n_R_x_E], n_sph_cd, 1.0, &E_x_R_x_E_PR[3 * n_sph_abcd],
                        n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs1_ab[ofs_ecoeffs1_ab + 1 * n_ecoeffs_ab], n_hermite_ab,
                        &R_x_E[ofs_R_x_E + 4 * n_R_x_E], n_sph_cd, 1.0, &E_x_R_x_E_PR[4 * n_sph_abcd],
                        n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs1_ab[ofs_ecoeffs1_ab + 2 * n_ecoeffs_ab], n_hermite_ab,
                        &R_x_E[ofs_R_x_E + 5 * n_R_x_E], n_sph_cd, 1.0, &E_x_R_x_E_PR[5 * n_sph_abcd],
                        n_sph_cd);

            // C
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 6 * n_R_x_E],
                        n_sph_cd, 1.0, eri4_batch[6].memptr(), n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 7 * n_R_x_E],
                        n_sph_cd, 1.0, eri4_batch[7].memptr(), n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 8 * n_R_x_E],
                        n_sph_cd, 1.0, eri4_batch[8].memptr(), n_sph_cd);

            // D
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 9 * n_R_x_E],
                        n_sph_cd, 1.0, eri4_batch[9].memptr(), n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 10 * n_R_x_E],
                        n_sph_cd, 1.0, eri4_batch[10].memptr(), n_sph_cd);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0_ab], n_hermite_ab, &R_x_E[ofs_R_x_E + 11 * n_R_x_E],
                        n_sph_cd, 1.0, eri4_batch[11].memptr(), n_sph_cd);

            // PR -> AB
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;
            cblas_daxpy(n_sph_abcd, (a / p), &E_x_R_x_E_PR[0 * n_sph_abcd], 1, eri4_batch[0].memptr(), 1);
            cblas_daxpy(n_sph_abcd, (a / p), &E_x_R_x_E_PR[1 * n_sph_abcd], 1, eri4_batch[1].memptr(), 1);
            cblas_daxpy(n_sph_abcd, (a / p), &E_x_R_x_E_PR[2 * n_sph_abcd], 1, eri4_batch[2].memptr(), 1);
            cblas_daxpy(n_sph_abcd, 1, &E_x_R_x_E_PR[3 * n_sph_abcd], 1, eri4_batch[0].memptr(), 1);
            cblas_daxpy(n_sph_abcd, 1, &E_x_R_x_E_PR[4 * n_sph_abcd], 1, eri4_batch[1].memptr(), 1);
            cblas_daxpy(n_sph_abcd, 1, &E_x_R_x_E_PR[5 * n_sph_abcd], 1, eri4_batch[2].memptr(), 1);

            cblas_daxpy(n_sph_abcd, (b / p), &E_x_R_x_E_PR[0 * n_sph_abcd], 1, eri4_batch[3].memptr(), 1);
            cblas_daxpy(n_sph_abcd, (b / p), &E_x_R_x_E_PR[1 * n_sph_abcd], 1, eri4_batch[4].memptr(), 1);
            cblas_daxpy(n_sph_abcd, (b / p), &E_x_R_x_E_PR[2 * n_sph_abcd], 1, eri4_batch[5].memptr(), 1);
            cblas_daxpy(n_sph_abcd, -1, &E_x_R_x_E_PR[3 * n_sph_abcd], 1, eri4_batch[3].memptr(), 1);
            cblas_daxpy(n_sph_abcd, -1, &E_x_R_x_E_PR[4 * n_sph_abcd], 1, eri4_batch[4].memptr(), 1);
            cblas_daxpy(n_sph_abcd, -1, &E_x_R_x_E_PR[5 * n_sph_abcd], 1, eri4_batch[5].memptr(), 1);
        }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab + 0];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd + 0];
    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
    const double *norms_a = &sp_data_ab.norms[ofs_norm_a];
    const double *norms_b = &sp_data_ab.norms[ofs_norm_b];
    const double *norms_c = &sp_data_cd.norms[ofs_norm_c];
    const double *norms_d = &sp_data_cd.norms[ofs_norm_d];

    for (int ideriv = 0; ideriv < 12; ideriv++)
        for (int mu = 0; mu < n_sph_a; mu++)
            for (int nu = 0; nu < n_sph_b; nu++)
                for (int ka = 0; ka < n_sph_c; ka++)
                    for (int ta = 0; ta < n_sph_d; ta++)
                    {
                        double norm_a = norms_a[mu];
                        double norm_b = norms_b[nu];
                        double norm_c = norms_c[ka];
                        double norm_d = norms_d[ta];

                        eri4_batch[ideriv](mu, nu, ka, ta) *= norm_a * norm_b * norm_c * norm_d;
                    }

    return eri4_batch;
}

array<lible::vec4d, 12> LIT::kernelERI4Deriv1Test(const int ipair_ab, const int ipair_cd,
                                                  const vector<double> &ecoeffs0_bra_dummy,
                                                  const vector<double> &ecoeffs1_bra_dummy,
                                                  const vector<double> &ecoeffs0_ket_dummy,
                                                  const vector<double> &ecoeffs1_ket_dummy,
                                                  const BoysGrid &boys_grid,
                                                  const ShellPairData &sp_data_ab,
                                                  const ShellPairData &sp_data_cd)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sp_data_cd.la;
    int ld = sp_data_cd.lb;
    int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
    int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
    int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
    int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
    int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
    int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
    // int eofs0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    // int eofs1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    // int eofs0_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];
    // int eofs1_cd = sp_data_cd.offsets_ecoeffs_deriv1[ipair_cd];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sp_data_cd.exps[cofs_c];
    const double *exps_d = &sp_data_cd.exps[cofs_d];
    const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
    const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
    const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];

    vector<double> ecoeffs0_bra = ecoeffsSHARKShellPair(ipair_ab, sp_data_ab);
    vector<double> ecoeffs1_bra = ecoeffsD1SHARKShellPair(ipair_ab, sp_data_ab);
    vector<double> ecoeffs0_ket = ecoeffsSHARKShellPair(ipair_cd, sp_data_cd);
    vector<double> ecoeffs1_ket = ecoeffsD1SHARKShellPair(ipair_cd, sp_data_cd);

    int lab = la + lb;
    int lcd = lc + ld;
    int labcd = lab + lcd;

    vector<array<int, 3>> idxs_tuv_ab = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> idxs_tuv_cd = getHermiteGaussianIdxs(lcd);

    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_c = numSphericals(lc);
    int n_sph_d = numSphericals(ld);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_sph_cd = n_sph_c * n_sph_d;
    int n_sph_abcd = n_sph_ab * n_sph_cd;
    int n_hermite_ab = numHermites(lab);
    int n_hermite_cd = numHermites(lcd);
    // int n_hermite_abcd = n_hermite_ab * n_hermite_cd;
    // int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    // int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;


    int n_r_rows = (cdepth_a * cdepth_b * n_hermite_ab);
    int n_r_cols = (cdepth_c * cdepth_d * n_hermite_cd);
    int n_rints = n_r_rows * n_r_cols;
    vector<double> ecoeffs0_bra_ap(n_sph_ab * n_r_rows);
    vector<double> ecoeffs0_ket_cq(n_sph_cd * n_r_cols);
    vector<double> rints(8 * n_rints);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;
            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double c = exps_c[ic];
                    double d = exps_d[id];
                    
                    double q = c + d;
                    double alpha = p * q / (p + q);

                    array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                           (a * xyz_a[1] + b * xyz_b[1]) / p,
                                           (a * xyz_a[2] + b * xyz_b[2]) / p};

                    array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                           (c * xyz_c[1] + d * xyz_d[1]) / q,
                                           (c * xyz_c[2] + d * xyz_d[2]) / q};

                    array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                            xyz_p[1] - xyz_q[1],
                                            xyz_p[2] - xyz_q[2]};

                    double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                    double x = alpha * (xx * xx + xy * xy + xz * xz);
                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));

                    vector<double> fnx = calcBoysF(labcd + 1, x, boys_grid);

                    int ofs_row = iab * n_hermite_ab;
                    int ofs_col = icd * n_hermite_cd;

                    calcRInts_ERI4_deriv1_test(labcd, fac, alpha, xyz_pq.data(), fnx.data(),
                                               n_rints, ofs_row, ofs_col, n_r_cols, n_r_rows,
                                               idxs_tuv_ab, idxs_tuv_cd, rints.data());
                }

            for (int munu = 0; munu < n_sph_ab; munu++)
            {
                int ofs = munu * n_r_rows + iab * n_hermite_ab;
                for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                    ecoeffs0_bra_ap[ofs + tuv] = (a / p) * ecoeffs0_bra[ofs + tuv];
            }
        }

    for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
        for (int id = 0; id < cdepth_d; id++, icd++)
        {
            double c = exps_c[ic];
            double d = exps_d[id];
            double q = c + d;

            for (int kata = 0; kata < n_sph_cd; kata++)
            {
                int ofs = kata * n_r_cols + icd * n_hermite_cd;
                for (int tuv = 0; tuv < n_hermite_cd; tuv++)
                    ecoeffs0_ket_cq[ofs + tuv] = (c / q) * ecoeffs0_ket[ofs + tuv];
            }
        }

    // SHARK integrals

    int n_R_x_E = n_r_rows * n_sph_cd;
    int n_E_x_R = n_sph_ab * n_r_cols;
    vector<double> R_x_E(4 * n_R_x_E, 0);
    vector<double> E_x_R(4 * n_E_x_R, 0);

    int m = 4 * n_r_rows;
    int n = n_sph_cd;
    int k = n_r_cols;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0, &rints[0 * n_rints], k,
                &ecoeffs0_ket[0], k, 1.0, &R_x_E[0], n);

    m = 4 * n_r_cols;
    n = n_sph_ab;
    k = n_r_rows;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0, &rints[4 * n_rints], k,
                &ecoeffs0_bra[0], k, 1.0, &E_x_R[0], n);

    vector<double> eri4_batch_raw(12 * n_sph_abcd, 0);

    // bra P
    m = n_sph_ab;
    n = n_sph_cd;
    k = n_r_rows;
    vector<double> eri4_batch_P(3 * n_sph_abcd, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra[0], k,
                &R_x_E[0 * n_R_x_E], n, 1.0, &eri4_batch_P[0 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra[0], k,
                &R_x_E[1 * n_R_x_E], n, 1.0, &eri4_batch_P[1 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra[0], k,
                &R_x_E[2 * n_R_x_E], n, 1.0, &eri4_batch_P[2 * n_sph_abcd], n);

    // bra A from P and R
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra_ap[0], k,
                &R_x_E[0 * n_R_x_E], n, 1.0, &eri4_batch_raw[0 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra_ap[0], k,
                &R_x_E[1 * n_R_x_E], n, 1.0, &eri4_batch_raw[1 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra_ap[0], k,
                &R_x_E[2 * n_R_x_E], n, 1.0, &eri4_batch_raw[2 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3 * m, n, k, 1.0, &ecoeffs1_bra[0], k,
                &R_x_E[3 * n_R_x_E], n, 1.0, &eri4_batch_raw[0 * n_sph_abcd], n);

    // B from P and A
    cblas_daxpy(3 * n_sph_abcd, 1.0, &eri4_batch_P[0], 1, &eri4_batch_raw[3 * n_sph_abcd], 1);                 // P
    cblas_daxpy(3 * n_sph_abcd, -1.0, &eri4_batch_raw[0 * n_sph_abcd], 1, &eri4_batch_raw[3 * n_sph_abcd], 1); // A

    // ket Q
    m = n_sph_cd;
    n = n_sph_ab;
    k = n_r_cols;
    vector<double> eri4_batch_Q(3 * n_sph_abcd, 0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket[0], k,
                &E_x_R[0 * n_E_x_R], n, 1.0, &eri4_batch_Q[0 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket[0], k,
                &E_x_R[1 * n_E_x_R], n, 1.0, &eri4_batch_Q[1 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket[0], k,
                &E_x_R[2 * n_E_x_R], n, 1.0, &eri4_batch_Q[2 * n_sph_abcd], n);

    // ket C from Q and S
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket_cq[0], k,
                &E_x_R[0 * n_E_x_R], n, 1.0, &eri4_batch_raw[6 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket_cq[0], k,
                &E_x_R[1 * n_E_x_R], n, 1.0, &eri4_batch_raw[7 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket_cq[0], k,
                &E_x_R[2 * n_E_x_R], n, 1.0, &eri4_batch_raw[8 * n_sph_abcd], n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3 * m, n, k, 1.0, &ecoeffs1_ket[0], k,
                &E_x_R[3 * n_E_x_R], n, 1.0, &eri4_batch_raw[6 * n_sph_abcd], n);

    // D from Q and C
    cblas_daxpy(3 * n_sph_abcd, 1.0, &eri4_batch_Q[0], 1, &eri4_batch_raw[9 * n_sph_abcd], 1);                 // Q
    cblas_daxpy(3 * n_sph_abcd, -1.0, &eri4_batch_raw[6 * n_sph_abcd], 1, &eri4_batch_raw[9 * n_sph_abcd], 1); // C

    array<vec4d, 12> eri4_batch;
    for (int ideriv = 0; ideriv < 12; ideriv++)
        eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab + 0];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd + 0];
    int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
    for (int ideriv = 0; ideriv < 12; ideriv++)
        for (int mu = 0; mu < n_sph_a; mu++)
            for (int nu = 0; nu < n_sph_b; nu++)
                for (int ka = 0; ka < n_sph_c; ka++)
                    for (int ta = 0; ta < n_sph_d; ta++)
                    {
                        double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                        double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                        double norm_c = sp_data_cd.norms[ofs_norm_c + ka];
                        double norm_d = sp_data_cd.norms[ofs_norm_d + ta];

                        int idx;
                        if (ideriv < 6)
                            idx = ideriv * n_sph_abcd + mu * (n_sph_b * n_sph_c * n_sph_d) +
                                  nu * (n_sph_c * n_sph_d) + ka * n_sph_d + ta;
                        else
                            idx = ideriv * n_sph_abcd + ka * (n_sph_d * n_sph_a * n_sph_b) +
                                  ta * (n_sph_a * n_sph_b) + mu * n_sph_b + nu;

                        eri4_batch[ideriv](mu, nu, ka, ta) = norm_a * norm_b * norm_c * norm_d *
                                                             eri4_batch_raw[idx];
                    }

    return eri4_batch;
}