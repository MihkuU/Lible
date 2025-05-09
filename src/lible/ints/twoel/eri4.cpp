#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/utils.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/utils.hpp>

#include <cstring>
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

    void kernelERI4Deriv1(const int lab, const int lcd, const int ipair_ab, const int ipair_cd,
                          const vector<array<int, 3>> &idxs_tuv_ab,
                          const vector<array<int, 3>> &idxs_tuv_cd,
                          const vector<double> &ecoeffs_ab,
                          const vector<double> &ecoeffs_deriv1_ab,
                          const vector<double> &ecoeffs_cd_tsp,
                          const vector<double> &ecoeffs_deriv1_cd_tsp,
                          const BoysGrid &boys_grid, const ShellPairData &sp_data_ab,
                          const ShellPairData &sp_data_cd)
    {
        int labcd = lab + lcd;

        int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
        int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
        int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
        int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
        int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
        int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
        int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
        int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

        const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab];
        const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];    
        const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd];
        const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];

        int n_sph_a = numSphericals(sp_data_ab.la);
        int n_sph_b = numSphericals(sp_data_ab.lb);
        int n_sph_c = numSphericals(sp_data_cd.la);
        int n_sph_d = numSphericals(sp_data_cd.lb);
        int n_tuv_ab = numHermites(lab);
        int n_tuv_cd = numHermites(lcd);
        int n_sph_ab = n_sph_a * n_sph_b;
        int n_sph_cd = n_sph_c * n_sph_d;
        int n_ecoeffs_ab = n_sph_ab * n_tuv_ab;
        int n_ecoeffs_cd = n_sph_cd * n_tuv_cd;
        int n_rints_x_ecoeffs = n_tuv_ab * n_sph_cd;
        
        vector<double> ecoeffs_ket(12 * n_ecoeffs_cd, 0);
        vector<double> T_PRP_R_(cdepth_a * cdepth_b * 12 * n_rints_x_ecoeffs, 0);
        vector<double> T_PRCD(cdepth_a * cdepth_b * 12 * n_rints_x_ecoeffs, 0);
        for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
            for (int id = 0; id < cdepth_d; id++, icd++)
            {                
                int ofs_ecoeffs_cd = sp_data_cd.offsets_ecoeffs[ipair_cd] + icd * n_ecoeffs_cd;
                for (int i = 0; i < 9; i++)
                    std::memcpy(ecoeffs_ket.data() + i * n_ecoeffs_cd, &ecoeffs_ket[ofs_ecoeffs_cd],
                                sizeof(double) * n_ecoeffs_cd);

                int ofs_ecoeffs_deriv1_cd = sp_data_cd.offsets_ecoeffs_deriv1[ipair_cd] + 3 * icd * n_ecoeffs_cd;
                std::memcpy(ecoeffs_ket.data() + 9 * n_ecoeffs_cd,
                            &ecoeffs_deriv1_cd_tsp[ofs_ecoeffs_deriv1_cd],
                            sizeof(double) * 3 * n_ecoeffs_cd);

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
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

                        array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],
                                                xyz_p[2] - xyz_q[2]};

                        double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                        double x = alpha * (xx * xx + xy * xy + xz * xz);                        

                        // Gota increment l by 1 for the R(t + t' + 1, u + u', v + v')...
                        vector<double> fnx = calcBoysF(labcd + 1, x, boys_grid); 

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));

                        vector<double> rints = calcRInts_ERI4_Deriv1(labcd, fac, alpha, xyz_pq.data(),
                                                                     fnx.data(), idxs_tuv_ab,
                                                                     idxs_tuv_cd);

                        // P
                        // R 
                        // P'
                        // R'

                        // // Trafo |P'R') -> |CD)
                        // for (int i = 0; i < n_sph_a; i++)
                        //     for (int j = 0; j < n_sph_b; j++)
                        //         for (size_t tuv = 0; tuv < idxs_tuv_cd.size(); tuv++)
                        //         {
                        //             int ofs_P_ = ;
                        //             int ofs_R_;
                        //         }
                    }
            }

        vector<double> ecoeffs_bra(12 * n_ecoeffs_ab, 0);
        vector<double> ints_PRCD;
        vector<double> ints_ABCD;        
        for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
            for (int ib = 0; ib < cdepth_b; ib++, iab++)
            {
                int ofs_ecoeffs_ab = sp_data_ab.offsets_ecoeffs[ipair_ab] + iab * n_ecoeffs_ab;
                for (int i = 0; i < 3; i++)
                    std::memcpy(ecoeffs_bra.data() + i * n_ecoeffs_ab, &ecoeffs_ab[ofs_ecoeffs_ab],
                                sizeof(double) * n_ecoeffs_ab);

                for (int i = 6; i < 12; i++)
                    std::memcpy(ecoeffs_bra.data() + i * n_ecoeffs_ab, &ecoeffs_ab[ofs_ecoeffs_ab],
                                sizeof(double) * n_ecoeffs_ab);

                int ofs_ecoeffs_deriv1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab] + 3 * iab * n_ecoeffs_ab;
                std::memcpy(ecoeffs_bra.data() + 3 * n_ecoeffs_ab,
                            &ecoeffs_deriv1_ab[ofs_ecoeffs_deriv1_ab],
                            sizeof(double) * 3 * n_ecoeffs_ab);

                // Trafo (PR| -> (AB|
            }
    }
}

lible::vec4d LIT::calcERI4(const Structure &structure)
{
    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

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
                                   eri4_shells_sph, fnx);

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
                                   eri4_shells_sph, fnx);

                        transferIntsERI4(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                         eri4_shells_sph, eri4);
                    }
        }

    return eri4;
}

lible::vec4d LIT::calcERI4New(const Structure &structure)
{
    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

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

    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

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
                                   eri4_shells_sph, fnx);

                        n_shells_abcd++;
                    }
            else
                for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                    for (int ipair_cd = 0; ipair_cd < n_pairs_cd; ipair_cd++)
                    {
                        vector<double> eri4_shells_sph(dim_ab_sph * dim_cd_sph, 0);
                        kernelERI4(lab, lcd, ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs_cd_tsp,
                                   idxs_tuv_ab, idxs_tuv_cd, sp_data_ab, sp_data_cd, boys_f,
                                   eri4_shells_sph, fnx);

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

    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

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
    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

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
                               sp_data_ab, eri4_shells_sph, fnx);

            transferIntsERI4Diag(ipair_ab, sp_data_ab, eri4_shells_sph, eri4_diagonal);
        }
    }

    return eri4_diagonal;
}

void LIT::kernelERI4Deriv1(const int la, const int lb, const int lc, const int ld,
                           const int cdepth_a, const int cdepth_b, const int cdepth_c,
                           const int cdepth_d, const double *exps_a, const double *exps_b,
                           const double *exps_c, const double *exps_d,
                           const double *xyz_a, const double *xyz_b,
                           const double *xyz_c, const double *xyz_d,
                           const double *ecoeffs_ab, const double *ecoeffs_deriv1_ab,
                           const double *ecoeffs_cd_tsp, const double *ecoeffs_deriv1_cd_tsp,
                           const double *norms_a, const double *norms_b,
                           const double *norms_c, const double *norms_d,
                           const BoysGrid &boys_grid, double *eri4_batch)
{
    int lab = la + lb;
    int lcd = lc + ld;
    int labcd = lab + lcd;

    vector<array<int, 3>> idxs_tuv_ab = returnHermiteGaussianIdxs(lab);
    vector<array<int, 3>> idxs_tuv_cd = returnHermiteGaussianIdxs(lcd);

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
    int n_ecoeffs_cd = n_hermite_ab * n_sph_cd; 

    std::fill(eri4_batch, eri4_batch + 12 * n_sph_abcd, 0);

    int n_R_x_E = n_hermite_ab * n_sph_cd;

    vector<double> R_x_E(12 * cdepth_a * cdepth_b * n_R_x_E, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a ; ia++)
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

                    // P
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[0 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 0 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[1 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 1 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[2 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 2 * n_R_x_E], n_sph_cd);

                    // R
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 3 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 4 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E[ofs_R_x_E + 5 * n_R_x_E], n_sph_cd);

                    // P'
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[4 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[0 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[5 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[1 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[6 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_cd_tsp[icd * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[2 * n_R_x_E], n_sph_cd);

                    // R'
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_deriv1_cd_tsp[(icd + 0) * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[3 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_deriv1_cd_tsp[(icd + 1) * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[4 * n_R_x_E], n_sph_cd);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab, n_sph_cd,
                                n_hermite_cd, 1.0, &rints[3 * n_hermite_abcd], n_hermite_cd,
                                &ecoeffs_deriv1_cd_tsp[(icd + 2) * n_ecoeffs_cd], n_sph_cd, 1.0,
                                &R_x_E_P_R_[5 * n_R_x_E], n_sph_cd);

                    // P'R' -> CD
                    for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                        for (int ka = 0, kata = 0; ka < n_sph_c; ka++)
                            for (int ta = 0; ta < n_sph_d; ta++, kata++)
                            {
                                int idx_e = tuv * n_sph_cd + kata;
 
                                int idx0 = 0 * n_R_x_E + idx_e;
                                int idx1 = 1 * n_R_x_E + idx_e;
                                int idx2 = 2 * n_R_x_E + idx_e;

                                int idx3 = 3 * n_R_x_E + idx_e;
                                int idx4 = 4 * n_R_x_E + idx_e;
                                int idx5 = 5 * n_R_x_E + idx_e;

                                int idx0_ = ofs_R_x_E + 6 * n_R_x_E + idx_e;
                                int idx1_ = ofs_R_x_E + 7 * n_R_x_E + idx_e;
                                int idx2_ = ofs_R_x_E + 8 * n_R_x_E + idx_e;
                                 
                                int idx3_ = ofs_R_x_E + 9 * n_R_x_E + idx_e;
                                int idx4_ = ofs_R_x_E + 10 * n_R_x_E + idx_e;
                                int idx5_ = ofs_R_x_E + 11 * n_R_x_E + idx_e;

                                // TODO: try BLAS here?

                                // A
                                R_x_E[idx0_] += (c / q) * R_x_E_P_R_[idx0] + R_x_E_P_R_[idx3];
                                R_x_E[idx1_] += (c / q) * R_x_E_P_R_[idx1] + R_x_E_P_R_[idx4];
                                R_x_E[idx2_] += (c / q) * R_x_E_P_R_[idx2] + R_x_E_P_R_[idx5];

                                // B
                                R_x_E[idx3_] += (d / q) * R_x_E_P_R_[idx0] - R_x_E_P_R_[idx3];
                                R_x_E[idx4_] += (d / q) * R_x_E_P_R_[idx1] - R_x_E_P_R_[idx4];
                                R_x_E[idx5_] += (d / q) * R_x_E_P_R_[idx2] - R_x_E_P_R_[idx5];
                            }
                }
        }

    vector<double> E_x_R_x_E_PR(6 * n_sph_ab * n_sph_cd, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {

        }
}