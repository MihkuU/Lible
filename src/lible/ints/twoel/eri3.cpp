#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>
#include <lible/ints/twoel/twoel_detail.hpp>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LIT = lible::ints::two;

using std::array, std::pair, std::vector;

lible::vec3d LIT::calcERI3(const Structure &structure)
{
    if (!structure.getUseRI())
        throw std::runtime_error("RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();
    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellData> sh_datas = shellDatasAux(l_max_aux, structure);
    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

    vector<vector<double>> ecoeffs_aux = ecoeffsSphericalShellDatas_Bra(sh_datas);
    vector<vector<double>> ecoeffs = ecoeffsSphericalSPDatas_Bra(sp_datas);

    size_t dim_ao = structure.getDimAO();
    size_t dim_ao_aux = structure.getDimAOAux();
    vec3d eri3(Fill(0), dim_ao, dim_ao, dim_ao_aux);
    for (size_t lalb = 0; lalb < l_pairs.size(); lalb++)
        for (int lc = 0; lc <= l_max_aux; lc++)
        {
            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sh_data_c = sh_datas[lc];

            auto [la, lb] = l_pairs[lalb];

            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);
            int n_sph_c = numSphericals(lc);

            // const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            // const vector<double> &ecoeffs_c = ecoeffs_aux[lc];
            // vector<double> ecoeffs_ab = ecoeffs[lalb];
            // vector<double> ecoeffs_c = ecoeffs_aux[lc];            
            // std::fill(ecoeffs_ab.begin(), ecoeffs_ab.end(), 1);
            // std::fill(ecoeffs_c.begin(), ecoeffs_c.end(), 1);

            // kernel_eri3_t kernel_eri3 = deployERI3Kernel(la, lb, lc);

            ERI3Kernel eri3_kernel = deployERI3Kernel(sp_data_ab, sh_data_c);

            for (int ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs; ipair_ab++)
                for (int ishell_c = 0; ishell_c < sh_data_c.n_shells; ishell_c++)
                {
                    // vec3d eri3_batch = kernel_eri3(ipair_ab, ishell_c, ecoeffs_ab, ecoeffs_c,
                    //                                sp_data_ab, sh_data_c);

                    vec3d eri3_batch = eri3_kernel(ipair_ab, ishell_c, sp_data_ab,
                                                   sh_data_c);

                    int ofs_a = sp_data_ab.offsets_sph[2 * ipair_ab];
                    int ofs_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
                    int ofs_c = sh_data_c.offsets_sph[ishell_c];
                    for (int ia = 0; ia < n_sph_a; ia++)
                        for (int ib = 0; ib < n_sph_b; ib++)
                            for (int ic = 0; ic < n_sph_c; ic++)
                            {                            
                                int mu = ofs_a + ia;
                                int nu = ofs_b + ib;
                                int ka = ofs_c + ic;

                                eri3(mu, nu, ka) = eri3_batch(ia, ib, ic);
                                eri3(nu, mu, ka) = eri3_batch(ia, ib, ic);
                            }
                }
        }

    return eri3;
}

array<lible::vec3d, 9> LIT::kernelERI3Deriv1(const int ipair_ab, const int ishell_c,
                                             const vector<double> &ecoeffs_ab,
                                             const vector<double> &ecoeffs1_ab,
                                             const vector<double> &ecoeffs_c,
                                             const BoysGrid &boys_grid,
                                             const ShellPairData &sp_data_ab,
                                             const ShellData &sh_data_c)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sh_data_c.l;
    int lab = la + lb;
    int labc = lab + lc;
    int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
    int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
    int cdepth_c = sh_data_c.cdepths[ishell_c];
    int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
    int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
    int cofs_c = sh_data_c.coffsets[ishell_c];
    int eofs_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
    int eofs1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
    int eofs_c = sh_data_c.offsets_ecoeffs[ishell_c];

    const double *exps_a = &sp_data_ab.exps[cofs_a];
    const double *exps_b = &sp_data_ab.exps[cofs_b];
    const double *exps_c = &sh_data_c.exps[cofs_c];
    const double *coords_a = &sp_data_ab.coords[6 * ipair_ab + 0];
    const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
    const double *coords_c = &sh_data_c.coords[3 * ishell_c];
    const double *pecoeffs_ab = &ecoeffs_ab[eofs_ab];
    const double *pecoeffs_deriv1_ab = &ecoeffs1_ab[eofs1_ab];
    const double *pecoeffs_c = &ecoeffs_c[eofs_c];

    vector<array<int, 3>> idxs_tuv_ab = getHermiteGaussianIdxs(lab);
    vector<array<int, 3>> idxs_tuv_c = getHermiteGaussianIdxs(lc);

    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_c = numSphericals(lc);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_sph_abc = n_sph_ab * n_sph_c;
    int n_hermite_ab = numHermites(lab);
    int n_hermite_c = numHermites(lc);
    int n_hermite_abc = n_hermite_ab * n_hermite_c;
    int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    int n_ecoeffs_c = n_sph_c * n_hermite_c;

    int n_R_x_E = n_hermite_ab * n_sph_c;
    // int n_E_x_R = n_sph_ab * n_hermite_c;

    vector<double> R_x_E(7 * cdepth_a * cdepth_b * n_R_x_E, 0);
    // vector<double> E_x_R(3 * cdepth_c * n_E_x_R, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
            for (int ic = 0; ic < cdepth_c; ic++)
            {
                double a = exps_a[ia];
                double b = exps_b[ib];
                double c = exps_c[ic];

                double p = a + b;
                double alpha = p * c / (p + c);

                std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                            (a * coords_a[1] + b * coords_b[1]) / p,
                                            (a * coords_a[2] + b * coords_b[2]) / p};

                std::array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                             xyz_p[1] - coords_c[1],
                                             xyz_p[2] - coords_c[2]};

                double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                double x = alpha * (dx * dx + dy * dy + dz * dz);
                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));

                vector<double> fnx = calcBoysF(labc + 1, x, boys_grid);

                vector<double> rints = calcRInts_ERI3_deriv1(labc, fac, alpha, xyz_pc.data(),
                                                             fnx.data(), idxs_tuv_ab, idxs_tuv_c);

                int ofs_R_x_E = 7 * iab * n_R_x_E;
                int ofs_E_c = ic * n_ecoeffs_c;

                // d/dP + d/dR
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 4 * n_hermite_ab, n_sph_c, n_hermite_c,
                            1.0, &rints[0], n_hermite_c, &pecoeffs_c[ofs_E_c], n_hermite_c,
                            1.0, &R_x_E[ofs_R_x_E + 0 * n_R_x_E], n_sph_c);

                // d/dC
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_hermite_ab, n_sph_c, n_hermite_c,
                            1.0, &rints[4 * n_hermite_abc], n_hermite_c, &pecoeffs_c[ofs_E_c],
                            n_hermite_c, 1.0, &R_x_E[ofs_R_x_E + 4 * n_R_x_E], n_sph_c);

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_hermite_ab, n_sph_c, n_hermite_c,
                            1.0, &rints[5 * n_hermite_abc], n_hermite_c, &pecoeffs_c[ofs_E_c],
                            n_hermite_c, 1.0, &R_x_E[ofs_R_x_E + 5 * n_R_x_E], n_sph_c);

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_hermite_ab, n_sph_c, n_hermite_c,
                            1.0, &rints[6 * n_hermite_abc], n_hermite_c, &pecoeffs_c[ofs_E_c],
                            n_hermite_c, 1.0, &R_x_E[ofs_R_x_E + 6 * n_R_x_E], n_sph_c);

                // // d/dP + d/dR
                // int ofs_E_c = ic * n_ecoeffs_c;
                // int ofs_R_x_E = 4 * iab * n_R_x_E;
                // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 4 * n_hermite_ab, n_sph_c,
                //             n_hermite_c, 1.0, &rints[0], n_hermite_c, &pecoeffs_c[ofs_E_c],
                //             n_hermite_c, 1.0, &R_x_E[ofs_R_x_E], n_sph_c);

                // // d/dC
                // int ofs_E_ab = iab * n_ecoeffs_ab;
                // int ofs_E_x_R = 3 * ic * n_E_x_R;
                // int ofs_R = 4 * n_hermite_abc;
                // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, 3 * n_hermite_c,
                //             n_hermite_ab, 1.0, &ecoeffs_ab[ofs_E_ab], n_hermite_ab, &rints[ofs_R],
                //             3 * n_hermite_c, 1.0, &E_x_R[ofs_E_x_R], 3 * n_hermite_c);
            }

    array<vec3d, 9> eri3_batch;
    for (int ideriv = 0; ideriv < 9; ideriv++)
        eri3_batch[ideriv] = vec3d(Fill(0), n_sph_a, n_sph_b, n_sph_c);

    vector<double> eri3_batch_PR(6 * n_sph_abc, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            std::fill(eri3_batch_PR.begin(), eri3_batch_PR.end(), 0);

            double a = exps_a[ia];
            double b = exps_b[ib];
            double p = a + b;

            int start = iab * 7 * n_R_x_E;
            int ofs0 = start + 0 * n_R_x_E;
            int ofs1 = start + 1 * n_R_x_E;
            int ofs2 = start + 2 * n_R_x_E;
            int ofs3 = start + 3 * n_R_x_E;
            int ofs4 = start + 4 * n_R_x_E;
            int ofs5 = start + 5 * n_R_x_E;
            int ofs6 = start + 6 * n_R_x_E;

            int ofs_ecoeffs0 = iab * n_ecoeffs_ab;
            int ofs_ecoeffs1 = 3 * iab * n_ecoeffs_ab;

            // P
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0], n_hermite_ab, &R_x_E[ofs0], n_sph_c, 1.0,
                        &eri3_batch_PR[0 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0], n_hermite_ab, &R_x_E[ofs1], n_sph_c, 1.0,
                        &eri3_batch_PR[1 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0], n_hermite_ab, &R_x_E[ofs2], n_sph_c, 1.0,
                        &eri3_batch_PR[2 * n_sph_abc], n_sph_c);

            // R
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3 * n_sph_ab, n_sph_c, n_hermite_ab,
                        1.0, &pecoeffs_deriv1_ab[ofs_ecoeffs1], n_hermite_ab, &R_x_E[ofs3], n_sph_c, 1.0,
                        &eri3_batch_PR[3 * n_sph_abc], n_sph_c);

            // C
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c, n_hermite_ab, 
                        1.0, &pecoeffs_ab[ofs_ecoeffs0], n_hermite_ab, &R_x_E[ofs4], n_sph_c, 1.0, 
                        eri3_batch[6].memptr(), n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c, n_hermite_ab, 
                        1.0, &pecoeffs_ab[ofs_ecoeffs0], n_hermite_ab, &R_x_E[ofs5], n_sph_c, 1.0, 
                        eri3_batch[7].memptr(), n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c, n_hermite_ab,
                        1.0, &pecoeffs_ab[ofs_ecoeffs0], n_hermite_ab, &R_x_E[ofs6], n_sph_c, 1.0,
                        eri3_batch[8].memptr(), n_sph_c);

            // // C
            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_sph_ab, n_sph_c, n_hermite_c,
            //             1.0, &R_x_E[ofs4], n_hermite_c, &R_x_E[ofs4], n_sph_c, 1.0,
            //             eri3_batch[6].getData(), n_sph_c);

            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_sph_ab, n_sph_c, n_hermite_c,
            //             1.0, &R_x_E[ofs5], n_hermite_c, &R_x_E[ofs5], n_sph_c, 1.0,
            //             eri3_batch[7].getData(), n_sph_c);

            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_sph_ab, n_sph_c, n_hermite_c,
            //             1.0, &R_x_E[ofs6], n_hermite_c, &R_x_E[ofs6], n_sph_c, 1.0,
            //             eri3_batch[8].getData(), n_sph_c);

            // PR -> AB
            for (int mu = 0; mu < n_sph_a; mu++)
                for (int nu = 0; nu < n_sph_b; nu++)
                    for (int ka = 0; ka < n_sph_c; ka++)
                    {
                        int munu = mu * n_sph_b + nu;
                        int munuka = munu * n_sph_c + ka;

                        int idx0 = 0 * n_sph_abc + munuka;
                        int idx1 = 1 * n_sph_abc + munuka;
                        int idx2 = 2 * n_sph_abc + munuka;

                        int idx3 = 3 * n_sph_abc + munuka;
                        int idx4 = 4 * n_sph_abc + munuka;
                        int idx5 = 5 * n_sph_abc + munuka;

                        // A
                        eri3_batch[0](mu, nu, ka) += (a / p) * eri3_batch_PR[idx0] + eri3_batch_PR[idx3];
                        eri3_batch[1](mu, nu, ka) += (a / p) * eri3_batch_PR[idx1] + eri3_batch_PR[idx4];
                        eri3_batch[2](mu, nu, ka) += (a / p) * eri3_batch_PR[idx2] + eri3_batch_PR[idx5];

                        // B
                        eri3_batch[3](mu, nu, ka) += (b / p) * eri3_batch_PR[idx0] - eri3_batch_PR[idx3];
                        eri3_batch[4](mu, nu, ka) += (b / p) * eri3_batch_PR[idx1] - eri3_batch_PR[idx4];
                        eri3_batch[5](mu, nu, ka) += (b / p) * eri3_batch_PR[idx2] - eri3_batch_PR[idx5];
                    }
        }

    // for (int ic = 0; ic < cdepth_c; ic++)
    // {
    //     int ofs_E_c = ic * n_ecoeffs_c;
    //     int ofs0 = (3 * ic + 0) * n_E_x_R;
    //     int ofs1 = (3 * ic + 1) * n_E_x_R;
    //     int ofs2 = (3 * ic + 2) * n_E_x_R;

    //     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_sph_ab, n_sph_c, n_hermite_c,
    //                 1.0, &E_x_R[ofs0], n_hermite_c, &pecoeffs_c[ofs_E_c], n_hermite_c, 1.0,
    //                 eri3_batch[6].getData(), n_sph_c);

    //     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_sph_ab, n_sph_c, n_hermite_c,
    //                 1.0, &E_x_R[ofs1], n_hermite_c, &pecoeffs_c[ofs_E_c], n_hermite_c, 1.0,
    //                 eri3_batch[7].getData(), n_sph_c);

    //     cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_sph_ab, n_sph_c, n_hermite_c,
    //                 1.0, &E_x_R[ofs2], n_hermite_c, &pecoeffs_c[ofs_E_c], n_hermite_c, 1.0,
    //                 eri3_batch[8].getData(), n_sph_c);

    //     // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_hermite_ab, n_sph_c, n_hermite_c,
    //     //             1.0, &rints[6 * n_hermite_abc], n_hermite_c, &pecoeffs_c[ofs_ecoeffs_c],
    //     //             n_hermite_c, 1.0, &R_x_E[ofs_R_x_E + 6 * n_R_x_E], n_sph_c);
    // }

    int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab + 0];
    int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
    int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];

    const double *norms_a = &sp_data_ab.norms[ofs_norm_a];
    const double *norms_b = &sp_data_ab.norms[ofs_norm_b];
    const double *norms_c = &sh_data_c.norms[ofs_norm_c];

    for (int ideriv = 0; ideriv < 9; ideriv++)
        for (int mu = 0; mu < n_sph_a; mu++)
            for (int nu = 0; nu < n_sph_b; nu++)
                for (int ka = 0; ka < n_sph_c; ka++)
                {
                    double norm_a = norms_a[mu];
                    double norm_b = norms_b[nu];
                    double norm_c = norms_c[ka];

                    eri3_batch[ideriv](mu, nu, ka) *= norm_a * norm_b * norm_c;
                }

    return eri3_batch;
}