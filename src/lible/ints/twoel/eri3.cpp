#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>

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

    vector<vector<double>> ecoeffs_aux = ecoeffsSphericalShellDatas_Bra(l_max_aux, sh_datas);
    vector<vector<double>> ecoeffs = ecoeffsSphericalSPDatas_Bra(l_pairs, sp_datas);

    size_t dim_ao = structure.getDimAO();
    size_t dim_ao_aux = structure.getDimAOAux();
    vec3d eri3(dim_ao, dim_ao, dim_ao_aux, 0);
    for (size_t lalb = 0; lalb < l_pairs.size(); lalb++)
        for (int lc = 0; lc <= l_max_aux; lc++)
        {
            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sh_data_c = sh_datas[lc];

            auto [la, lb] = l_pairs[lalb];

            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);
            int n_sph_c = numSphericals(lc);
            int n_sph_ab = n_sph_a * n_sph_b;

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_c = ecoeffs_aux[lc];

            kernel_eri3_t kernel_eri3 = deployERI3Kernel(la, lb, lc);

            vector<double> eri3_batch(n_sph_ab * n_sph_c, 0);
            for (int ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs; ipair_ab++)
                for (int ishell_c = 0; ishell_c < sh_data_c.n_shells; ishell_c++)
                {
                    int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
                    int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
                    int pos_c = sh_data_c.coffsets[ishell_c];

                    kernel_eri3(sp_data_ab.cdepths[2 * ipair_ab],
                                sp_data_ab.cdepths[2 * ipair_ab + 1],
                                sh_data_c.cdepths[ishell_c],
                                &sp_data_ab.exps[pos_a],
                                &sp_data_ab.exps[pos_b],
                                &sh_data_c.exps[pos_c],
                                &sp_data_ab.coords[6 * ipair_ab],
                                &sp_data_ab.coords[6 * ipair_ab + 3],
                                &sh_data_c.coords[3 * ishell_c],
                                &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]],
                                &ecoeffs_c[sh_data_c.offsets_ecoeffs[ishell_c]],
                                &eri3_batch[0]);

                    transferIntsERI3(ipair_ab, ishell_c, sh_data_c, sp_data_ab, eri3_batch, eri3);
                }
        }

    return eri3;
}

void LIT::kernelERI3Deriv1(const int la, const int lb, const int lc,
                           const int cdepth_a, const int cdepth_b, const int cdepth_c,
                           const double *exps_a, const double *exps_b, const double *exps_c,
                           const double *coords_a, const double *coords_b, const double *coords_c,
                           const double *ecoeffs_ab, const double *ecoeffs_deriv1_ab,
                           const double *ecoeffs_c, const double *norms_a,
                           const double *norms_b, const double *norms_c,
                           const BoysGrid &boys_grid, double *eri3_batch)
{
    int lab = la + lb;
    int labc = lab + lc;

    vector<array<int, 3>> idxs_tuv_ab = returnHermiteGaussianIdxs(lab);
    vector<array<int, 3>> idxs_tuv_c = returnHermiteGaussianIdxs(lc);

    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_c = numSphericals(lc);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_sph_abc = n_sph_ab * n_sph_c;
    int n_hermite_ab = numHermites(lab);
    int n_hermite_c = numHermites(lc);
    int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
    int n_ecoeffs_c = n_sph_c * n_hermite_c;

    std::fill(eri3_batch, eri3_batch + 9 * n_sph_abc, 0);

    int n_R_x_E = n_hermite_ab * n_sph_c;
    vector<double> R_x_E(7 * cdepth_a * cdepth_b * n_R_x_E, 0);
    for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
        for (int ib = 0; ib < cdepth_b; ib++, iab++)
        {
            int ofs_R_x_E = 7 * iab * n_R_x_E;
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

                int ofs_ecoeffs_c = ic * n_ecoeffs_c;
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 7 * n_hermite_ab, n_sph_c,
                            n_hermite_c, 1.0, &rints[0], n_hermite_c, &ecoeffs_c[ofs_ecoeffs_c],
                            n_hermite_c, 1.0, &R_x_E[ofs_R_x_E], n_sph_c);
            }
        }

    vector<double> eri3_batch_PR(6 * n_sph_ab * n_sph_c, 0);
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

            int ofs_ecoeffs = iab * n_ecoeffs_ab;
            int ofs_ecoeffs_deriv1 = 3 * iab * n_ecoeffs_ab;

            // P
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_ab[ofs_ecoeffs], n_hermite_ab, &R_x_E[ofs0],
                        n_sph_c, 1.0, &eri3_batch_PR[0 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_ab[ofs_ecoeffs], n_hermite_ab, &R_x_E[ofs1],
                        n_sph_c, 1.0, &eri3_batch_PR[1 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_ab[ofs_ecoeffs], n_hermite_ab, &R_x_E[ofs2],
                        n_sph_c, 1.0, &eri3_batch_PR[2 * n_sph_abc], n_sph_c);

            // R
            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3 * n_sph_ab, n_sph_c,
            //             n_hermite_ab, 1.0, &ecoeffs_deriv1_ab[ofs_ecoeffs_deriv1], n_hermite_ab,
            //             &R_x_E[ofs3], n_sph_c, 1.0, &eri3_batch_PR[3 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_deriv1_ab[ofs_ecoeffs_deriv1 + 0 * n_ecoeffs_ab], n_hermite_ab,
                        &R_x_E[ofs3], n_sph_c, 1.0, &eri3_batch_PR[3 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_deriv1_ab[ofs_ecoeffs_deriv1 + 1 * n_ecoeffs_ab], n_hermite_ab,
                        &R_x_E[ofs3], n_sph_c, 1.0, &eri3_batch_PR[4 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_deriv1_ab[ofs_ecoeffs_deriv1 + 2 * n_ecoeffs_ab], n_hermite_ab,
                        &R_x_E[ofs3], n_sph_c, 1.0, &eri3_batch_PR[5 * n_sph_abc], n_sph_c);

            // C
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_ab[ofs_ecoeffs], n_hermite_ab, &R_x_E[ofs4],
                        n_sph_c, 1.0, &eri3_batch[6 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_ab[ofs_ecoeffs], n_hermite_ab, &R_x_E[ofs5],
                        n_sph_c, 1.0, &eri3_batch[7 * n_sph_abc], n_sph_c);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                        n_hermite_ab, 1.0, &ecoeffs_ab[ofs_ecoeffs], n_hermite_ab, &R_x_E[ofs6],
                        n_sph_c, 1.0, &eri3_batch[8 * n_sph_abc], n_sph_c);

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

                        // TODO: try BLAS here?

                        // A
                        eri3_batch[idx0] += (a / p) * eri3_batch_PR[idx0] + eri3_batch_PR[idx3];
                        eri3_batch[idx1] += (a / p) * eri3_batch_PR[idx1] + eri3_batch_PR[idx4];
                        eri3_batch[idx2] += (a / p) * eri3_batch_PR[idx2] + eri3_batch_PR[idx5];

                        // B
                        eri3_batch[idx3] += (b / p) * eri3_batch_PR[idx0] - eri3_batch_PR[idx3];
                        eri3_batch[idx4] += (b / p) * eri3_batch_PR[idx1] - eri3_batch_PR[idx4];
                        eri3_batch[idx5] += (b / p) * eri3_batch_PR[idx2] - eri3_batch_PR[idx5];
                    }
        }

    for (int ideriv = 0; ideriv < 9; ideriv++)
    {
        int ofs = ideriv * n_sph_abc;
        for (int a = 0, ab = 0; a < n_sph_a; a++)
            for (int b = 0; b < n_sph_b; b++, ab++)
                for (int c = 0; c < n_sph_c; c++)
                {
                    int abc = ab * n_sph_c + c;
                    int idx = ofs + abc;

                    double norm_a = norms_a[a];
                    double norm_b = norms_b[b];
                    double norm_c = norms_c[c];

                    eri3_batch[idx] *= norm_a * norm_b * norm_c;
                }
    }
}