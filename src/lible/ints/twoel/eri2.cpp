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

using std::array, std::vector;

namespace lible::ints::two
{
    void kernelERI2Diagonal(const int ishell, const int la, const vector<double> &ecoeffs_a,
                            const vector<double> &ecoeffs_a_tsp,
                            const vector<array<int, 3>> &idxs_tuv_a, const BoysF &boys_f,
                            const ShellData &sh_data_a, vector<double> &fnx,
                            vector<double> &eri2_shells_sph)
    {
        int laa = la + la;
        int dim_a = sh_data_a.cdepths[ishell];
        int pos_a = sh_data_a.coffsets[ishell];

        array<double, 3> xyz_aa{0, 0, 0};

        int dim_sph_a = numSphericals(la);
        int dim_tuv_a = numHermites(la);
        int dim_ecoeffs_a = dim_sph_a * dim_tuv_a;
        int dim_rints_x_ecoeffs = dim_tuv_a * dim_sph_a;
        vector<double> rints_x_ecoeffs(dim_a * dim_rints_x_ecoeffs, 0);

        for (int ia = 0; ia < dim_a; ia++)
        {
            int pos_rints_x_ecoeffs = ia * dim_rints_x_ecoeffs;
            for (int ib = 0; ib < dim_a; ib++)
            {
                double a = sh_data_a.exps[pos_a + ia];
                double b = sh_data_a.exps[pos_a + ib];

                double alpha = a * b / (a + b);
                double x = 0;
                boys_f.calcFnx(laa, x, fnx);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));

                vector<double> rints = calcRIntsMatrix(laa, fac, alpha, xyz_aa.data(), fnx.data(),
                                                       idxs_tuv_a, idxs_tuv_a);

                int pos_ecoeffs_b = sh_data_a.offsets_ecoeffs[ishell] + ib * dim_ecoeffs_a;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_tuv_a, dim_sph_a,
                            dim_tuv_a, 1.0, &rints[0], dim_tuv_a, &ecoeffs_a_tsp[pos_ecoeffs_b],
                            dim_sph_a, 1.0, &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_a);
            }
        }

        for (int ia = 0; ia < dim_a; ia++)
        {
            int pos_rints_x_ecoeffs = ia * dim_rints_x_ecoeffs;
            int pos_ecoeffs_a = sh_data_a.offsets_ecoeffs[ishell] + ia * dim_ecoeffs_a;

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_a, dim_sph_a,
                        dim_tuv_a, 1.0, &ecoeffs_a[pos_ecoeffs_a], dim_tuv_a,
                        &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_a, 1.0, &eri2_shells_sph[0],
                        dim_sph_a);
        }
    }
}

lible::vec2d LIT::calcERI2(const Structure &structure)
{
    if (!structure.getUseRI())
        throw std::runtime_error("RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();
    vector<ShellData> sh_datas = shellDatasAux(l_max_aux, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalShellDatas_BraKet(l_max_aux, sh_datas);

    size_t dim_ao_aux = structure.getDimAOAux();
    vec2d eri2(dim_ao_aux, dim_ao_aux, 0);
    for (int la = 0; la <= l_max_aux; la++)
        for (int lb = 0; lb <= la; lb++)
        {
            const auto &sh_data_a = sh_datas[la];
            const auto &sh_data_b = sh_datas[lb];

            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);

            const vector<double> &ecoeffs_a = ecoeffs[la];
            const vector<double> &ecoeffs_b_tsp = ecoeffs_tsp[lb];

            kernel_eri2_t kernel_eri2 = deployERI2Kernel(la, lb);

            vector<double> eri2_batch(n_sph_a * n_sph_b, 0);
            for (int ishell_a = 0; ishell_a < sh_data_a.n_shells; ishell_a++)
            {
                int bound_b = (la == lb) ? ishell_a + 1 : sh_data_b.n_shells;
                for (int ishell_b = 0; ishell_b < bound_b; ishell_b++)
                {
                    int pos_a = sh_data_a.coffsets[ishell_a];
                    int pos_b = sh_data_b.coffsets[ishell_b];

                    kernel_eri2(sh_data_a.cdepths[ishell_a], sh_data_b.cdepths[ishell_b],
                                &sh_data_a.exps[pos_a], &sh_data_b.exps[pos_b],
                                &sh_data_a.coords[3 * ishell_a], &sh_data_b.coords[3 * ishell_b],
                                &ecoeffs_a[sh_data_a.offsets_ecoeffs[ishell_a]],
                                &ecoeffs_b_tsp[sh_data_b.offsets_ecoeffs[ishell_b]],
                                &eri2_batch[0]);

                    transferIntsERI2(ishell_a, ishell_b, sh_data_a, sh_data_b,
                                     eri2_batch, eri2);
                }
            }
        }

    return eri2;
}

vector<double> LIT::calcERI2Diagonal(const Structure &structure)
{
    if (!structure.getUseRI())
        throw std::runtime_error("RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();
    vector<ShellData> sh_datas = shellDatasAux(l_max_aux, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalShellDatas_BraKet(l_max_aux, sh_datas);

    size_t dim_ao_aux = structure.getDimAOAux();
    vector<double> eri2_diagonal(dim_ao_aux, 0);
    for (int la = 0; la <= l_max_aux; la++)
    {
        const auto &sh_data_a = sh_datas[la];

        int n_shells_a = sh_data_a.n_shells;

        int dim_sph_a = numSphericals(la);

        int laa = la + la;
        BoysF boys_f(laa);

        const vector<double> &ecoeffs_a = ecoeffs[la];
        const vector<double> &ecoeffs_a_tsp = ecoeffs_tsp[la];

        vector<array<int, 3>> idxs_tuv_a = returnHermiteGaussianIdxs(la);

        vector<double> fnx(laa + 1, 0);

        for (int ishell = 0; ishell < n_shells_a; ishell++)
        {
            vector<double> eri2_shells_sph(dim_sph_a * dim_sph_a, 0);

            kernelERI2Diagonal(ishell, la, ecoeffs_a, ecoeffs_a_tsp, idxs_tuv_a, boys_f, sh_data_a,
                               fnx, eri2_shells_sph);

            transferIntsERI2Diag(ishell, sh_data_a, eri2_shells_sph, eri2_diagonal);
        }
    }

    return eri2_diagonal;
}

void LIT::kernelERI2Deriv1(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                           const double *exps_a, const double *exps_b, const double *coords_a,
                           const double *coords_b, const double *ecoeffs_a,
                           const double *ecoeffs_b_tsp, const double *norms_a,
                           const double *norms_b, const BoysGrid &boys_grid,
                           double *eri2_batch)
{
    int lab = la + lb;

    vector<array<int, 3>> idxs_tuv_a = returnHermiteGaussianIdxs(la);
    vector<array<int, 3>> idxs_tuv_b = returnHermiteGaussianIdxs(lb);

    int n_sph_a = numSphericals(la);
    int n_sph_b = numSphericals(lb);
    int n_sph_ab = n_sph_a * n_sph_b;
    int n_hermite_a = numHermites(la);
    int n_hermite_b = numHermites(lb);
    int n_ecoeffs_b = n_sph_b * n_hermite_b;

    std::fill(eri2_batch, eri2_batch + 6 * n_sph_ab, 0);

    array<double, 3> xyz_ab{coords_a[0] - coords_b[0], coords_a[1] - coords_b[1],
                            coords_a[2] - coords_b[2]};
    double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
    double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

    int n_R_x_E = n_hermite_a * n_sph_b;
    vector<double> R_x_E(6 * cdepth_a * n_R_x_E, 0);
    for (int ia = 0; ia < cdepth_a; ia++)
    {
        int ofs_R_x_E = 6 * ia * n_R_x_E;
        for (int ib = 0; ib < cdepth_b; ib++)
        {
            double a = exps_a[ia];
            double b = exps_b[ib];

            double alpha = a * b / (a + b);
            double x = alpha * xyz_ab_dot;
            double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));

            vector<double> fnx = calcBoysF(lab + 1, x, boys_grid);

            vector<double> rints = calcRInts_ERI2_deriv1(lab, fac, alpha, xyz_ab.data(),
                                                         fnx.data(), idxs_tuv_a, idxs_tuv_b);

            int ofs_ecoeffs_b = ib * n_ecoeffs_b;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6 * n_hermite_a,
                        n_sph_b, n_hermite_b, 1.0, &rints[0], n_hermite_b,
                        &ecoeffs_b_tsp[ofs_ecoeffs_b], n_sph_b, 1.0,
                        &R_x_E[ofs_R_x_E], n_sph_b);
        }
    }

    for (int ia = 0; ia < cdepth_a; ia++)
    {
        int start = ia * 6 * n_R_x_E;
        int ofs0 = start + 0 * n_R_x_E;
        int ofs1 = start + 1 * n_R_x_E;
        int ofs2 = start + 2 * n_R_x_E;
        int ofs3 = start + 3 * n_R_x_E;
        int ofs4 = start + 4 * n_R_x_E;
        int ofs5 = start + 5 * n_R_x_E;

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                    1.0, &ecoeffs_a[0], n_hermite_a, &R_x_E[ofs0], n_sph_b, 1.0,
                    &eri2_batch[0 * n_sph_ab], n_sph_b);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                    1.0, &ecoeffs_a[0], n_hermite_a, &R_x_E[ofs1], n_sph_b, 1.0,
                    &eri2_batch[1 * n_sph_ab], n_sph_b);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                    1.0, &ecoeffs_a[0], n_hermite_a, &R_x_E[ofs2], n_sph_b, 1.0,
                    &eri2_batch[2 * n_sph_ab], n_sph_b);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                    1.0, &ecoeffs_a[0], n_hermite_a, &R_x_E[ofs3], n_sph_b, 1.0,
                    &eri2_batch[3 * n_sph_ab], n_sph_b);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                    1.0, &ecoeffs_a[0], n_hermite_a, &R_x_E[ofs4], n_sph_b, 1.0,
                    &eri2_batch[4 * n_sph_ab], n_sph_b);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                    1.0, &ecoeffs_a[0], n_hermite_a, &R_x_E[ofs5], n_sph_b, 1.0,
                    &eri2_batch[5 * n_sph_ab], n_sph_b);
    }

    for (int ideriv = 0; ideriv < 6; ideriv++)
    {
        for (int a = 0; a < n_sph_a; a++)
            for (int b = 0; b < n_sph_b; b++)
            {
                int ab = a * n_sph_b + b;
                int idx = ideriv * n_sph_ab + ab;

                double norm_a = norms_a[a];
                double norm_b = norms_b[b];

                eri2_batch[idx] *= norm_a * norm_b;
            }
    }
}