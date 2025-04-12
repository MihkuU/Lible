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
                            const ShellData &sh_data_a, vector<double> &fnx, vector<double> &rints,
                            vec4d &rints_tmp, vector<double> &eri2_shells_sph)
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
                calcRInts(la, la, fac, alpha, xyz_aa, fnx, idxs_tuv_a, idxs_tuv_a, rints_tmp,
                          rints);

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
    vector<ShellData> sh_datas = constructShellDatasAux(l_max_aux, structure);

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
    vector<ShellData> sh_datas = constructShellDatasAux(l_max_aux, structure);

    auto [ecoeffs, ecoeffs_tsp] = ecoeffsSphericalShellDatas_BraKet(l_max_aux, sh_datas);

    size_t dim_ao_aux = structure.getDimAOAux();
    vector<double> eri2_diagonal(dim_ao_aux, 0);
    for (int la = 0; la <= l_max_aux; la++)
    {
        const auto &sh_data_a = sh_datas[la];

        int n_shells_a = sh_data_a.n_shells;

        int dim_sph_a = numSphericals(la);
        int dim_tuv_a = numHermites(la);

        int laa = la + la;
        BoysF boys_f(laa);

        const vector<double> &ecoeffs_a = ecoeffs[la];
        const vector<double> &ecoeffs_a_tsp = ecoeffs_tsp[la];

        vector<array<int, 3>> idxs_tuv_a = returnHermiteGaussianIdxs(la);

        vector<double> rints(std::pow(dim_tuv_a, 2), 0);
        vector<double> fnx(laa + 1, 0);
        vec4d rints_tmp(laa + 1, 0);

        for (int ishell = 0; ishell < n_shells_a; ishell++)
        {
            vector<double> eri2_shells_sph(dim_sph_a * dim_sph_a, 0);

            kernelERI2Diagonal(ishell, la, ecoeffs_a, ecoeffs_a_tsp, idxs_tuv_a, boys_f, sh_data_a,
                               fnx, rints, rints_tmp, eri2_shells_sph);

            transferIntsERI2Diag(ishell, sh_data_a, eri2_shells_sph, eri2_diagonal);
        }
    }

    return eri2_diagonal;
}