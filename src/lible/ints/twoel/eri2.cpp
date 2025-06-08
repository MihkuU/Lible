#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LI = lible::ints;

using std::array, std::vector;

namespace lible::ints
{
    vector<double> eri2Diagonal(const Structure &structure);

    vec2d eri2(const Structure &structure);    

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

    void transferIntsERI2Diag(const int ishell, const ShellData &sh_data,
                              const vector<double> &eri2_shells_sph,
                              vector<double> &eri2_diagonal)
    {
        int dim_a = numSphericals(sh_data.l);
        int pos_a = sh_data.offsets_sph[ishell];
        int pos_norm_a = sh_data.offsets_norms[ishell];

        for (int mu = 0; mu < dim_a; mu++)
        {
            int munu = mu * dim_a + mu;

            double norm_a = sh_data.norms[pos_norm_a];
            double normalized_int = norm_a * norm_a * eri2_shells_sph[munu];

            int a = pos_a + mu;

            eri2_diagonal[a] = normalized_int;
        }
    }
}

lible::vec2d LI::eri2(const Structure &structure)
{
    if (!structure.getUseRI())
        throw std::runtime_error("RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();
    vector<ShellData> sh_datas = shellDatasAux(l_max_aux, structure);

    size_t dim_ao_aux = structure.getDimAOAux();
    vec2d eri2(Fill(0), dim_ao_aux, dim_ao_aux);
    for (int la = 0; la <= l_max_aux; la++)
        for (int lb = 0; lb <= la; lb++)
        {
            const auto &sh_data_a = sh_datas[la];
            const auto &sh_data_b = sh_datas[lb];

            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);

            ERI2Kernel eri2_kernel = deployERI2Kernel(sh_data_a, sh_data_b);

            for (int ishell_a = 0; ishell_a < sh_data_a.n_shells; ishell_a++)
            {
                int bound_b = (la == lb) ? ishell_a + 1 : sh_data_b.n_shells;
                for (int ishell_b = 0; ishell_b < bound_b; ishell_b++)
                {
                    vec2d eri2_batch = eri2_kernel(ishell_a, ishell_b, sh_data_a, sh_data_b);

                    int ofs_a = sh_data_a.offsets_sph[ishell_a];
                    int ofs_b = sh_data_b.offsets_sph[ishell_b];

                    for (int ia = 0; ia < n_sph_a; ia++)
                        for (int ib = 0; ib < n_sph_b; ib++)
                        {
                            int mu = ofs_a + ia;
                            int nu = ofs_b + ib;
                            eri2(mu, nu) = eri2_batch(ia, ib);
                            eri2(nu, mu) = eri2_batch(ia, ib);
                        }
                }
            }
        }

    return eri2;
}

vector<double> LI::eri2Diagonal(const Structure &structure)
{
    if (!structure.getUseRI())
        throw std::runtime_error("RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();
    vector<ShellData> sh_datas = shellDatasAux(l_max_aux, structure);

    size_t dim_ao_aux = structure.getDimAOAux();
    vector<double> eri2_diagonal(dim_ao_aux, 0);
    for (int la = 0; la <= l_max_aux; la++)
    {
        const auto &sh_data_a = sh_datas[la];

        int n_shells_a = sh_data_a.n_shells;

        int dim_sph_a = numSphericals(la);

        int laa = la + la;
        BoysF boys_f(laa);

        auto [ecoeffs_a, ecoeffs_a_tsp] = ecoeffsSphericalShellData_BraKet(sh_datas[la]);

        vector<array<int, 3>> idxs_tuv_a = getHermiteGaussianIdxs(la);

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