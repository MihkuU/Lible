#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/utils.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>

#include <fmt/core.h>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LIT = lible::ints::two;

using std::array, std::vector;

namespace lible::ints::two
{
    void kernelERI2(const int ishell_a, const int ishell_b, const int la, const int lb,
                    const vector<double> &ecoeffs_a, const vector<double> &ecoeffs_b_tsp,
                    const vector<array<int, 3>> &idxs_tuv_a,
                    const vector<array<int, 3>> &idxs_tuv_b,
                    const BoysF &boys_f, const ShellData &sh_data_a, const ShellData &sh_data_b,
                    vector<double> &fnx, vector<double> &rints, vec4d &rints_tmp,
                    vector<double> &eri2_shells_sph)
    {
        int lab = la + lb;
        int dim_a = sh_data_a.cdepths[ishell_a];
        int dim_b = sh_data_b.cdepths[ishell_b];
        int pos_a = sh_data_a.coffsets[ishell_a];
        int pos_b = sh_data_b.coffsets[ishell_b];

        array<double, 3> xyz_a{sh_data_a.coords[3 * ishell_a],
                               sh_data_a.coords[3 * ishell_a + 1],
                               sh_data_a.coords[3 * ishell_a + 2]};

        array<double, 3> xyz_b{sh_data_b.coords[3 * ishell_b],
                               sh_data_b.coords[3 * ishell_b + 1],
                               sh_data_b.coords[3 * ishell_b + 2]};

        array<double, 3> xyz_ab{xyz_a[0] - xyz_b[0], xyz_a[1] - xyz_b[1],
                                xyz_a[2] - xyz_b[2]};

        double xx{xyz_ab[0]}, xy{xyz_ab[1]}, xz{xyz_ab[2]};
        double xyz_ab_dot = xx * xx + xy * xy + xz * xz;

        int dim_sph_a = dimSphericals(la);
        int dim_sph_b = dimSphericals(lb);
        int dim_tuv_a = dimHermiteGaussians(la);
        int dim_tuv_b = dimHermiteGaussians(lb);
        int dim_ecoeffs_a = dim_sph_a * dim_tuv_a;
        int dim_ecoeffs_b = dim_sph_b * dim_tuv_b;
        int dim_rints_x_ecoeffs = dim_tuv_a * dim_sph_b;
        vector<double> rints_x_ecoeffs(dim_a * dim_rints_x_ecoeffs, 0);

        for (int ia = 0; ia < dim_a; ia++)
        {
            int pos_rints_x_ecoeffs = ia * dim_rints_x_ecoeffs;
            for (int ib = 0; ib < dim_b; ib++)
            {
                double a = sh_data_a.exps[pos_a + ia];
                double b = sh_data_b.exps[pos_b + ib];

                double alpha = a * b / (a + b);
                double x = alpha * xyz_ab_dot;
                boys_f.calcFnx(lab, x, fnx);

                double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
                calcRInts(la, lb, fac, alpha, xyz_ab, fnx, idxs_tuv_a, idxs_tuv_b, rints_tmp,
                          rints);

                int pos_ecoeffs_b = sh_data_b.offsets_ecoeffs[ishell_b] + ib * dim_ecoeffs_b;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_tuv_a, dim_sph_b,
                            dim_tuv_b, 1.0, &rints[0], dim_tuv_b, &ecoeffs_b_tsp[pos_ecoeffs_b],
                            dim_sph_b, 1.0, &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_b);
            }
        }

        for (int ia = 0; ia < dim_a; ia++)
        {
            int pos_rints_x_ecoeffs = ia * dim_rints_x_ecoeffs;
            int pos_ecoeffs_a = sh_data_a.offsets_ecoeffs[ishell_a] + ia * dim_ecoeffs_a;

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_a, dim_sph_b,
                        dim_tuv_a, 1.0, &ecoeffs_a[pos_ecoeffs_a], dim_tuv_a,
                        &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_b, 1.0, &eri2_shells_sph[0],
                        dim_sph_b);
        }
    }

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

        int dim_sph_a = dimSphericals(la);
        int dim_tuv_a = dimHermiteGaussians(la);
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
    assert(structure.getUseRI()); // TODO: use throw exception here instead

    palPrint(fmt::format("Lible::{:<40}", "SHARK ERI2..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max_aux = structure.getMaxLAux();

    vector<ShellData> sh_datas;
    for (int l = 0; l <= l_max_aux; l++)
        sh_datas.emplace_back(constructShellDataAux(l, structure));

    vector<vector<double>> ecoeffs(sh_datas.size());
    vector<vector<double>> ecoeffs_tsp(sh_datas.size());
    for (int l = 0; l <= l_max_aux; l++)
    {
        int n_ecoeffs = dimSphericals(l) * dimHermiteGaussians(l) * sh_datas[l].n_primitives;

        vector<double> ecoeffs_l(n_ecoeffs, 0);
        vector<double> ecoeffs_tsp_l(n_ecoeffs, 0);
        ecoeffsShellsSpherical(l, sh_datas[l], ecoeffs_l, ecoeffs_tsp_l);

        ecoeffs[l] = std::move(ecoeffs_l);
        ecoeffs_tsp[l] = std::move(ecoeffs_tsp_l);
    }

    size_t dim_ao_aux = structure.getDimAOAux();
    vec2d eri2(dim_ao_aux, dim_ao_aux, 0);
    for (int la = 0; la <= l_max_aux; la++)
        for (int lb = 0; lb <= la; lb++)
        {
            const auto &sh_data_a = sh_datas[la];
            const auto &sh_data_b = sh_datas[lb];

            int n_shells_a = sh_data_a.n_shells;
            int n_shells_b = sh_data_b.n_shells;
            int dim_sph_a = dimSphericals(la);
            int dim_sph_b = dimSphericals(lb);
            int dim_tuv_a = dimHermiteGaussians(la);
            int dim_tuv_b = dimHermiteGaussians(lb);

            int lab = la + lb;
            BoysF boys_f(lab);

            const vector<double> &ecoeffs_a = ecoeffs[la];
            const vector<double> &ecoeffs_b_tsp = ecoeffs_tsp[lb];

            vector<array<int, 3>> idxs_tuv_a = returnHermiteGaussianIdxs(la);
            vector<array<int, 3>> idxs_tuv_b = returnHermiteGaussianIdxs(lb);

            vector<double> rints(dim_tuv_a * dim_tuv_b, 0);
            vector<double> fnx(lab + 1, 0);
            vec4d rints_tmp(lab + 1, 0);

            if (la == lb)
            {
                for (int ishell_a = 0; ishell_a < n_shells_a; ishell_a++)
                    for (int ishell_b = 0; ishell_b <= ishell_a; ishell_b++)
                    {
                        vector<double> eri2_shells_sph(dim_sph_a * dim_sph_b, 0);

                        kernelERI2(ishell_a, ishell_b, la, lb, ecoeffs_a, ecoeffs_b_tsp,
                                   idxs_tuv_a, idxs_tuv_b, boys_f, sh_data_a, sh_data_b, fnx,
                                   rints, rints_tmp, eri2_shells_sph);

                        transferIntegrals(ishell_a, ishell_b, sh_data_a, sh_data_b,
                                          eri2_shells_sph, eri2);
                    }
            }
            else
            {
                for (int ishell_a = 0; ishell_a < n_shells_a; ishell_a++)
                    for (int ishell_b = 0; ishell_b < n_shells_b; ishell_b++)
                    {
                        vector<double> eri2_shells_sph(dim_sph_a * dim_sph_b, 0);

                        kernelERI2(ishell_a, ishell_b, la, lb, ecoeffs_a, ecoeffs_b_tsp,
                                   idxs_tuv_a, idxs_tuv_b, boys_f, sh_data_a, sh_data_b, fnx,
                                   rints, rints_tmp, eri2_shells_sph);

                        transferIntegrals(ishell_a, ishell_b, sh_data_a, sh_data_b,
                                          eri2_shells_sph, eri2);
                    }
            }
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));

    return eri2;
}

vector<double> LIT::calcERI2Diagonal(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}", "SHARK ERI2-diagonal..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max_aux = structure.getMaxLAux();

    vector<ShellData> sh_datas;
    for (int l = 0; l <= l_max_aux; l++)
        sh_datas.emplace_back(constructShellDataAux(l, structure));

    vector<vector<double>> ecoeffs(sh_datas.size());
    vector<vector<double>> ecoeffs_tsp(sh_datas.size());
    for (int l = 0; l <= l_max_aux; l++)
    {
        int n_ecoeffs = dimSphericals(l) * dimHermiteGaussians(l) * sh_datas[l].n_primitives;

        vector<double> ecoeffs_l(n_ecoeffs, 0);
        vector<double> ecoeffs_tsp_l(n_ecoeffs, 0);
        ecoeffsShellsSpherical(l, sh_datas[l], ecoeffs_l, ecoeffs_tsp_l);

        ecoeffs[l] = std::move(ecoeffs_l);
        ecoeffs_tsp[l] = std::move(ecoeffs_tsp_l);
    }

    size_t dim_ao_aux = structure.getDimAOAux();
    vector<double> eri2_diagonal(dim_ao_aux, 0);
    for (int la = 0; la <= l_max_aux; la++)
    {
        const auto &sh_data_a = sh_datas[la];

        int n_shells_a = sh_data_a.n_shells;

        int dim_sph_a = dimSphericals(la);
        int dim_tuv_a = dimHermiteGaussians(la);

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

            transferIntegrals(ishell, sh_data_a, eri2_shells_sph, eri2_diagonal);
        }
    }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));

    return eri2_diagonal;
}