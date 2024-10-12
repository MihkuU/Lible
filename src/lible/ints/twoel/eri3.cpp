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

using std::array, std::pair, std::vector;

namespace lible::ints::two
{
    void kernelERI3(const int lab, const int lc, const int ipair_ab, const int ishell_c,
                    const vector<double> &ecoeffs_ab, const vector<double> &ecoeffs_c,
                    const vector<array<int, 3>> &idxs_tuv_ab,
                    const vector<array<int, 3>> &idxs_tuv_c,
                    const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                    const BoysF &boys_f, vector<double> &eri3_shells_sph, vector<double> &rints,
                    vector<double> &fnx, vec4d &rints_tmp)
    {
        int labc = lab + lc;

        int dim_a = sp_data_ab.cdepths[2 * ipair_ab];
        int dim_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
        int dim_c = sh_data_c.cdepths[ishell_c];
        int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
        int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
        int pos_c = sh_data_c.coffsets[ishell_c];

        array<double, 3> xyz_a{sp_data_ab.coords[6 * ipair_ab],
                               sp_data_ab.coords[6 * ipair_ab + 1],
                               sp_data_ab.coords[6 * ipair_ab + 2]};

        array<double, 3> xyz_b{sp_data_ab.coords[6 * ipair_ab + 3],
                               sp_data_ab.coords[6 * ipair_ab + 4],
                               sp_data_ab.coords[6 * ipair_ab + 5]};

        array<double, 3> xyz_c{sh_data_c.coords[3 * ishell_c],
                               sh_data_c.coords[3 * ishell_c + 1],
                               sh_data_c.coords[3 * ishell_c + 2]};

        int dim_sph_a = dimSphericals(sp_data_ab.la);
        int dim_sph_b = dimSphericals(sp_data_ab.lb);
        int dim_sph_c = dimSphericals(sh_data_c.l);
        int dim_sph_ab = dim_sph_a * dim_sph_b;
        int dim_tuv_ab = dimHermiteGaussians(lab);
        int dim_tuv_c = dimHermiteGaussians(lc);
        int dim_ecoeffs_ab = dim_sph_ab * dim_tuv_ab;
        int dim_ecoeffs_c = dim_sph_c * dim_tuv_c;
        int dim_rints_x_ecoeffs = dim_tuv_ab * dim_sph_c;
        vector<double> rints_x_ecoeffs(dim_a * dim_b * dim_rints_x_ecoeffs, 0);

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                int pos_rints_x_ecoeffs = iab * dim_rints_x_ecoeffs;
                for (int ic = 0; ic < dim_c; ic++)
                {
                    double a = sp_data_ab.exps[pos_a + ia];
                    double b = sp_data_ab.exps[pos_b + ib];
                    double c = sh_data_c.exps[pos_c + ic];

                    double p = a + b;
                    double alpha = p * c / (p + c);

                    array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                           (a * xyz_a[1] + b * xyz_b[1]) / p,
                                           (a * xyz_a[2] + b * xyz_b[2]) / p};

                    array<double, 3> xyz_pc{xyz_p[0] - xyz_c[0], xyz_p[1] - xyz_c[1],
                                            xyz_p[2] - xyz_c[2]};

                    double xx{xyz_pc[0]}, xy{xyz_pc[1]}, xz{xyz_pc[2]};
                    double x = alpha * (xx * xx + xy * xy + xz * xz);
                    boys_f.calcFnx(labc, x, fnx);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                    calcRInts(lab, lc, fac, alpha, xyz_pc, fnx, idxs_tuv_ab, idxs_tuv_c,
                              rints_tmp, rints);

                    int pos_ecoeffs_c = sh_data_c.offsets_ecoeffs[ishell_c] + ic * dim_ecoeffs_c;

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, dim_tuv_ab,
                                dim_sph_c, dim_tuv_c, 1.0, &rints[0], dim_tuv_c,
                                &ecoeffs_c[pos_ecoeffs_c], dim_tuv_c, 1.0,
                                &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_c);
                }
            }

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++, iab++)
            {
                int pos_rints_x_ecoeffs = iab * dim_rints_x_ecoeffs;
                int pos_ecoeffs_ab = sp_data_ab.offsets_ecoeffs[ipair_ab] + iab * dim_ecoeffs_ab;

                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim_sph_ab, dim_sph_c,
                            dim_tuv_ab, 1.0, &ecoeffs_ab[pos_ecoeffs_ab], dim_tuv_ab,
                            &rints_x_ecoeffs[pos_rints_x_ecoeffs], dim_sph_c, 1.0,
                            &eri3_shells_sph[0], dim_sph_c);
            }
    }
}

lible::vec3d LIT::calcERI3(const Structure &structure)
{
    assert(structure.getUseRI()); // TODO: use throw exception here instead

    palPrint(fmt::format("Lible::{:<40}", "SHARK ERI3..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();
    int l_max_aux = structure.getMaxLAux();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> sp_datas;
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        sp_datas.emplace_back(constructShellPairData(la, lb, structure));
    }

    vector<vector<double>> ecoeffs(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        int lab = la + lb;
        int n_ecoeffs_sph = dimSphericals(la) * dimSphericals(lb) * dimHermiteGaussians(lab) *
                            sp_datas[ipair].n_prim_pairs;

        vector<double> ecoeffs_ipair(n_ecoeffs_sph, 0);
        ecoeffsSPsSpherical(la, lb, sp_datas[ipair], ecoeffs_ipair);

        ecoeffs[ipair] = std::move(ecoeffs_ipair);
    }

    vector<ShellData> sh_datas;
    for (int l = 0; l <= l_max_aux; l++)
        sh_datas.emplace_back(constructShellDataAux(l, structure));

    vector<vector<double>> ecoeffs_aux(sh_datas.size());
    for (int l = 0; l <= l_max_aux; l++)
    {
        int n_ecoeffs = dimSphericals(l) * dimHermiteGaussians(l) * sh_datas[l].n_primitives;

        vector<double> ecoeffs_l(n_ecoeffs, 0);
        ecoeffsShellsSpherical(l, sh_datas[l], ecoeffs_l);

        ecoeffs_aux[l] = std::move(ecoeffs_l);
    }

    size_t dim_ao = structure.getDimAO();
    size_t dim_ao_aux = structure.getDimAOAux();
    vec3d eri3(dim_ao, dim_ao, dim_ao_aux, 0);
    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
        for (int lc = 0; lc <= l_max_aux; lc++)
        {
            const auto &sp_data_ab = sp_datas[lalb];
            const auto &sh_data_c = sh_datas[lc];

            auto [la, lb] = l_pairs[lalb];
            int lab = la + lb;

            int dim_a_sph = dimSphericals(la);
            int dim_b_sph = dimSphericals(lb);
            int dim_c_sph = dimSphericals(lc);
            int dim_ab_sph = dim_a_sph * dim_b_sph;
            int dim_tuv_ab = dimHermiteGaussians(lab);
            int dim_tuv_c = dimHermiteGaussians(lc);

            int n_pairs_ab = sp_data_ab.n_pairs;
            int n_shells_c = sh_data_c.n_shells;

            int labc = lab + lc;
            BoysF boys_f(labc);

            const vector<double> &ecoeffs_ab = ecoeffs[lalb];
            const vector<double> &ecoeffs_c = ecoeffs_aux[lc];

            vector<array<int, 3>> idxs_tuv_ab = returnHermiteGaussianIdxs(lab);
            vector<array<int, 3>> idxs_tuv_c = returnHermiteGaussianIdxs(lc);

            vector<double> rints(dim_tuv_ab * dim_tuv_c, 0);
            vector<double> fnx(labc + 1, 0);
            vec4d rints_tmp(labc + 1, 0);

            for (int ipair_ab = 0; ipair_ab < n_pairs_ab; ipair_ab++)
                for (int ishell_c = 0; ishell_c < n_shells_c; ishell_c++)
                {
                    vector<double> eri3_shells_sph(dim_ab_sph * dim_c_sph, 0);

                    kernelERI3(lab, lc, ipair_ab, ishell_c, ecoeffs_ab, ecoeffs_c, idxs_tuv_ab,
                               idxs_tuv_c, sp_data_ab, sh_data_c, boys_f, eri3_shells_sph, rints,
                               fnx, rints_tmp);

                    transferIntegrals(ipair_ab, ishell_c, sh_data_c, sp_data_ab, eri3_shells_sph,
                                      eri3);
                }
        }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));

    return eri3;
}