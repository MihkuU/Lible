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
    void kernelERI4Diagonal(const int lab, const int ipair_ab, const vector<double> &ecoeffs_ab,
                            const vector<double> &ecoeffs_ab_tsp, 
                            const vector<array<int, 3>> &idxs_tuv,
                            const BoysF &boys_f, const ShellPairData &sp_data_ab,
                            vec4d &rints_tmp, vector<double> &eri4_shells_sph, vector<double> &fnx,
                            vector<double> &rints)
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

        arma::vec::fixed<3> xyz_a{sp_data_ab.coords[6 * ipair_ab],
                                  sp_data_ab.coords[6 * ipair_ab + 1],
                                  sp_data_ab.coords[6 * ipair_ab + 2]};

        arma::vec::fixed<3> xyz_b{sp_data_ab.coords[6 * ipair_ab + 3],
                                  sp_data_ab.coords[6 * ipair_ab + 4],
                                  sp_data_ab.coords[6 * ipair_ab + 5]};

        int dim_a_sph = dimSphericals(sp_data_ab.la);
        int dim_b_sph = dimSphericals(sp_data_ab.lb);
        int dim_tuv_ab = dimHermiteGaussians(lab);
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
                        double alpha = p * q / (p + q);

                        arma::vec::fixed<3> xyz_p = (a * xyz_a + b * xyz_b) / p;
                        arma::vec::fixed<3> xyz_q = (c * xyz_a + d * xyz_b) / q;
                        arma::vec::fixed<3> xyz_pq = xyz_p - xyz_q;

                        double x = alpha * arma::dot(xyz_pq, xyz_pq);
                        boys_f.calcFnx(labab, x, fnx);

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                        calcRIntsDiagonal(lab, fac, alpha, xyz_pq, fnx, idxs_tuv, rints_tmp,
                                          rints);

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
}

lible::vec2d LIT::calcERI4Diagonal(const Structure &structure)
{
    palPrint(fmt::format("Lible::{:<40}", "SHARK ERI4-diagonal..."));

    auto start{std::chrono::steady_clock::now()};

    int l_max = structure.getMaxL();

    vector<pair<int, int>> l_pairs = returnLPairs(l_max);

    vector<ShellPairData> sp_datas;    
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];
        sp_datas.emplace_back(constructShellPairData(la, lb, structure));
    }

    vector<vector<double>> ecoeffs(l_pairs.size());
    vector<vector<double>> ecoeffs_tsp(l_pairs.size());
    for (size_t ipair = 0; ipair < l_pairs.size(); ipair++)
    {
        auto [la, lb] = l_pairs[ipair];

        int lab = la + lb;
        int n_ecoeffs_sph = dimSphericals(la) * dimSphericals(lb) * dimHermiteGaussians(lab) *
                            sp_datas[ipair].n_prim_pairs;

        vector<double> ecoeffs_ipair(n_ecoeffs_sph, 0);
        vector<double> ecoeffs_tsp_ipair(n_ecoeffs_sph, 0);
        calcECoeffsSpherical(la, lb, sp_datas[ipair], ecoeffs_ipair, ecoeffs_tsp_ipair);

        ecoeffs[ipair] = std::move(ecoeffs_ipair);
        ecoeffs_tsp[ipair] = std::move(ecoeffs_tsp_ipair);
    }

    size_t dim_ao = structure.getDimAO();
    vec2d eri4_diagonal(dim_ao, dim_ao, 0);
    for (int lalb = 0; lalb < (int)l_pairs.size(); lalb++)
    {
        const auto &sp_data_ab = sp_datas[lalb];

        auto [la, lb] = l_pairs[lalb];

        int lab = la + lb;

        int dim_a_sph = dimSphericals(la);
        int dim_b_sph = dimSphericals(lb);
        int dim_ab_sph = dim_a_sph * dim_b_sph;
        int dim_tuv_ab = dimHermiteGaussians(lab);

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
                               sp_data_ab, rints_tmp, eri4_shells_sph, fnx, rints);

            transferIntegrals(ipair_ab, sp_data_ab, eri4_shells_sph, eri4_diagonal);
        }
    }

    auto end{std::chrono::steady_clock::now()};
    std::chrono::duration<double> duration{end - start};
    palPrint(fmt::format(" {:.2e} s\n", duration.count()));

    return eri4_diagonal;
}