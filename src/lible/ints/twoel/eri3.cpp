#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/spherical_trafo.hpp>

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
                    vector<double> &fnx, vec3d &rints_tmp)
    {
        int labc = lab + lc;

        int dim_a = sp_data_ab.cdepths[2 * ipair_ab];
        int dim_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
        int dim_c = sh_data_c.cdepths[ishell_c];
        int pos_a = sp_data_ab.coffsets[2 * ipair_ab];
        int pos_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
        int pos_c = sh_data_c.coffsets[ishell_c];

        arma::vec::fixed<3> xyz_a{sp_data_ab.coords[6 * ipair_ab],
                                  sp_data_ab.coords[6 * ipair_ab + 1],
                                  sp_data_ab.coords[6 * ipair_ab + 2]};

        arma::vec::fixed<3> xyz_b{sp_data_ab.coords[6 * ipair_ab + 3],
                                  sp_data_ab.coords[6 * ipair_ab + 4],
                                  sp_data_ab.coords[6 * ipair_ab + 5]};

        arma::vec::fixed<3> xyz_c{sh_data_c.coords[3 * ishell_c],
                                  sh_data_c.coords[3 * ishell_c + 1],
                                  sh_data_c.coords[3 * ishell_c + 2]};

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++)
            {
                for (int ic = 0; ic < dim_c; ic++)
                {
                }
            }

        for (int ia = 0, iab = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++)
            {
            }
    }
}

lible::vec3d LIT::calcERI3(const Structure &structure)
{
    assert(structure.getUseRI());

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
        calcECoeffsSpherical(la, lb, sp_datas[ipair], ecoeffs_ipair);        

        ecoeffs[ipair] = std::move(ecoeffs_ipair);
    }

    vector<ShellData> sh_datas;
    for (int l = 0; l <= l_max_aux; l++)
        sh_datas.emplace_back(constructShellDataAux(l, structure));

    vector<vector<double>> ecoeffs_aux(l_max_aux);
    for (int l = 0; l <= l_max_aux; l++)
    {
        int n_ecoeffs = dimSphericals(l) * dimHermiteGaussians(l) * sh_datas[l].n_primitives;

        vector<double> ecoeffs_l(n_ecoeffs, 0);
        calcECoeffsSpherical(l, sh_datas[l], ecoeffs_l);

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
            vec3d rints_tmp(labc + 1, 0);

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

    return eri3;
}