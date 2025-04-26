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