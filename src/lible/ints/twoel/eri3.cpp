#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

namespace LI = lible::ints;

using std::array, std::pair, std::vector;

namespace lible::ints
{
    vec3d eri3(const Structure &structure);
}

lible::vec3d LI::eri3(const Structure &structure)
{
    if (!structure.getUseRI())
        throw std::runtime_error("RI approximation is not enabled!");

    int l_max_aux = structure.getMaxLAux();
    vector<pair<int, int>> l_pairs = getLPairsSymm(structure.getMaxL());

    vector<ShellData> sh_datas = shellDatasAux(l_max_aux, structure);
    vector<ShellPairData> sp_datas = shellPairDatasSymm(l_pairs, structure);

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

            ERI3Kernel eri3_kernel = deployERI3Kernel(sp_data_ab, sh_data_c);

            for (int ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs; ipair_ab++)
                for (int ishell_c = 0; ishell_c < sh_data_c.n_shells; ishell_c++)
                {
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