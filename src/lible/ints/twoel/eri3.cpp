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
    if (structure.getUseRI() == false)
        throw std::runtime_error("RI approximation is not enabled!");

    vector<ShellData> sh_datas = shellDataAux(structure);
    vector<ShellPairData> sp_data = shellPairData(true, structure);

    int dim_ao = structure.getDimAO();
    int dim_ao_aux = structure.getDimAOAux();
    vec3d eri3(Fill(0), dim_ao, dim_ao, dim_ao_aux);
    for (size_t ispdata_ab = 0; ispdata_ab < sp_data.size(); ispdata_ab++)
        for (size_t ishdata_c = 0; ishdata_c < sh_datas.size(); ishdata_c++)
        {
            const auto &sp_data_ab = sp_data[ispdata_ab];
            const auto &sh_data_c = sh_datas[ishdata_c];

            ERI3Kernel eri3_kernel = deployERI3Kernel(sp_data_ab, sh_data_c);

            for (int ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs; ipair_ab++)
                for (int ishell_c = 0; ishell_c < sh_data_c.n_shells; ishell_c++)
                {
                    vec3d eri3_batch = eri3_kernel(ipair_ab, ishell_c, sp_data_ab,
                                                   sh_data_c);

                    int ofs_a = sp_data_ab.offsets_sph[2 * ipair_ab];
                    int ofs_b = sp_data_ab.offsets_sph[2 * ipair_ab + 1];
                    int ofs_c = sh_data_c.offsets_sph[ishell_c];
                    for (size_t ia = 0; ia < eri3_batch.dim<0>(); ia++)
                        for (size_t ib = 0; ib < eri3_batch.dim<1>(); ib++)
                            for (size_t ic = 0; ic < eri3_batch.dim<2>(); ic++)
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