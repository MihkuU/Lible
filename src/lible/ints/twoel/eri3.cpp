#include <lible/ints/ints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

namespace lints = lible::ints;

lible::vec3d lints::eri3(const Structure &structure)
{
    if (structure.getUseRI() == false)
        throw std::runtime_error("RI approximation is not enabled!");

    std::vector<ShellData> sh_datas = shellDataAux(structure);
    std::vector<ShellPairData> sp_data = shellPairData(true, structure);

    size_t dim_ao = structure.getDimAO();
    size_t dim_ao_aux = structure.getDimAOAux();
    vec3d eri3(Fill(0), dim_ao, dim_ao, dim_ao_aux);
    for (const auto &sp_data_ab : sp_data)
        for (const auto &sh_data_c : sh_datas)
        {
            ERI3Kernel eri3_kernel(sp_data_ab, sh_data_c);

#pragma omp parallel for
            for (size_t ipair_ab = 0; ipair_ab < sp_data_ab.n_pairs_; ipair_ab++)
                for (size_t ishell_c = 0; ishell_c < sh_data_c.n_shells_; ishell_c++)
                {
                    vec3d eri3_batch = eri3_kernel(ipair_ab, ishell_c, sp_data_ab,
                                                   sh_data_c);

                    size_t ofs_a = sp_data_ab.offsets_sph_[2 * ipair_ab];
                    size_t ofs_b = sp_data_ab.offsets_sph_[2 * ipair_ab + 1];
                    size_t ofs_c = sh_data_c.offsets_sph_[ishell_c];
                    for (size_t ia = 0; ia < eri3_batch.dim<0>(); ia++)
                        for (size_t ib = 0; ib < eri3_batch.dim<1>(); ib++)
                            for (size_t ic = 0; ic < eri3_batch.dim<2>(); ic++)
                            {
                                size_t mu = ofs_a + ia;
                                size_t nu = ofs_b + ib;
                                size_t ka = ofs_c + ic;

                                eri3(mu, nu, ka) = eri3_batch(ia, ib, ic);
                                eri3(nu, mu, ka) = eri3_batch(ia, ib, ic);
                            }
                }
        }

    return eri3;
}