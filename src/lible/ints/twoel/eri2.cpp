#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/spherical_trafo.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

namespace lints = lible::ints;

using std::array, std::vector;

namespace lible::ints
{
    vector<double> eri2Diagonal(const Structure &structure);

    vec2d eri2(const Structure &structure);

    void transferIntsERI2Diag(const int ishell, const ShellData &sh_data, const vec2d &eri2_batch,
                              vector<double> &eri2_diagonal)
    {
        const int dim_a = numSphericals(sh_data.l);
        const int ofs_a = sh_data.offsets_sph[ishell];

        for (int a = 0; a < dim_a; a++)
        {
            const int mu = ofs_a + a;
            eri2_diagonal[mu] = eri2_batch(a, a);
        }
    }
}

lible::vec2d lints::eri2(const Structure &structure)
{
    if (structure.getUseRI() == false)
        throw std::runtime_error("eri2(): RI approximation is not enabled");

    const vector<ShellData> sh_datas = shellDataAux(structure);

    const int dim_ao_aux = structure.getDimAOAux();
    vec2d eri2(Fill(0), dim_ao_aux, dim_ao_aux);
    for (size_t ishdata_a = 0; ishdata_a < sh_datas.size(); ishdata_a++)
        for (size_t ishdata_b = 0; ishdata_b <= ishdata_a; ishdata_b++)
        {
            const auto &sh_data_a = sh_datas[ishdata_a];
            const auto &sh_data_b = sh_datas[ishdata_b];

            const int la = sh_data_a.l;
            const int lb = sh_data_b.l;
            const int n_sph_a = numSphericals(la);
            const int n_sph_b = numSphericals(lb);

            ERI2Kernel eri2_kernel = deployERI2Kernel(sh_data_a, sh_data_b);

#pragma omp parallel for
            for (int ishell_a = 0; ishell_a < sh_data_a.n_shells; ishell_a++)
            {
                const int bound_b = (la == lb) ? ishell_a + 1 : sh_data_b.n_shells;
                for (int ishell_b = 0; ishell_b < bound_b; ishell_b++)
                {
                    vec2d eri2_batch = eri2_kernel(ishell_a, ishell_b, sh_data_a, sh_data_b);

                    const int ofs_a = sh_data_a.offsets_sph[ishell_a];
                    const int ofs_b = sh_data_b.offsets_sph[ishell_b];
                    for (int ia = 0; ia < n_sph_a; ia++)
                        for (int ib = 0; ib < n_sph_b; ib++)
                        {
                            const int mu = ofs_a + ia;
                            const int nu = ofs_b + ib;
                            eri2(mu, nu) = eri2_batch(ia, ib);
                            eri2(nu, mu) = eri2_batch(ia, ib);
                        }
                }
            }
        }

    return eri2;
}

vector<double> lints::eri2Diagonal(const Structure &structure)
{
    if (structure.getUseRI() == false)
        throw std::runtime_error("eri2Diagonal(): RI approximation is not enabled");

    const vector<ShellData> sh_datas = shellDataAux(structure);

    vector<double> eri2_diagonal(structure.getDimAOAux(), 0);
    for (const ShellData &sh_data_a: sh_datas)
    {
        ERI2Kernel eri2_kernel = deployERI2Kernel(sh_data_a, sh_data_a);

#pragma omp parallel for
        for (int ishell = 0; ishell < sh_data_a.n_shells; ishell++)
        {
            vec2d eri2_batch = eri2_kernel(ishell, ishell, sh_data_a, sh_data_a);

            transferIntsERI2Diag(ishell, sh_data_a, eri2_batch, eri2_diagonal);
        }
    }

    return eri2_diagonal;
}