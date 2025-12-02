#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

namespace lints = lible::ints;

namespace lible::ints
{
    /// Copies the integrals from a shell batch to the target container.
    void transferIntsERI2Diag(size_t ishell, const ShellData &sh_data, const vec2d &eri2_batch,
                              std::vector<double> &eri2_diagonal);
}

void lints::transferIntsERI2Diag(const size_t ishell, const ShellData &sh_data, const vec2d &eri2_batch,
                                 std::vector<double> &eri2_diagonal)
{
    int dim_a = numSphericals(sh_data.l_);
    size_t ofs_a = sh_data.offsets_sph_[ishell];

    for (int a = 0; a < dim_a; a++)
    {
        size_t mu = ofs_a + a;
        eri2_diagonal[mu] = eri2_batch(a, a);
    }
}

lible::vec2d lints::eri2(const Structure &structure)
{
    if (structure.getUseRI() == false)
        throw std::runtime_error("eri2(): RI approximation is not enabled");

    std::vector<ShellData> sh_datas = shellDataAux(structure);

    size_t dim_ao_aux = structure.getDimAOAux();
    vec2d eri2(Fill(0), dim_ao_aux, dim_ao_aux);
    for (size_t ishdata_a = 0; ishdata_a < sh_datas.size(); ishdata_a++)
        for (size_t ishdata_b = 0; ishdata_b <= ishdata_a; ishdata_b++)
        {
            ShellData &sh_data_a = sh_datas[ishdata_a];
            ShellData &sh_data_b = sh_datas[ishdata_b];

            int la = sh_data_a.l_;
            int lb = sh_data_b.l_;
            int n_sph_a = numSphericals(la);
            int n_sph_b = numSphericals(lb);

            ERI2Kernel eri2_kernel(sh_data_a, sh_data_b);

#pragma omp parallel for
            for (size_t ishell_a = 0; ishell_a < sh_data_a.n_shells_; ishell_a++)
            {
                size_t bound_b = (la == lb) ? ishell_a + 1 : sh_data_b.n_shells_;
                for (size_t ishell_b = 0; ishell_b < bound_b; ishell_b++)
                {
                    vec2d eri2_batch = eri2_kernel(ishell_a, ishell_b, sh_data_a, sh_data_b);

                    size_t ofs_a = sh_data_a.offsets_sph_[ishell_a];
                    size_t ofs_b = sh_data_b.offsets_sph_[ishell_b];
                    for (int ia = 0; ia < n_sph_a; ia++)
                        for (int ib = 0; ib < n_sph_b; ib++)
                        {
                            size_t mu = ofs_a + ia;
                            size_t nu = ofs_b + ib;
                            eri2(mu, nu) = eri2_batch(ia, ib);
                            eri2(nu, mu) = eri2_batch(ia, ib);
                        }
                }
            }
        }

    return eri2;
}

std::vector<double> lints::eri2Diagonal(const Structure &structure)
{
    if (structure.getUseRI() == false)
        throw std::runtime_error("eri2Diagonal(): RI approximation is not enabled");

    std::vector<ShellData> sh_datas = shellDataAux(structure);

    std::vector<double> eri2_diagonal(structure.getDimAOAux(), 0);
    for (const ShellData &sh_data_a : sh_datas)
    {
        ERI2Kernel eri2_kernel(sh_data_a, sh_data_a);

#pragma omp parallel for
        for (size_t ishell = 0; ishell < sh_data_a.n_shells_; ishell++)
        {
            vec2d eri2_batch = eri2_kernel(ishell, ishell, sh_data_a, sh_data_a);

            transferIntsERI2Diag(ishell, sh_data_a, eri2_batch, eri2_diagonal);
        }
     }

     return eri2_diagonal;
}
