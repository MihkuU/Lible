#include <hip/hip_runtime.h>
#include "structure_gpu.h"

using namespace Lible;
using namespace Lible::Ints;

using std::vector;

#define HIPCHECK(cmd)                                                                                           \
    {                                                                                                           \
        hipError_t error = cmd;                                                                                 \
        if (error != hipSuccess)                                                                                \
        {                                                                                                       \
            fprintf(stderr, "error: '%s'(%d) at %s:%d\n", hipGetErrorString(error), error, __FILE__, __LINE__); \
            exit(EXIT_FAILURE);                                                                                 \
        }                                                                                                       \
    }

StructureGPU::~StructureGPU()
{
    HIPCHECK(hipFree(shells_data.contr_exps));
    HIPCHECK(hipFree(shells_data.contr_coeffs));
    HIPCHECK(hipFree(shells_data.norm_coeffs));
    HIPCHECK(hipFree(shells_data.cart_exps));

    HIPCHECK(hipFree(shells_info.angmom));
    HIPCHECK(hipFree(shells_info.atomic_nr));
    HIPCHECK(hipFree(shells_info.dim_cart_exps));
    HIPCHECK(hipFree(shells_info.dim_contraction));
    HIPCHECK(hipFree(shells_info.dim_normalization));
    HIPCHECK(hipFree(shells_info.offset_cartexps));
    HIPCHECK(hipFree(shells_info.offset_contraction));
    HIPCHECK(hipFree(shells_info.offset_normalization));
    HIPCHECK(hipFree(shells_info.n_cartesian));
    HIPCHECK(hipFree(shells_info.n_spherical));
    HIPCHECK(hipFree(shells_info.pos_ao));
}

StructureGPU::StructureGPU(const std::unique_ptr<Structure> &structure) // Reconsider, maybe take only shells?
{
    shell_counts.resize(structure->shells.size());
    shell_offsets.resize(structure->shells.size());

    vector<double> contr_exps_flat, contr_coeffs_flat;
    vector<double> norm_coeffs_flat;
    vector<int> cart_exps_flat;

    vector<int> angmoms, atomic_nrs;
    vector<int> n_cartesians, n_sphericals;
    vector<int> pos_aos;
    vector<int> dim_cart_exps, offset_cart_exps;
    vector<int> dim_contraction, offset_contraction;
    vector<int> dim_normalization, offset_normalization;

    int ipos = 0;
    for (const auto &[angmom, shells_per_angmom] : structure->shells)
    {
        shell_counts[ipos] = shells_per_angmom.size();
        if (ipos == 0)
            shell_offsets[ipos] = 0;
        else
            shell_offsets[ipos] = shell_offsets[ipos - 1] + shell_counts[ipos - 1];
        ipos++;

        for (const auto &shell : shells_per_angmom)
        {
            contr_exps_flat.insert(contr_exps_flat.end(),
                                   shell.contraction_exps.begin(),
                                   shell.contraction_exps.end());
            contr_coeffs_flat.insert(contr_coeffs_flat.end(),
                                     shell.contraction_coeffs.begin(),
                                     shell.contraction_coeffs.end());
            norm_coeffs_flat.insert(norm_coeffs_flat.end(),
                                    shell.normalization.begin(),
                                    shell.normalization.end());
            for (const auto &trio : shell.cartesian_exps)
                cart_exps_flat.insert(cart_exps_flat.end(), trio.begin(), trio.end());

            angmoms.push_back(shell.angular_momentum);
            atomic_nrs.push_back(shell.atomic_number);
            n_cartesians.push_back(shell.dim_cartesian);
            n_sphericals.push_back(shell.dim_spherical);
            pos_aos.push_back(shell.pos);

            if (offset_cart_exps.size() == 0)
                offset_cart_exps.push_back(0);
            else
                offset_cart_exps.push_back(offset_cart_exps.back() + dim_cart_exps.back());
            dim_cart_exps.push_back(3 * shell.cartesian_exps.size());

            if (offset_contraction.size() == 0)
                offset_contraction.push_back(0);
            else
                offset_contraction.push_back(offset_contraction.back() + dim_contraction.back());
            dim_contraction.push_back(shell.contraction_exps.size());

            if (offset_normalization.size() == 0)
                offset_normalization.push_back(0);
            else
                offset_normalization.push_back(offset_normalization.back() + dim_normalization.back());
            dim_normalization.push_back(shell.normalization.size());
        }
    }

    HIPCHECK(hipMalloc(&shells_data.contr_exps, contr_exps_flat.size() * sizeof(double)));
    HIPCHECK(hipMalloc(&shells_data.contr_coeffs, contr_coeffs_flat.size() * sizeof(double)));
    HIPCHECK(hipMalloc(&shells_data.norm_coeffs, norm_coeffs_flat.size() * sizeof(double)));
    HIPCHECK(hipMalloc(&shells_data.cart_exps, cart_exps_flat.size() * sizeof(int)));

    HIPCHECK(hipMalloc(&shells_info.angmom, angmoms.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.atomic_nr, atomic_nrs.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.n_cartesian, n_cartesians.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.n_spherical, n_sphericals.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.pos_ao, pos_aos.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.dim_cart_exps, dim_cart_exps.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.dim_contraction, dim_contraction.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.dim_normalization, dim_normalization.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.offset_cartexps, offset_cart_exps.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.offset_contraction, offset_contraction.size() * sizeof(int)));
    HIPCHECK(hipMalloc(&shells_info.offset_normalization, offset_contraction.size() * sizeof(int)));

    HIPCHECK(hipMemcpy(shells_data.contr_exps, contr_exps_flat.data(),
                       contr_exps_flat.size() * sizeof(double), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_data.contr_coeffs, contr_coeffs_flat.data(),
                       contr_coeffs_flat.size() * sizeof(double), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_data.norm_coeffs, norm_coeffs_flat.data(),
                       norm_coeffs_flat.size() * sizeof(double), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_data.cart_exps, cart_exps_flat.data(),
                       cart_exps_flat.size() * sizeof(int), hipMemcpyHostToDevice));

    HIPCHECK(hipMemcpy(shells_info.angmom, angmoms.data(),
                       angmoms.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.atomic_nr, atomic_nrs.data(),
                       atomic_nrs.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.n_cartesian, n_cartesians.data(),
                       n_cartesians.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.n_spherical, n_sphericals.data(),
                       n_sphericals.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.pos_ao, pos_aos.data(),
                       pos_aos.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.dim_cart_exps, dim_cart_exps.data(),
                       dim_cart_exps.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.dim_contraction, dim_contraction.data(),
                       dim_contraction.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.dim_normalization, dim_normalization.data(),
                       dim_normalization.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.offset_cartexps, offset_cart_exps.data(),
                       offset_cart_exps.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.offset_contraction, offset_contraction.data(),
                       offset_contraction.size() * sizeof(int), hipMemcpyHostToDevice));
    HIPCHECK(hipMemcpy(shells_info.offset_normalization, offset_normalization.data(),
                       offset_normalization.size() * sizeof(int), hipMemcpyHostToDevice));
}

double *StructureGPU::allocateOneElInts(const size_t &n_ao)
{
    vector<double> ints_init(n_ao * n_ao, 0);
    double *ints_gpu;

    HIPCHECK(hipMalloc(&ints_gpu, n_ao * n_ao * sizeof(double)));
    HIPCHECK(hipMemcpy(ints_gpu, ints_init.data(), n_ao * n_ao * sizeof(double),
                       hipMemcpyHostToDevice));

    return ints_gpu;
}

void StructureGPU::deallocateOneElInts(double *one_el_ints)
{
    HIPCHECK(hipFree(one_el_ints));
}

void StructureGPU::copyOneElInts(const size_t &n_ao, double *one_el_ints_gpu,
                                 double *one_el_ints_cpu)
{
    HIPCHECK(hipMemcpy(one_el_ints_cpu, one_el_ints_gpu, n_ao * sizeof(double),
                       hipMemcpyDeviceToHost));
}