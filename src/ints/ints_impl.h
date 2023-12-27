#pragma once

#include "kernels_1el_ints.h"
#include "spherical_trafo.h"
#include "structure.h"
#include "experimental/structure_gpu.h"
#include "ints.h"

namespace Lible
{
    // template <Ints::Option1El option>
    // std::vector<double> Ints::calcOneElInts()
    // {
    //     int max_l = structure->max_angular_momentum;
    //     std::size_t n_ao = structure->n_atomic_orbitals;

    //     std::vector<double> one_el_ints(n_ao * n_ao, 0);
    //     for (int la = max_l; la >= 0; la--)
    //         for (int lb = la; lb >= 0; lb--)
    //             Kernels1El::oneElIntKernel<option>(la, lb, n_ao,
    //                                                structure->shell_pairs.at(std::make_pair(la, lb)),
    //                                                one_el_ints);

    //     return one_el_ints;
    // }

    // template <Ints::Option1El option>
    // std::vector<double> Ints::calcOneElIntsGPU()
    // {
    //     int max_l = structure->max_angular_momentum;
    //     std::size_t n_ao = structure->n_atomic_orbitals;

    //     double *one_el_ints_gpu = structure_gpu->allocateOneElInts(n_ao);

    //     for (int la = max_l; la >= 0; la--)
    //         for (int lb = la; lb >= 0; lb--)
    //             Kernels1El::oneElIntKernelGPU<option>(la, lb, structure_gpu, one_el_ints_gpu);

    //     std::vector<double> one_el_ints(n_ao * n_ao);
    //     structure_gpu->copyOneElInts(n_ao, one_el_ints_gpu, one_el_ints.data());
    //     structure_gpu->deallocateOneElInts(one_el_ints_gpu);

    //     return one_el_ints;
    // }
}