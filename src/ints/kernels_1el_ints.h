#pragma once

#include <functional>

#include "ints.h"
#include "shell.h"
#include "experimental/structure_gpu.h"

namespace lible
{
    // namespace Kernels1El
    // {
    //     template <ints::Option1El option>
    //     void oneElIntKernel(const int &la, const int &lb, const std::size_t &n_ao,
    //                         const std::vector<Shells::ShellPair> &shell_pairs,
    //                         std::vector<double> &one_el_ints);
    //     // TODO: change this to smth like:
    //     // void oneElIntKernel(const int &la, const int &lb,        
    //     //                     const std::unique_ptr<ints::Structure> &structure,
    //     //                     std::vector<double> &one_el_ints);

    //     template <ints::Option1El option>
    //     void oneElIntKernelGPU(const int &la, const int &lb,
    //                            const std::unique_ptr<ints::StructureGPU> &structure_gpu,
    //                            double *one_el_ints_gpu);

    //     namespace ObaraSaika
    //     {
    //         // TODO: move this somewhere else prolly
    //         void overlapX(const int &angmom_a, const int &angmom_b, const double &one_over_2p,
    //                       const double &x_pa, const double &x_pb, const double &overlap_00,
    //                       std::vector<double> &overlap_x);
    //     }
    // }
}