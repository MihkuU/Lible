#include <hip/hip_runtime.h>
#include <stdio.h>

#include "../kernels_1el_ints.h"

using namespace Lible;
using std::vector;

__device__ void overlapGPUInner()
{
}

__global__ void overlapGPU(int count_shells_a, int count_shells_b,
                           int offset_shells_a, int offset_shells_b,
                           double *overlap_ints)
{
    for (int ishella = 0; ishella < count_shells_a; ishella++)
    {
        for (int ishellb = ishella; ishellb < count_shells_b; ishellb++)
        {
        }
    }
}

__global__ void overlapGPU_sameL()
{
}

// template <>
// void Kernels1El::oneElIntKernelGPU<Ints::OVERLAP>(const int &la, const int &lb,
//                                                   const std::unique_ptr<Ints::StructureGPU> &structure_gpu,
//                                                   double *one_el_ints_gpu)
// {
//     if (la != lb)
//     {
//         // overlapGPU<<<1,1>>>();
//     }
//     else
//     {
//         // overlapGPU_sameL();
//     }
//     // overlapGPU_sameL<<<1, 1>>>(structure_gpu->getShellCounts()[la], structure_gpu->getShellCounts()[lb],
//     //                      structure_gpu->getShellOffsets()[la], structure_gpu->getShellOffsets()[lb]);
// }