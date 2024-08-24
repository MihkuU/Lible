#pragma once

#include <hip/hip_runtime.h>

#define hipCheck(call)                                                                              \
    do                                                                                              \
    {                                                                                               \
        hipError_t gpuErr = call;                                                                   \
        if (hipSuccess != gpuErr)                                                                   \
        {                                                                                           \
            printf("GPU API Error - %s:%d: '%s'\n", __FILE__, __LINE__, hipGetErrorString(gpuErr)); \
            exit(1);                                                                                \
        }                                                                                           \
    } while (0)

namespace lible
{
    namespace ints
    {
        // namespace gpu
        // {
        //     /** */
        //     // extern "C" 
        //     __device__ int idxE(int d3, int d23, int i, int j, int t, int pos);

        //     /** */
        //     // extern "C"
        //     __device__ int dimCartIdxs(int l);
        // }
    }
}