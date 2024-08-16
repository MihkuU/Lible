#pragma once

#include <hip/hip_runtime.h>

namespace lible::ints::gpu
{

    __device__ inline int idxE(int d3, int d23, int i, int j, int t, int pos)
    {
        return pos + i * d23 + j * d3 + t;
    }

    __device__ inline int dimCartIdxs(int l)
    {
        return (l + 1) * (l + 1) / 2;
    }
}