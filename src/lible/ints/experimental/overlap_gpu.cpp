#include <lible/ints/experimental/oneel_detail_gpu.hpp>
#include <lible/util.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <fmt/core.h>
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

namespace LIOG = lible::ints::one_gpu;

using std::vector;

lible::vec2d LIOG::calculateS_L0(const Structure &structure)
{
    printf("LIOG::calculateS_L0\n");
    int l_max = structure.getMaxL();
    size_t dim_ao = structure.getDimAO();

    ShellPairData shell_pair_data = ShellPairData(0, 0, structure);

    vec2d overlap_ints(dim_ao, dim_ao, 0);

    // host data
    size_t n_pairs = shell_pair_data.n_pairs;
    size_t n_pairs_2x = 2 * n_pairs;
    size_t n_coords = 6 * n_pairs;

    size_t n_ccoeffs = 0, n_norms = 0;
    for (size_t ipair = 0; ipair < n_pairs; ipair++)
    {
        const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
        const auto &[norms_a, norms_b] = shell_pair_data.norms[ipair];
        n_ccoeffs += ccoeffs_a.size();
        n_ccoeffs += ccoeffs_b.size();
        n_norms += norms_a.size();
        n_norms += norms_b.size();
    }

    printf("n_pairs   = %zu\n", n_pairs);
    printf("n_ccoeffs = %zu\n", n_ccoeffs);

    // vector<double> ccoeffs(n_ccoeffs);
    // vector<double> exps(n_ccoeffs);
    // vector<double> coords(6 * n_pairs);

    vector<double> ccoeffs(n_ccoeffs), coords(n_coords), exps(n_ccoeffs), norms(n_norms);
    vector<int> dims(n_pairs_2x), offsets(n_pairs_2x), poss_cntrs(n_pairs_2x), poss_norms(n_pairs_2x);

    int pos_cntrs = 0, pos_norms = 0;
    for (size_t ipair = 0; ipair < n_pairs; ipair++)
    {
        const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
        const auto &[coords_a, coords_b] = shell_pair_data.coords[ipair];
        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
        const auto &[norms_a, norms_b] = shell_pair_data.norms[ipair];
        const auto &[offset_a, offset_b] = shell_pair_data.offsets[ipair];

        dims[2 * ipair] = exps_a.size();
        dims[2 * ipair + 1] = exps_b.size();

        poss_cntrs[2 * ipair] = pos_cntrs;

        for (size_t i = 0; i < exps_a.size(); i++)
        {
            exps[pos_cntrs] = exps_a[i];
            ccoeffs[pos_cntrs] = ccoeffs_a[i];
            pos_cntrs++;
        }

        poss_cntrs[2 * ipair + 1] = pos_cntrs;

        for (size_t i = 0; i < exps_b.size(); i++)
        {
            exps[pos_cntrs] = exps_b[i];
            ccoeffs[pos_cntrs] = ccoeffs_b[i];
            pos_cntrs++;
        }

        coords[6 * ipair] = coords_a[0];
        coords[6 * ipair + 1] = coords_a[1];
        coords[6 * ipair + 2] = coords_a[2];
        coords[6 * ipair + 3] = coords_b[0];
        coords[6 * ipair + 4] = coords_b[1];
        coords[6 * ipair + 5] = coords_b[2];

        poss_norms[2 * ipair] = pos_norms;

        for (size_t i = 0; i < norms_a.size(); i++)
        {
            norms[pos_norms] = norms_a[i];
            pos_norms++;
        }

        poss_norms[2 * ipair + 1] = pos_norms;

        for (size_t i = 0; i < norms_b.size(); i++)
        {
            norms[pos_norms] = norms_b[i];
            pos_norms++;
        }

        offsets[2 * ipair] = offset_a;
        offsets[2 * ipair + 1] = offset_b;
    }

    // device data
    double *dev_ccoeffs, *dev_exps, *dev_coords, *dev_norms;
    int *dev_dims, *dev_offsets, *dev_poss_cntrs, *dev_poss_norms;

    hipCheck(hipMalloc(&dev_ccoeffs, n_ccoeffs * sizeof(double)));
    hipCheck(hipMalloc(&dev_exps, n_ccoeffs * sizeof(double)));
    hipCheck(hipMalloc(&dev_coords, n_coords * sizeof(double)));
    hipCheck(hipMalloc(&dev_norms, n_norms * sizeof(double)));
    hipCheck(hipMalloc(&dev_dims, n_pairs_2x * sizeof(int)));
    hipCheck(hipMalloc(&dev_offsets, n_pairs_2x * sizeof(int)));
    hipCheck(hipMalloc(&dev_poss_cntrs, n_pairs_2x * sizeof(int)));
    hipCheck(hipMalloc(&dev_poss_norms, n_pairs_2x * sizeof(int)));

    hipCheck(hipMemcpy(dev_ccoeffs, ccoeffs.data(), n_ccoeffs * sizeof(double),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_exps, exps.data(), n_ccoeffs * sizeof(double),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_coords, coords.data(), n_coords * sizeof(double),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_norms, norms.data(), n_norms * sizeof(double),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_dims, dims.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_offsets, offsets.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_poss_cntrs, poss_cntrs.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_poss_norms, poss_norms.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));

    // device cleanup
    hipCheck(hipFree(dev_ccoeffs));
    hipCheck(hipFree(dev_exps));
    hipCheck(hipFree(dev_coords));
    hipCheck(hipFree(dev_norms));
    hipCheck(hipFree(dev_dims));
    hipCheck(hipFree(dev_offsets));
    hipCheck(hipFree(dev_poss_cntrs));
    hipCheck(hipFree(dev_poss_norms));

    return overlap_ints;
}
