#include <lible/oneel_detail_gpu.hpp>
#include <lible/shell_pair_data.hpp>
#include <lible/util.hpp>

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

    // vector<int>
    size_t n_pairs = shell_pair_data.n_pairs;
    size_t n_pairs_2x = 2 * n_pairs;

    vector<double> ccoeffs(n_pairs_2x);
    vector<double> norms(n_pairs_2x);
    double *dev_norms;
    
    vector<double> dev_ccoeffs(n_pairs_2x);
    vector<double> dev_cexps(n_pairs_2x);
    vector<int> dev_offsets(n_pairs_2x);
    vector<int> dev_primitive_dimensions(n_pairs_2x);

    hipCheck(hipMalloc(&dev_norms, n_pairs_2x * sizeof(double)));

    hipCheck(hipMemcpy(dev_norms, norms.data(), n_pairs_2x * sizeof(double), hipMemcpyHostToDevice));

    size_t n_primitive_pairs = 0;
    for (size_t ipair = 0; ipair < shell_pair_data.n_pairs; ipair++)
    {
        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
        n_primitive_pairs += exps_a.size() * exps_b.size();
    }    

    printf("   n_primitive_pairs = %d\n", n_primitive_pairs);


    hipCheck(hipFree(dev_norms));

    return overlap_ints;
}