#include <lible/ints/gpu/gpuints.hpp>
#include <lible/util.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>

#include <chrono>
#include <iostream>
#include <math.h>

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

namespace LIG = lible::ints::gpu;

using std::vector;

// Monolithic kernel
__global__ void devKernelS00(const int dim_ao, const int n_shell_pairs, const double *ccoeffs,
                             const double *exps, const double *coords, const double *norms,
                             const int *cdepths, const int *offsets, const int *poss_cntrs,
                             const int *poss_norms, double *sints_out)
{
    int ipair = blockDim.x * blockIdx.x + threadIdx.x;
    if (ipair < n_shell_pairs)
    {
        int dim_a = cdepths[2 * ipair];
        int dim_b = cdepths[2 * ipair + 1];

        int pos_cntr_a = poss_cntrs[2 * ipair];
        int pos_cntr_b = poss_cntrs[2 * ipair + 1];

        double s000_sum = 0;
        for (int ia = 0; ia < dim_a; ia++)
            for (int ib = 0; ib < dim_b; ib++)
            {
                double xa = coords[6 * ipair];
                double ya = coords[6 * ipair + 1];
                double za = coords[6 * ipair + 2];
                double xb = coords[6 * ipair + 3];
                double yb = coords[6 * ipair + 4];
                double zb = coords[6 * ipair + 5];
                double r2_ab = pow(xa - xb, 2) + pow(ya - yb, 2) + pow(za - zb, 2);                

                double a = exps[pos_cntr_a + ia];
                double b = exps[pos_cntr_b + ib];
                double p = a + b;
                double mu = a * b / p;

                double fac = exp(-mu * r2_ab);
                double s000 = fac * pow(M_PI / p, 1.5);

                s000_sum += ccoeffs[pos_cntr_a + ia] * ccoeffs[pos_cntr_b + ib] * s000;
            }

        int pos_a = offsets[2 * ipair];
        int pos_b = offsets[2 * ipair + 1];
        int pos_ab = pos_a * dim_ao + pos_b;
        int pos_ba = pos_b * dim_ao + pos_a;

        int pos_norm_a = poss_norms[2 * ipair];
        int pos_norm_b = poss_norms[2 * ipair + 1];

        s000_sum = norms[pos_norm_a] * norms[pos_norm_b] * s000_sum;

        sints_out[pos_ab] = s000_sum;
        sints_out[pos_ba] = s000_sum;
    }
}

// // Grid-strided loop kernel
// __global__ void devKernelS00(const int dim_ao, const int n_shell_pairs, const double *ccoeffs,
//                              const double *exps, const double *coords, const double *norms,
//                              const int *cdepths, const int *offsets, const int *poss_cntrs,
//                              const int *poss_norms, double *sints_out)
// {
//     for (int ipair = blockIdx.x * blockDim.x + threadIdx.x;
//          ipair < n_shell_pairs;
//          ipair += blockDim.x * gridDim.x)
//     {
//         int dim_a = cdepths[2 * ipair];
//         int dim_b = cdepths[2 * ipair + 1];

//         int pos_cntr_a = poss_cntrs[2 * ipair];
//         int pos_cntr_b = poss_cntrs[2 * ipair + 1];

//         double s000_sum = 0;
//         for (int ia = 0; ia < dim_a; ia++)
//             for (int ib = 0; ib < dim_b; ib++)
//             {
//                 double xa = coords[6 * ipair];
//                 double ya = coords[6 * ipair + 1];
//                 double za = coords[6 * ipair + 2];
//                 double xb = coords[6 * ipair + 3];
//                 double yb = coords[6 * ipair + 4];
//                 double zb = coords[6 * ipair + 5];
//                 double r2_ab = pow(xa - xb, 2) + pow(ya - yb, 2) + pow(za - zb, 2);

//                 double a = exps[pos_cntr_a + ia];
//                 double b = exps[pos_cntr_b + ib];
//                 double p = a + b;
//                 double mu = a * b / p;

//                 double fac = exp(-mu * r2_ab);
//                 double s000 = fac * pow(M_PI / p, 1.5);

//                 s000_sum += ccoeffs[pos_cntr_a + ia] * ccoeffs[pos_cntr_b + ib] * s000;
//             }

//         int pos_a = offsets[2 * ipair];
//         int pos_b = offsets[2 * ipair + 1];
//         int pos_ab = pos_a * dim_ao + pos_b;
//         int pos_ba = pos_b * dim_ao + pos_a;

//         int pos_norm_a = poss_norms[2 * ipair];
//         int pos_norm_b = poss_norms[2 * ipair + 1];

//         s000_sum = norms[pos_norm_a] * norms[pos_norm_b] * s000_sum;

//         sints_out[pos_ab] = s000_sum;
//         sints_out[pos_ba] = s000_sum;
//     }
// }

lible::vec2d LIG::calculateS_L0(const Structure &structure)
{
    printf("LIG::calculateS_L0\n");

    auto t0_total = std::chrono::steady_clock::now();

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

    vector<double> ccoeffs(n_ccoeffs), coords(n_coords), exps(n_ccoeffs), norms(n_norms);
    vector<int> cdepths(n_pairs_2x), offsets(n_pairs_2x), poss_cntrs(n_pairs_2x), poss_norms(n_pairs_2x);

    int pos_cntrs = 0, pos_norms = 0;
    for (size_t ipair = 0; ipair < n_pairs; ipair++)
    {
        const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
        const auto &[coords_a, coords_b] = shell_pair_data.coords[ipair];
        const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
        const auto &[norms_a, norms_b] = shell_pair_data.norms[ipair];
        const auto &[offset_a, offset_b] = shell_pair_data.offsets[ipair];

        cdepths[2 * ipair] = exps_a.size();
        cdepths[2 * ipair + 1] = exps_b.size();

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
    int *dev_cdepths, *dev_offsets, *dev_poss_cntrs, *dev_poss_norms;

    hipCheck(hipMalloc(&dev_ccoeffs, n_ccoeffs * sizeof(double)));
    hipCheck(hipMalloc(&dev_exps, n_ccoeffs * sizeof(double)));
    hipCheck(hipMalloc(&dev_coords, n_coords * sizeof(double)));
    hipCheck(hipMalloc(&dev_norms, n_norms * sizeof(double)));
    hipCheck(hipMalloc(&dev_cdepths, n_pairs_2x * sizeof(int)));
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
    hipCheck(hipMemcpy(dev_cdepths, cdepths.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_offsets, offsets.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_poss_cntrs, poss_cntrs.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));
    hipCheck(hipMemcpy(dev_poss_norms, poss_norms.data(), n_pairs_2x * sizeof(int),
                       hipMemcpyHostToDevice));

    // print some bs
    hipDeviceProp_t dev_props;
    hipCheck(hipGetDeviceProperties(&dev_props, 0));

    std::cout << "dev_props.name                = " << dev_props.name << std::endl;
    std::cout << "dev_props.multiProcessorCount = " << dev_props.multiProcessorCount << std::endl;
    std::cout << "dev_props.major               = " << dev_props.major << std::endl;
    std::cout << "dev_props.minor               = " << dev_props.minor << std::endl;

    // kernel call
    auto t0_kernel = std::chrono::steady_clock::now();

    double *dev_overlap_ints;
    hipCheck(hipMalloc(&dev_overlap_ints, std::pow(dim_ao, 2) * sizeof(double)));
    hipCheck(hipMemcpy(dev_overlap_ints, overlap_ints.getData(),
                       std::pow(dim_ao, 2) * sizeof(double), hipMemcpyHostToDevice));

    auto t0_gpu = std::chrono::steady_clock::now();

    int n_thr_per_blk = 128; // Choosing this for now, it can vary depending on GPU...
    int n_blk_per_grd = std::ceil(double(n_pairs) / 128);

    devKernelS00<<<n_blk_per_grd, n_thr_per_blk>>>(structure.getDimAO(), n_pairs, dev_ccoeffs,
                                                   dev_exps, dev_coords, dev_norms, dev_cdepths,
                                                   dev_offsets, dev_poss_cntrs, dev_poss_norms,
                                                   dev_overlap_ints);

    // devKernelS00<<<6, 128>>>(structure.getDimAO(), n_pairs, dev_ccoeffs, dev_exps, dev_coords,
    //                          dev_norms, dev_cdepths, dev_offsets, dev_poss_cntrs, dev_poss_norms,
    //                          dev_overlap_ints);

    hipCheck(hipDeviceSynchronize());
    auto t1_gpu = std::chrono::steady_clock::now();

    hipCheck(hipMemcpy(overlap_ints.getData(), dev_overlap_ints, std::pow(dim_ao, 2) * sizeof(double),
                       hipMemcpyDeviceToHost));

    auto t1_kernel = std::chrono::steady_clock::now();

    // device cleanup
    hipCheck(hipFree(dev_ccoeffs));
    hipCheck(hipFree(dev_exps));
    hipCheck(hipFree(dev_coords));
    hipCheck(hipFree(dev_norms));
    hipCheck(hipFree(dev_cdepths));
    hipCheck(hipFree(dev_offsets));
    hipCheck(hipFree(dev_poss_cntrs));
    hipCheck(hipFree(dev_poss_norms));

    auto t1_total = std::chrono::steady_clock::now();

    std::chrono::duration<double> duration = t1_gpu - t0_gpu;
    palPrint(fmt::format("t(kernel): {:.2e} s\n", duration.count()));

    duration = t1_kernel - t0_kernel;
    palPrint(fmt::format("t(gpu):    {:.2e} s\n", duration.count()));

    duration = t1_total - t0_total;
    palPrint(fmt::format("t(total):  {:.2e} s\n", duration.count()));

    return overlap_ints;
}

// __device__ int idx(int la, int lb, int i, int j, int t)
// {
//     int dim2 = la + lb + 1;
//     int dim1 = (lb + 1) * dim2;

//     return i * dim1 * dim2 + j * dim2 + t;
// }

// __device__ int idx(int la, int lb, int i, int j, int t, int pos)
// {
//     int dim2 = la + lb + 1;
//     int dim1 = (lb + 1) * dim2;

//     return pos + i * dim1 * dim2 + j * dim2 + t;
// }

//////////////////////////////////////////////////////////////////////////////////////////////

// void LI::coeffs(const double a, const double b, const double PA, const double PB,
//                 const double one_o_2p, const int la, const int lb, vec3d &E)
// {
//     for (int i = 1; i <= la; i++)
//     {
//         E(i, 0, 0) = PA * E(i - 1, 0, 0) + E(i - 1, 0, 1);

//         for (int t = 1; t < i; t++)
//             E(i, 0, t) = one_o_2p * E(i - 1, 0, t - 1) +
//                          PA * E(i - 1, 0, t) +
//                          (t + 1) * E(i - 1, 0, t + 1);

//         E(i, 0, i) = one_o_2p * E(i - 1, 0, i - 1) + PA * E(i - 1, 0, i);
//     }

//     for (int j = 1; j <= lb; j++)
//         for (int i = 0; i <= la; i++)
//         {
//             E(i, j, 0) = PB * E(i, j - 1, 0) + E(i, j - 1, 1);

//             for (int t = 1; t < i + j; t++)
//                 E(i, j, t) = one_o_2p * E(i, j - 1, t - 1) +
//                              PB * E(i, j - 1, t) +
//                              (t + 1) * E(i, j - 1, t + 1);

//             E(i, j, i + j) = one_o_2p * E(i, j - 1, i + j - 1) + PB * E(i, j - 1, i + j);
//         }
// }

// __device__ int idx(int d3, int d23, int i, int j, int t, int pos)
// {
//     return pos + i * d23 + j * d3 + t;
// }

// __device__ int dimCartIdxs(int l)
// {
//     return (l + 1) * (l + 1) / 2;
// }

// __device__ void devKernelECoeffs_1(const int la, const int lb, const int pos, const double PA,
//                                    const double PB, const double one_o_2p, double *E)
// {
//     int d2 = lb + 1;
//     int d3 = la + lb + 1;
//     int d23 = d2 * d3;

//     E[idx(d3, d23, 0, 0, 0, pos)] = 0;

//     for (int i = 1; i <= la; i++)
//     {
//         E[idx(d3, d23, i, 0, 0, pos)] = PA * E[idx(d3, d23, i - 1, 0, 0, pos)] +
//                                         E[idx(d3, d23, i - 1, 0, 1, pos)];

//         for (int t = 1; t < i; t++)
//             E[idx(d3, d23, i, 0, t, pos)] = one_o_2p * E[idx(d3, d23, i - 1, 0, t - 1, pos)] +
//                                             PA * E[idx(d3, d23, i - 1, 0, t, pos)] +
//                                             (t + 1) * E[idx(d3, d23, i - 1, 0, t + 1, pos)];

//         E[idx(d3, d23, i, 0, i, pos)] = one_o_2p * E[idx(d3, d23, i - 1, 0, i - 1, pos)] +
//                                         PA * E[idx(d3, d23, i - 1, 0, i, pos)];
//     }

//     for (int j = 1; j <= lb; j++)
//         for (int i = 0; i <= la; i++)
//         {
//             E[idx(d3, d23, i, j, 0, pos)] = PB * E[idx(d3, d23, i, j - 1, 0, pos)] +
//                                             E[idx(d3, d23, i, j - 1, 1, pos)];

//             for (int t = 1; t < i + j; t++)
//                 E[idx(d3, d23, i, j, t, pos)] = one_o_2p * E[idx(d3, d23, i, j - 1, t - 1, pos)] +
//                                                 PB * E[idx(d3, d23, i, j - 1, t, pos)] +
//                                                 (t + 1) * E[idx(d3, d23, i, j - 1, t + 1, pos)];

//             E[idx(d3, d23, i, j, i + j, pos)] = one_o_2p * E[idx(d3, d23, i, j - 1, i + j - 1, pos)] +
//                                                 PB * E[idx(d3, d23, i, j - 1, i + j, pos)];
//         }
// }

// __global__ void devKernelECoeffs(const int la, const int lb, const int n_shell_pairs,
//                                  const int *cdepths, const int *eoffsets, const int *poss_cntrs,
//                                  const double *coords, const double *exps, double *ecoeffs)
// {
//     int id = blockDim.x * blockIdx.x + threadIdx.x;
//     if (id < n_shell_pairs)
//     {
//         int ipair = id;

//         int n_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1);
//         int n_ecoeffs_3x = 3 * n_ecoeffs;

//         int dim_cntr_a = cdepths[2 * ipair];
//         int dim_cntr_b = cdepths[2 * ipair + 1];
//         int pos_cntr_a = poss_cntrs[2 * ipair];
//         int pos_cntr_b = poss_cntrs[2 * ipair + 1];

//         double xa = coords[6 * ipair];
//         double ya = coords[6 * ipair + 1];
//         double za = coords[6 * ipair + 2];
//         double xb = coords[6 * ipair + 3];
//         double yb = coords[6 * ipair + 4];
//         double zb = coords[6 * ipair + 5];
//         double rab2[3] = {std::pow(xa - ya, 2), std::pow(xb - yb, 2), std::pow(za - zb, 2)};

//         for (int ia = 0, iab = 0; ia < dim_cntr_a; ia++)
//             for (int ib = 0; ib < dim_cntr_b; ib++, iab++)
//             {
//                 double a = exps[pos_cntr_a + ia];
//                 double b = exps[pos_cntr_b + ib];

//                 double p = a + b;
//                 double one_o_2p = 1.0 / (2 * p);

//                 double P[3] = {(a * xa + b * xb) / p, (a * ya + b * yb) / p, (a * za + b * zb) / p};
//                 double PA[3] = {P[0] - xa, P[1] - ya, P[2] - za};
//                 double PB[3] = {P[0] - xb, P[1] - yb, P[2] - zb};

//                 double mu = a * b / (a + b);
//                 double Kab[3] = {std::exp(-mu * rab2[0]), std::exp(-mu * rab2[1]), std::exp(-mu * rab2[2])};

//                 int pos_ecoeffs = eoffsets[id] + iab * n_ecoeffs_3x;
//                 int pos_x = pos_ecoeffs;
//                 int pos_y = pos_ecoeffs + n_ecoeffs;
//                 int pos_z = pos_ecoeffs + 2 * n_ecoeffs;

//                 devKernelECoeffs_1(la, lb, pos_x, PA[0], PB[0], one_o_2p, ecoeffs);
//                 devKernelECoeffs_1(la, lb, pos_y, PA[1], PB[1], one_o_2p, ecoeffs);
//                 devKernelECoeffs_1(la, lb, pos_z, PA[2], PB[2], one_o_2p, ecoeffs);
//             }
//     }
// }

// __global__ void devKernelOverlap(const int dim_ao_cart, const int la, const int lb,
//                                  const int n_shell_pairs, const int *cart_idxs_poss_a,
//                                  const int *cart_idxs_poss_b, const int *cdepths,
//                                  const int *eoffsets, const int *offsets_cart,
//                                  const int *poss_cntrs, const double *ccoeffs,
//                                  const double *exps, const double *ecoeffs,
//                                  double *sints_cart)
// {
//     int id = blockDim.x * blockIdx.x + threadIdx.x;
//     if (id < n_shell_pairs)
//     {
//         int ipair = id;

//         int d2 = lb + 1;
//         int d3 = la + lb + 1;
//         int d23 = d2 * d3;
//         int n_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1);
//         int n_ecoeffs_3x = 3 * n_ecoeffs;

//         int dim_cart_idxs_a = dimCartIdxs(la);
//         int dim_cart_idxs_b = dimCartIdxs(lb);
//         int dim_cntr_a = cdepths[2 * ipair];
//         int dim_cntr_b = cdepths[2 * ipair + 1];

//         int offset_a = offsets_cart[2 * ipair];
//         int offset_b = offsets_cart[2 * ipair + 1];

//         int pos_cntr_a = poss_cntrs[2 * ipair];
//         int pos_cntr_b = poss_cntrs[2 * ipair + 1];

//         for (int mu = 0; mu < dim_cart_idxs_a; mu++)
//             for (int nu = 0; nu < dim_cart_idxs_b; nu++)
//             {
//                 int pos_a = offset_a + mu;
//                 int pos_b = offset_b + nu;
//                 int pos_ab = pos_a * dim_ao_cart + pos_b;
//                 sints_cart[pos_ab] = 0;
//             }

//         for (int ia = 0, iab = 0; ia < dim_cntr_a; ia++)
//             for (int ib = 0; ib < dim_cntr_b; ib++, iab++)
//             {
//                 double a = exps[pos_cntr_a + ia];
//                 double b = exps[pos_cntr_b + ib];
//                 double p = a + b;

//                 int pos_ecoeffs = eoffsets[id] + iab * n_ecoeffs_3x;
//                 int pos_x = pos_ecoeffs;
//                 int pos_y = pos_ecoeffs + n_ecoeffs;
//                 int pos_z = pos_ecoeffs + 2 * n_ecoeffs;

//                 double val = pow(M_PI / p, 1.5);

//                 for (int mu = 0; mu < dim_cart_idxs_a; mu++)
//                 {
//                     int i = cart_idxs_poss_a[3 * mu];
//                     int j = cart_idxs_poss_a[3 * mu + 1];
//                     int k = cart_idxs_poss_a[3 * mu + 2];

//                     for (int nu = 0; nu < dim_cart_idxs_b; nu++)
//                     {
//                         int i_ = cart_idxs_poss_b[3 * nu];
//                         int j_ = cart_idxs_poss_b[3 * nu + 1];
//                         int k_ = cart_idxs_poss_b[3 * nu + 2];

//                         double ecoeff = ecoeffs[idx(d3, d23, i, i_, 0, pos_x)] *
//                                         ecoeffs[idx(d3, d23, j, j_, 0, pos_y)] *
//                                         ecoeffs[idx(d3, d23, k, k_, 0, pos_z)];

//                         int mu_ = offset_a + mu;
//                         int nu_ = offset_b + nu;
//                         int idx = mu_ * dim_ao_cart + nu_;
//                         sints_cart[idx] += ccoeffs[pos_cntr_a + ia] * ccoeffs[pos_cntr_b + ib] *
//                                            ecoeff * val;
//                     }
//                 }
//             }
//     }
// }

// lible::vec2d LIG::calculateS(const Structure &structure)
// {
//     // Make shell pair datas
//     // Collect info about dimensions
//     // Preallocate CPU/GPU data
//     // Copy stuff on GPU and call kernels

//     printf("LIG::calculateS\n");

//     auto t0_total = std::chrono::steady_clock::now();

//     int l_max = structure.getMaxL();
//     int l_pairs = (l_max + 1) * (l_max + 2) / 2;

//     vector<ShellPairData> shell_pair_datas(l_pairs);
//     for (int la = l_max, idx = 0; la >= 0; la--)
//         for (int lb = la; lb >= 0; lb--, idx++)
//             shell_pair_datas[idx] = ShellPairData(la, lb, structure);

//     size_t n_pairs_max{0}, n_pairs_2x_max{0}, n_coords_max{0}, n_ccoeffs_max{0}, n_norms_max{0};
//     size_t n_ecoeffs_max{0};
//     for (int la = l_max, idx = 0; la >= 0; la--)
//         for (int lb = la; lb >= 0; lb--, idx++)
//         {
//             const ShellPairData &shell_pair_data{shell_pair_datas[idx]};

//             size_t n_pairs = shell_pair_data.n_pairs;
//             if (n_pairs > n_pairs_max)
//                 n_pairs_max = n_pairs;

//             size_t n_pairs_2x = 2 * n_pairs;
//             if (n_pairs_2x > n_pairs_2x_max)
//                 n_pairs_2x_max = n_pairs_2x;

//             size_t n_coords = 6 * n_pairs;
//             if (n_coords > n_coords_max)
//                 n_coords_max = n_coords;

//             size_t dim_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1) * 3;

//             size_t n_ccoeffs{0}, n_norms{0}, n_ecoeffs{0};
//             for (size_t ipair = 0; ipair < n_pairs; ipair++)
//             {
//                 const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
//                 const auto &[norms_a, norms_b] = shell_pair_data.norms[ipair];
//                 n_ccoeffs += ccoeffs_a.size() + ccoeffs_b.size();
//                 n_norms += norms_a.size() + norms_b.size();
//                 n_ecoeffs += ccoeffs_a.size() * ccoeffs_b.size() * dim_ecoeffs;
//             }

//             if (n_ccoeffs > n_ccoeffs_max)
//                 n_ccoeffs_max = n_ccoeffs;

//             if (n_norms > n_norms_max)
//                 n_norms_max = n_norms;
            
//             if (n_ecoeffs > n_ecoeffs_max)
//                 n_ecoeffs_max = n_ecoeffs;
//         }

//     // Allocate CPU data
//     vector<double> ccoeffs(n_ccoeffs_max), coords(n_coords_max), exps(n_ccoeffs_max),
//         norms(n_norms_max);

//     vector<int> cdepths(n_pairs_2x_max), offsets_cart(n_pairs_2x_max),
//         offsets_sph(n_pairs_2x_max), poss_cntrs(n_pairs_2x_max), poss_norms(n_pairs_2x_max);

//     vector<int> offsets_ecoeffs(n_pairs_max);

//     // Allocate GPU data
//     double *dev_ccoeffs, *dev_exps, *dev_coords, *dev_norms;
//     int *dev_cdepths, *dev_offsets_cart, *dev_offsets_sph, *dev_poss_cntrs, *dev_poss_norms;

//     double *dev_ecoeffs;
//     int *dev_offsets_ecoeffs; 

//     size_t dim_ao = structure.getDimAO();
//     size_t dim_ao_cart = structure.getDimAOCart();

//     size_t bytes_gpu = n_ccoeffs_max * sizeof(double);
//     bytes_gpu += n_ccoeffs_max * sizeof(double);
//     bytes_gpu += n_coords_max * sizeof(double);
//     bytes_gpu += n_norms_max * sizeof(double);
//     bytes_gpu += n_pairs_2x_max * sizeof(int);
//     bytes_gpu += n_pairs_2x_max * sizeof(int);
//     bytes_gpu += n_pairs_2x_max * sizeof(int);
//     bytes_gpu += n_pairs_2x_max * sizeof(int);
//     bytes_gpu += n_pairs_2x_max * sizeof(int);
//     bytes_gpu += n_ecoeffs_max * sizeof(double);
//     bytes_gpu += n_pairs_max * sizeof(int);
//     bytes_gpu += std::pow(dim_ao, 2) * sizeof(double);
//     bytes_gpu += std::pow(dim_ao_cart, 2) * sizeof(double);

//     hipCheck(hipMalloc(&dev_ccoeffs, n_ccoeffs_max * sizeof(double)));
//     hipCheck(hipMalloc(&dev_exps, n_ccoeffs_max * sizeof(double)));
//     hipCheck(hipMalloc(&dev_coords, n_coords_max * sizeof(double)));
//     hipCheck(hipMalloc(&dev_norms, n_norms_max * sizeof(double)));
//     hipCheck(hipMalloc(&dev_cdepths, n_pairs_2x_max * sizeof(int)));
//     hipCheck(hipMalloc(&dev_offsets_cart, n_pairs_2x_max * sizeof(int)));
//     hipCheck(hipMalloc(&dev_offsets_sph, n_pairs_2x_max * sizeof(int)));
//     hipCheck(hipMalloc(&dev_poss_cntrs, n_pairs_2x_max * sizeof(int)));
//     hipCheck(hipMalloc(&dev_poss_norms, n_pairs_2x_max * sizeof(int)));

//     hipCheck(hipMalloc(&dev_ecoeffs, n_ecoeffs_max * sizeof(double)));
//     hipCheck(hipMalloc(&dev_offsets_ecoeffs, n_pairs_max * sizeof(int)));

//     vec2d overlap_ints(dim_ao, dim_ao, 0);

//     double *dev_overlap_ints_cart, *dev_overlap_ints_sph;
//     hipCheck(hipMalloc(&dev_overlap_ints_sph, std::pow(dim_ao, 2) * sizeof(double)));
//     hipCheck(hipMalloc(&dev_overlap_ints_cart, std::pow(dim_ao_cart, 2) * sizeof(double)));    

//     // Copy data to GPU and call kernels
//     for (int la = l_max, idx = 0; la >= 0; la--)
//         for (int lb = la; lb >= 0; lb--, idx++)
//         {
//             const ShellPairData &shell_pair_data{shell_pair_datas[idx]};

//             size_t n_pairs = shell_pair_data.n_pairs;
//             size_t n_pairs_2x = 2 * n_pairs;

//             size_t n_ccoeffs = 0;
//             size_t n_coords = 0;
//             size_t n_norms = 0;
//             size_t pos_cntrs = 0;
//             size_t pos_ecoeffs = 0;
//             size_t pos_norms = 0;

//             // Set up host data
//             for (size_t ipair = 0; ipair < n_pairs; ipair++)
//             {
//                 const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
//                 const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
//                 const auto &[coords_a, coords_b] = shell_pair_data.coords[ipair];
//                 const auto &[norms_a, norms_b] = shell_pair_data.norms[ipair];
//                 const auto &[offset_a, offset_b] = shell_pair_data.offsets[ipair];
//                 const auto &[offset_a_cart, offset_b_cart] = shell_pair_data.offsets_cart[ipair];

//                 n_ccoeffs += exps_a.size() + exps_b.size();

//                 cdepths[2 * ipair] = exps_a.size();
//                 cdepths[2 * ipair + 1] = exps_b.size();

//                 poss_cntrs[2 * ipair] = pos_cntrs;
//                 for (size_t i = 0; i < exps_a.size(); i++)
//                 {
//                     ccoeffs[pos_cntrs] = ccoeffs_a[i];
//                     exps[pos_cntrs] = exps_a[i];
//                     pos_cntrs++;
//                 }

//                 poss_cntrs[2 * ipair + 1] = pos_cntrs;
//                 for (size_t i = 0; i < exps_b.size(); i++)
//                 {
//                     ccoeffs[pos_cntrs] = ccoeffs_b[i];
//                     exps[pos_cntrs] = exps_b[i];
//                     pos_cntrs++;
//                 }

//                 n_coords += 6;
//                 coords[6 * ipair] = coords_a[0];
//                 coords[6 * ipair + 1] = coords_a[1];
//                 coords[6 * ipair + 2] = coords_a[2];
//                 coords[6 * ipair] = coords_b[0];
//                 coords[6 * ipair + 1] = coords_b[1];
//                 coords[6 * ipair + 2] = coords_b[2];

//                 n_norms += norms_a.size() + norms_b.size();

//                 poss_norms[2 * ipair] = pos_norms;
//                 for (size_t i = 0; i < norms_a.size(); i++)
//                 {
//                     norms[pos_norms] = norms_a[i];
//                     pos_norms++;
//                 }

//                 poss_norms[2 * ipair + 1] = pos_norms;
//                 for (size_t i = 0; i < norms_b.size(); i++)
//                 {
//                     norms[pos_norms] = norms_b[i];
//                     pos_norms++;
//                 }

//                 offsets_sph[2 * ipair] = offset_a;
//                 offsets_sph[2 * ipair + 1] = offset_b;

//                 offsets_cart[2 * ipair] = offset_a_cart;
//                 offsets_cart[2 * ipair + 1] = offset_b_cart;                            
                        
//                 offsets_ecoeffs[ipair] = pos_ecoeffs;

//                 size_t dim_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1) * 3;        
//                 pos_ecoeffs += dim_ecoeffs * ccoeffs_a.size() * ccoeffs_b.size();
//             }

//             printf("n_ccoeffs = %zu\n", n_ccoeffs);

//             // Copy data from host to CPU
//             hipCheck(hipMemcpy(dev_ccoeffs, ccoeffs.data(), n_ccoeffs * sizeof(double), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_exps, exps.data(), n_ccoeffs * sizeof(double), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_cdepths, cdepths.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_coords, coords.data(), n_coords * sizeof(double), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_norms, norms.data(), n_norms * sizeof(double), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_poss_norms, poss_norms.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_offsets_sph, offsets_sph.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_offsets_cart, offsets_cart.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
//             hipCheck(hipMemcpy(dev_offsets_ecoeffs, offsets_ecoeffs.data(), n_pairs * sizeof(int), hipMemcpyHostToDevice));

//             // Call the GPU kernels   
//             int n_thr_per_blk = 128; // Choosing this for now, in the future, will use a more flexible approach
//             int n_blk_per_grd = std::ceil(double(n_pairs) / n_thr_per_blk);

//             devKernelECoeffs<<<n_blk_per_grd, n_thr_per_blk>>>(la, lb, n_pairs, dev_cdepths,
//                                                                dev_offsets_ecoeffs, dev_poss_cntrs,
//                                                                dev_coords, dev_exps, dev_ecoeffs);
            
//             // 1) overlap ints
//             // devKernelOverlap<<<n_blk_per_grd, n_thr_per_blk>>>(dim_ao_cart, la, lb, n_pairs, );
//             // 2) spherical trafo
//         }    

//     // Release GPU data
//     hipCheck(hipFree(dev_ccoeffs));
//     hipCheck(hipFree(dev_exps));
//     hipCheck(hipFree(dev_coords));
//     hipCheck(hipFree(dev_norms));
//     hipCheck(hipFree(dev_cdepths));
//     hipCheck(hipFree(dev_offsets_cart));
//     hipCheck(hipFree(dev_offsets_sph));
//     hipCheck(hipFree(dev_poss_cntrs));
//     hipCheck(hipFree(dev_poss_norms));

//     hipCheck(hipFree(dev_ecoeffs));
//     hipCheck(hipFree(dev_offsets_ecoeffs));

//     auto t1_total = std::chrono::steady_clock::now();
//     std::chrono::duration<double> duration = t1_total - t0_total;

//     printf("bytes_gpu = %zu\n", bytes_gpu);
//     palPrint(fmt::format("t(total): {:.2e}\n", duration.count()));

//     return overlap_ints;
// }

// lible::vec2d LIG::calculateS(const Structure &structure)
// {
//     palPrint("LIG::calculateS\n");

//     auto t0_total = std::chrono::steady_clock::now();

//     int l_max = structure.getMaxL();
//     int l_pairs = (l_max + 1) * (l_max + 2) / 2;
//     size_t dim_ao = structure.getDimAO();

//     vector<ShellPairData> shell_pair_datas(l_pairs);
//     for (int la = l_max, idx = 0; la >= 0; la--)
//         for (int lb = la; lb >= 0; lb--, idx++)
//             shell_pair_datas[idx] = ShellPairData(la, lb, structure);

//     size_t max_n_pairs = 0;
//     size_t max_n_ecoeffs = 0;
//     size_t max_n_ecoeffs_raw = 0;
//     vector<vector<int>> ecoeff_offsets_all(l_pairs);
//     for (int ispd = 0; ispd < l_pairs; ispd++)
//     {
//         const auto &shell_pair_data = shell_pair_datas[ispd];
//         int la = shell_pair_data.la;
//         int lb = shell_pair_data.lb;
//         int dim_cart_a = dimCartesians(la);
//         int dim_cart_b = dimCartesians(lb);
//         int dim_cart_axb = dim_cart_a * dim_cart_b;
//         int dim_ecoeffs_raw = (la + 1) * (lb + 1) * (la + lb + 1) * 3;

//         size_t n_pairs = shell_pair_data.n_pairs;

//         size_t n_ecoeffs = 0;
//         size_t n_ecoeffs_raw = 0;

//         int ecoeff_offset = 0;
//         vector<int> ecoeff_offsets(n_pairs);
//         for (size_t ipair = 0; ipair < n_pairs; ipair++)
//         {
//             ecoeff_offsets[ipair] = ecoeff_offset;

//             const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
//             size_t cdepth_a = exps_a.size();
//             size_t cdepth_b = exps_b.size();

//             size_t n_prim_pairs = cdepth_a * cdepth_b;
//             n_ecoeffs += dim_cart_axb * n_prim_pairs;
//             n_ecoeffs_raw += (la + 1) * (lb + 1) * (la + lb + 1) * 3 * n_prim_pairs;

//             ecoeff_offset += dim_ecoeffs_raw * n_prim_pairs;
//         }

//         if (n_pairs > max_n_pairs)
//             max_n_pairs = n_pairs;

//         if (n_ecoeffs > max_n_ecoeffs)
//             max_n_ecoeffs = n_ecoeffs;

//         if (n_ecoeffs_raw > max_n_ecoeffs_raw)
//             max_n_ecoeffs_raw = n_ecoeffs_raw;

//         ecoeff_offsets_all[ispd] = ecoeff_offsets;

//         printf("la = %d, lb = %d, n_pairs = %zu, n_ecoeffs = %zu, n_ecoeffs_raw = %zu\n",
//                la, lb, n_pairs, n_ecoeffs, n_ecoeffs_raw);
//     }

//     printf("max_n_ecoeffs     = %zu\n", max_n_ecoeffs);
//     printf("max_n_ecoeffs_raw = %zu\n", max_n_ecoeffs_raw);

//     vec2d overlap_ints(dim_ao, dim_ao, 0);

//     // Allocate memory on the device
//     double *dev_ecoeffs;
//     int *dev_ecoeff_offsets;

//     hipCheck(hipMalloc(&dev_ecoeffs, max_n_ecoeffs_raw * sizeof(double)));
//     hipCheck(hipMalloc(&dev_ecoeff_offsets, max_n_pairs * sizeof(int)));

//     // Copy data and call the kernels to calculate ints
//     for (int la = l_max, ispd = 0; la >= 0; la--)
//         for (int lb = la; lb >= 0; lb--, ispd++)        
//         {
//             const ShellPairData &shell_pair_data = shell_pair_datas[ispd];

//             size_t n_pairs = shell_pair_data.n_pairs;

//             size_t n_prims = 0;
//             for (size_t ipair = 0; ipair < n_pairs; ipair++)
//             {
//                 const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
//                 n_prims += exps_a.size(); 
//                 n_prims += exps_b.size();
//             }

//             vector<int> &ecoeff_offsets = ecoeff_offsets_all[ispd];

//             hipCheck(hipMemcpy(dev_ecoeff_offsets, ecoeff_offsets.data(),
//                                ecoeff_offsets.size() * sizeof(int), hipMemcpyHostToDevice));
//         }

//     hipCheck(hipFree(dev_ecoeffs));
//     hipCheck(hipFree(dev_ecoeff_offsets));

//     auto t1_total = std::chrono::steady_clock::now();

//     std::chrono::duration<double> duration = t1_total - t0_total;
//     palPrint(fmt::format("t(total):  {:.2e} s\n", duration.count()));

//     return overlap_ints;
// }