#include <lible/ints/gpu/gpuints.hpp>
#include <lible/ints/gpu/dev_ecoeffs.hpp>
#include <lible/ints/gpu/utils.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/util.hpp>

#include <hip/hip_runtime.h>
#include <fmt/core.h>

namespace LIG = lible::ints::gpu;

using std::vector;

template <>
lible::vec2d LIG::calculate<LIG::Option::overlap>(const Structure &structure)
{
    // Make shell pair datas
    // Collect info about dimensions
    // Preallocate CPU/GPU data
    // Copy stuff on GPU and call kernels

    printf("LIG::calculateS\n");

    auto t0_total = std::chrono::steady_clock::now();

    int l_max = structure.getMaxL();
    int n_l_pairs = (l_max + 1) * (l_max + 2) / 2;

    vector<ShellPairData> shell_pair_datas(n_l_pairs);
    for (int la = l_max, idx = 0; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--, idx++)
            shell_pair_datas[idx] = ShellPairData(la, lb, structure);

    size_t n_pairs_max{0}, n_pairs_2x_max{0}, n_coords_max{0}, n_ccoeffs_max{0}, n_norms_max{0};
    size_t n_ecoeffs_max{0};
    for (int la = l_max, idx = 0; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--, idx++)
        {
            const ShellPairData &shell_pair_data{shell_pair_datas[idx]};

            size_t n_pairs = shell_pair_data.n_pairs;
            if (n_pairs > n_pairs_max)
                n_pairs_max = n_pairs;

            size_t n_pairs_2x = 2 * n_pairs;
            if (n_pairs_2x > n_pairs_2x_max)
                n_pairs_2x_max = n_pairs_2x;

            size_t n_coords = 6 * n_pairs;
            if (n_coords > n_coords_max)
                n_coords_max = n_coords;

            size_t dim_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1) * 3;
            size_t n_ccoeffs{0}, n_norms{0}, n_ecoeffs{0};
            for (size_t ipair = 0; ipair < n_pairs; ipair++)
            {
                const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
                const auto &[norms_a, norms_b] = shell_pair_data.norms[ipair];
                n_ccoeffs += ccoeffs_a.size() + ccoeffs_b.size();
                n_norms += norms_a.size() + norms_b.size();
                n_ecoeffs += ccoeffs_a.size() * ccoeffs_b.size() * dim_ecoeffs;
            }

            if (n_ccoeffs > n_ccoeffs_max)
                n_ccoeffs_max = n_ccoeffs;

            if (n_norms > n_norms_max)
                n_norms_max = n_norms;

            if (n_ecoeffs > n_ecoeffs_max)
                n_ecoeffs_max = n_ecoeffs;
        }

    // Allocate CPU data
    vector<double> ccoeffs(n_ccoeffs_max), coords(n_coords_max), exps(n_ccoeffs_max),
        norms(n_norms_max);

    vector<int> cdepths(n_pairs_2x_max), offsets_cart(n_pairs_2x_max),
        offsets_sph(n_pairs_2x_max), poss_cntrs(n_pairs_2x_max), poss_norms(n_pairs_2x_max);

    vector<int> offsets_ecoeffs(n_pairs_max);

    // Allocate GPU data
    double *dev_ccoeffs, *dev_exps, *dev_coords, *dev_norms;
    int *dev_cdepths, *dev_offsets_cart, *dev_offsets_sph, *dev_poss_cntrs, *dev_poss_norms;

    double *dev_ecoeffs;
    int *dev_offsets_ecoeffs;

    size_t dim_ao = structure.getDimAO();
    size_t dim_ao_cart = structure.getDimAOCart();

    size_t bytes_gpu = n_ccoeffs_max * sizeof(double);
    bytes_gpu += n_ccoeffs_max * sizeof(double);
    bytes_gpu += n_coords_max * sizeof(double);
    bytes_gpu += n_norms_max * sizeof(double);
    bytes_gpu += n_pairs_2x_max * sizeof(int);
    bytes_gpu += n_pairs_2x_max * sizeof(int);
    bytes_gpu += n_pairs_2x_max * sizeof(int);
    bytes_gpu += n_pairs_2x_max * sizeof(int);
    bytes_gpu += n_pairs_2x_max * sizeof(int);
    bytes_gpu += n_ecoeffs_max * sizeof(double);
    bytes_gpu += n_pairs_max * sizeof(int);
    bytes_gpu += std::pow(dim_ao, 2) * sizeof(double);
    bytes_gpu += std::pow(dim_ao_cart, 2) * sizeof(double);

    hipCheck(hipMalloc(&dev_ccoeffs, n_ccoeffs_max * sizeof(double)));
    hipCheck(hipMalloc(&dev_exps, n_ccoeffs_max * sizeof(double)));
    hipCheck(hipMalloc(&dev_coords, n_coords_max * sizeof(double)));
    hipCheck(hipMalloc(&dev_norms, n_norms_max * sizeof(double)));
    hipCheck(hipMalloc(&dev_cdepths, n_pairs_2x_max * sizeof(int)));
    hipCheck(hipMalloc(&dev_offsets_cart, n_pairs_2x_max * sizeof(int)));
    hipCheck(hipMalloc(&dev_offsets_sph, n_pairs_2x_max * sizeof(int)));
    hipCheck(hipMalloc(&dev_poss_cntrs, n_pairs_2x_max * sizeof(int)));
    hipCheck(hipMalloc(&dev_poss_norms, n_pairs_2x_max * sizeof(int)));

    hipCheck(hipMalloc(&dev_ecoeffs, n_ecoeffs_max * sizeof(double)));
    hipCheck(hipMalloc(&dev_offsets_ecoeffs, n_pairs_max * sizeof(int)));

    vec2d overlap_ints(dim_ao, dim_ao, 0);

    double *dev_overlap_ints_cart, *dev_overlap_ints_sph;
    hipCheck(hipMalloc(&dev_overlap_ints_sph, std::pow(dim_ao, 2) * sizeof(double)));
    hipCheck(hipMalloc(&dev_overlap_ints_cart, std::pow(dim_ao_cart, 2) * sizeof(double)));

    // Copy data to GPU and call kernels
    for (int la = l_max, idx = 0; la >= 0; la--)
        for (int lb = la; lb >= 0; lb--, idx++)
        {
            const ShellPairData &shell_pair_data{shell_pair_datas[idx]};

            size_t n_pairs = shell_pair_data.n_pairs;
            size_t n_pairs_2x = 2 * n_pairs;

            size_t n_ccoeffs = 0;
            size_t n_coords = 0;
            size_t n_norms = 0;
            size_t pos_cntrs = 0;
            size_t pos_ecoeffs = 0;
            size_t pos_norms = 0;

            // Set up host data
            for (size_t ipair = 0; ipair < n_pairs; ipair++)
            {
                const auto &[ccoeffs_a, ccoeffs_b] = shell_pair_data.ccoeffs[ipair];
                const auto &[exps_a, exps_b] = shell_pair_data.exps[ipair];
                const auto &[coords_a, coords_b] = shell_pair_data.coords[ipair];
                const auto &[norms_a, norms_b] = shell_pair_data.norms[ipair];
                const auto &[offset_a, offset_b] = shell_pair_data.offsets[ipair];
                const auto &[offset_a_cart, offset_b_cart] = shell_pair_data.offsets_cart[ipair];

                n_ccoeffs += exps_a.size() + exps_b.size();

                cdepths[2 * ipair] = exps_a.size();
                cdepths[2 * ipair + 1] = exps_b.size();

                poss_cntrs[2 * ipair] = pos_cntrs;
                for (size_t i = 0; i < exps_a.size(); i++)
                {
                    ccoeffs[pos_cntrs] = ccoeffs_a[i];
                    exps[pos_cntrs] = exps_a[i];
                    pos_cntrs++;
                }

                poss_cntrs[2 * ipair + 1] = pos_cntrs;
                for (size_t i = 0; i < exps_b.size(); i++)
                {
                    ccoeffs[pos_cntrs] = ccoeffs_b[i];
                    exps[pos_cntrs] = exps_b[i];
                    pos_cntrs++;
                }

                n_coords += 6;
                coords[6 * ipair] = coords_a[0];
                coords[6 * ipair + 1] = coords_a[1];
                coords[6 * ipair + 2] = coords_a[2];
                coords[6 * ipair] = coords_b[0];
                coords[6 * ipair + 1] = coords_b[1];
                coords[6 * ipair + 2] = coords_b[2];

                n_norms += norms_a.size() + norms_b.size();

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

                offsets_sph[2 * ipair] = offset_a;
                offsets_sph[2 * ipair + 1] = offset_b;

                offsets_cart[2 * ipair] = offset_a_cart;
                offsets_cart[2 * ipair + 1] = offset_b_cart;

                offsets_ecoeffs[ipair] = pos_ecoeffs;
                size_t dim_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1) * 3;
                pos_ecoeffs += dim_ecoeffs * ccoeffs_a.size() * ccoeffs_b.size();
            }

            printf("n_ccoeffs = %zu\n", n_ccoeffs);

            // Copy data from host to CPU
            hipCheck(hipMemcpy(dev_ccoeffs, ccoeffs.data(), n_ccoeffs * sizeof(double), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_exps, exps.data(), n_ccoeffs * sizeof(double), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_cdepths, cdepths.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_coords, coords.data(), n_coords * sizeof(double), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_norms, norms.data(), n_norms * sizeof(double), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_poss_norms, poss_norms.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_offsets_sph, offsets_sph.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_offsets_cart, offsets_cart.data(), n_pairs_2x * sizeof(int), hipMemcpyHostToDevice));
            hipCheck(hipMemcpy(dev_offsets_ecoeffs, offsets_ecoeffs.data(), n_pairs * sizeof(int), hipMemcpyHostToDevice));

            // Call the GPU kernels
            int n_thr_per_blk = 128; // Choosing this for now, in the future, will use a more flexible approach
            int n_blk_per_grd = std::ceil(double(n_pairs) / n_thr_per_blk);

            calcECoeffs<<<n_blk_per_grd, n_thr_per_blk>>>(la, lb, n_pairs, dev_cdepths,
                                                          dev_offsets_ecoeffs, dev_poss_cntrs,
                                                          dev_coords, dev_exps, dev_ecoeffs);

            // 1) overlap ints
            // devKernelOverlap<<<n_blk_per_grd, n_thr_per_blk>>>(dim_ao_cart, la, lb, n_pairs, );
            // 2) spherical trafo
        }

    // // Release GPU data
    // hipCheck(hipFree(dev_ccoeffs));
    // hipCheck(hipFree(dev_exps));
    // hipCheck(hipFree(dev_coords));
    // hipCheck(hipFree(dev_norms));
    // hipCheck(hipFree(dev_cdepths));
    // hipCheck(hipFree(dev_offsets_cart));
    // hipCheck(hipFree(dev_offsets_sph));
    // hipCheck(hipFree(dev_poss_cntrs));
    // hipCheck(hipFree(dev_poss_norms));

    // hipCheck(hipFree(dev_ecoeffs));
    // hipCheck(hipFree(dev_offsets_ecoeffs));

    // auto t1_total = std::chrono::steady_clock::now();
    // std::chrono::duration<double> duration = t1_total - t0_total;

    // printf("bytes_gpu = %zu\n", bytes_gpu);
    // palPrint(fmt::format("t(total): {:.2e}\n", duration.count()));

    return overlap_ints;
}