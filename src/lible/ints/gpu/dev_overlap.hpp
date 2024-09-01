#pragma once

#include <hip/hip_runtime.h>

namespace lible::ints::gpu
{
    __global__ inline void devCalcOverlap(const int dim_ao_cart, const int la, const int lb,
                                          const int n_shell_pairs, const int *cart_idxs_poss_a,
                                          const int *cart_idxs_poss_b, const int *cdepths,
                                          const int *eoffsets, const int *offsets_cart,
                                          const int *poss_cntrs, const double *ccoeffs,
                                          const double *exps, const double *ecoeffs,
                                          double *sints_cart)
    {
        // int id = blockDim.x * blockIdx.x + threadIdx.x;
        // if (id < n_shell_pairs)
        // {
        //     int ipair = id;

        //     int d2 = lb + 1;
        //     int d3 = la + lb + 1;
        //     int d23 = d2 * d3;
        //     int n_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1);
        //     int n_ecoeffs_3x = 3 * n_ecoeffs;

        //     int dim_cart_idxs_a = dimCartIdxs(la);
        //     int dim_cart_idxs_b = dimCartIdxs(lb);
        //     int dim_cntr_a = cdepths[2 * ipair];
        //     int dim_cntr_b = cdepths[2 * ipair + 1];

        //     int offset_a = offsets_cart[2 * ipair];
        //     int offset_b = offsets_cart[2 * ipair + 1];

        //     int pos_cntr_a = poss_cntrs[2 * ipair];
        //     int pos_cntr_b = poss_cntrs[2 * ipair + 1];

        //     for (int mu = 0; mu < dim_cart_idxs_a; mu++)
        //         for (int nu = 0; nu < dim_cart_idxs_b; nu++)
        //         {
        //             int pos_a = offset_a + mu;
        //             int pos_b = offset_b + nu;
        //             int pos_ab = pos_a * dim_ao_cart + pos_b;
        //             sints_cart[pos_ab] = 0;
        //         }

        //     for (int ia = 0, iab = 0; ia < dim_cntr_a; ia++)
        //         for (int ib = 0; ib < dim_cntr_b; ib++, iab++)
        //         {
        //             double a = exps[pos_cntr_a + ia];
        //             double b = exps[pos_cntr_b + ib];
        //             double p = a + b;

        //             int pos_ecoeffs = eoffsets[id] + iab * n_ecoeffs_3x;
        //             int pos_x = pos_ecoeffs;
        //             int pos_y = pos_ecoeffs + n_ecoeffs;
        //             int pos_z = pos_ecoeffs + 2 * n_ecoeffs;

        //             double val = pow(M_PI / p, 1.5);

        //             for (int mu = 0; mu < dim_cart_idxs_a; mu++)
        //             {
        //                 int i = cart_idxs_poss_a[3 * mu];
        //                 int j = cart_idxs_poss_a[3 * mu + 1];
        //                 int k = cart_idxs_poss_a[3 * mu + 2];

        //                 for (int nu = 0; nu < dim_cart_idxs_b; nu++)
        //                 {
        //                     int i_ = cart_idxs_poss_b[3 * nu];
        //                     int j_ = cart_idxs_poss_b[3 * nu + 1];
        //                     int k_ = cart_idxs_poss_b[3 * nu + 2];

        //                     double ecoeff = ecoeffs[idxE(d3, d23, i, i_, 0, pos_x)] *
        //                                     ecoeffs[idxE(d3, d23, j, j_, 0, pos_y)] *
        //                                     ecoeffs[idxE(d3, d23, k, k_, 0, pos_z)];

        //                     int mu_ = offset_a + mu;
        //                     int nu_ = offset_b + nu;
        //                     int idx = mu_ * dim_ao_cart + nu_;
        //                     sints_cart[idx] += ccoeffs[pos_cntr_a + ia] * ccoeffs[pos_cntr_b + ib] *
        //                                        ecoeff * val;
        //                 }
        //             }
        //         }
        // }
    }
}