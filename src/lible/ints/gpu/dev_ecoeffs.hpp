#pragma once

#include <lible/ints/gpu/dev_utils.hpp>

#include <hip/hip_runtime.h>

namespace lible::ints::gpu
{
    __device__ inline void calcECoeffsX(const int la, const int lb, const int pos, const double PA,
                                        const double PB, const double one_o_2p, double *ecoeffs)
    {
        int d2 = lb + 1;
        int d3 = la + lb + 1;
        int d23 = d2 * d3;

        ecoeffs[idxE(d3, d23, 0, 0, 0, pos)] = 0;

        for (int i = 1; i <= la; i++)
        {
            ecoeffs[idxE(d3, d23, i, 0, 0, pos)] = PA * ecoeffs[idxE(d3, d23, i - 1, 0, 0, pos)] +
                                                   ecoeffs[idxE(d3, d23, i - 1, 0, 1, pos)];

            for (int t = 1; t < i; t++)
                ecoeffs[idxE(d3, d23, i, 0, t, pos)] = one_o_2p * ecoeffs[idxE(d3, d23, i - 1, 0, t - 1, pos)] +
                                                       PA * ecoeffs[idxE(d3, d23, i - 1, 0, t, pos)] +
                                                       (t + 1) * ecoeffs[idxE(d3, d23, i - 1, 0, t + 1, pos)];

            ecoeffs[idxE(d3, d23, i, 0, i, pos)] = one_o_2p * ecoeffs[idxE(d3, d23, i - 1, 0, i - 1, pos)] +
                                                   PA * ecoeffs[idxE(d3, d23, i - 1, 0, i, pos)];
        }

        for (int j = 1; j <= lb; j++)
            for (int i = 0; i <= la; i++)
            {
                ecoeffs[idxE(d3, d23, i, j, 0, pos)] = PB * ecoeffs[idxE(d3, d23, i, j - 1, 0, pos)] +
                                                       ecoeffs[idxE(d3, d23, i, j - 1, 1, pos)];

                for (int t = 1; t < i + j; t++)
                    ecoeffs[idxE(d3, d23, i, j, t, pos)] = one_o_2p * ecoeffs[idxE(d3, d23, i, j - 1, t - 1, pos)] +
                                                           PB * ecoeffs[idxE(d3, d23, i, j - 1, t, pos)] +
                                                           (t + 1) * ecoeffs[idxE(d3, d23, i, j - 1, t + 1, pos)];

                ecoeffs[idxE(d3, d23, i, j, i + j, pos)] = one_o_2p * ecoeffs[idxE(d3, d23, i, j - 1, i + j - 1, pos)] +
                                                           PB * ecoeffs[idxE(d3, d23, i, j - 1, i + j, pos)];
            }
    }

    __global__ void calcECoeffs(const int la, const int lb, const int n_shell_pairs,
                                const int *cdepths, const int *eoffsets,
                                const int *poss_cntrs, const double *coords,
                                const double *exps, double *ecoeffs)
    {
        int id = blockDim.x * blockIdx.x + threadIdx.x;
        if (id < n_shell_pairs)
        {
            int ipair = id;

            int n_ecoeffs = (la + 1) * (lb + 1) * (la + lb + 1);
            int n_ecoeffs_3x = 3 * n_ecoeffs;

            int dim_cntr_a = cdepths[2 * ipair];
            int dim_cntr_b = cdepths[2 * ipair + 1];
            int pos_cntr_a = poss_cntrs[2 * ipair];
            int pos_cntr_b = poss_cntrs[2 * ipair + 1];

            double xa = coords[6 * ipair];
            double ya = coords[6 * ipair + 1];
            double za = coords[6 * ipair + 2];
            double xb = coords[6 * ipair + 3];
            double yb = coords[6 * ipair + 4];
            double zb = coords[6 * ipair + 5];
            double rab2[3] = {std::pow(xa - ya, 2), std::pow(xb - yb, 2), std::pow(za - zb, 2)};

            for (int ia = 0, iab = 0; ia < dim_cntr_a; ia++)
                for (int ib = 0; ib < dim_cntr_b; ib++, iab++)
                {
                    double a = exps[pos_cntr_a + ia];
                    double b = exps[pos_cntr_b + ib];

                    double p = a + b;
                    double one_o_2p = 1.0 / (2 * p);

                    double P[3] = {(a * xa + b * xb) / p, (a * ya + b * yb) / p, (a * za + b * zb) / p};
                    double PA[3] = {P[0] - xa, P[1] - ya, P[2] - za};
                    double PB[3] = {P[0] - xb, P[1] - yb, P[2] - zb};

                    double mu = a * b / (a + b);
                    double Kab[3] = {std::exp(-mu * rab2[0]), std::exp(-mu * rab2[1]), std::exp(-mu * rab2[2])};

                    int pos_ecoeffs = eoffsets[id] + iab * n_ecoeffs_3x;
                    int pos_x = pos_ecoeffs;
                    int pos_y = pos_ecoeffs + n_ecoeffs;
                    int pos_z = pos_ecoeffs + 2 * n_ecoeffs;

                    // calcECoeffsX(la, lb, pos_x, PA[0], PB[0], one_o_2p, ecoeffs);
                    // calcECoeffsX(la, lb, pos_y, PA[1], PB[1], one_o_2p, ecoeffs);
                    // calcECoeffsX(la, lb, pos_z, PA[2], PB[2], one_o_2p, ecoeffs);
                }
        }
    }
}