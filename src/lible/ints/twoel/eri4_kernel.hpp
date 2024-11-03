#pragma once

#include <lible/ints/boys_function.hpp>
#include <lible/ints/rints_meta.hpp>

#include <array>
#include <vector>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace lible
{
    namespace ints
    {
        namespace two
        {
            template <int la, int lb, int lc, int ld>
            void eri4Kernel(const int cdepth_a, const int cdepth_b,
                            const int cdepth_c, const int cdepth_d,
                            const double *exps_a, const double *exps_b,
                            const double *exps_c, const double *exps_d,
                            const double *coords_a, const double *coords_b,
                            const double *coords_c, const double *coords_d,
                            const double *ecoeffs_ab,
                            const double *ecoeffs_cd_tsp,
                            double *eri4_batch)
            {
                constexpr int lab = la + lb;
                constexpr int lcd = lc + ld;
                constexpr int labcd = lab + lcd;

                constexpr int n_sph_a = numSphericalsC(la);
                constexpr int n_sph_b = numSphericalsC(lb);
                constexpr int n_sph_c = numSphericalsC(lc);
                constexpr int n_sph_d = numSphericalsC(ld);
                constexpr int n_hermite_ab = numHermitesC(lab);
                constexpr int n_hermite_cd = numHermitesC(lcd);
                constexpr int n_sph_ab = n_sph_a * n_sph_b;
                constexpr int n_sph_cd = n_sph_c * n_sph_d;
                constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
                constexpr int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;

                std::fill(eri4_batch, eri4_batch + n_sph_ab * n_sph_cd, 0);

                std::array<double, labcd + 1> fnx;
                BoysF2<labcd> boys_f;

                constexpr int n_hermites_abcd = numHermitesC(lab) * numHermitesC(lcd);
                std::array<double, n_hermites_abcd> rints;

                constexpr int n_rints_x_ecoeffs = n_sph_cd * n_hermite_ab;
                std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs);

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
                        for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                            for (int id = 0; id < cdepth_d; id++, icd++)
                            {
                                double a = exps_a[ia];
                                double b = exps_b[ib];
                                double c = exps_c[ic];
                                double d = exps_d[id];

                                double p = a + b;
                                double q = c + d;

                                std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                                            (a * coords_a[1] + b * coords_b[1]) / p,
                                                            (a * coords_a[2] + b * coords_b[2]) / p};

                                std::array<double, 3> xyz_q{(c * coords_c[0] + d * coords_d[0]) / q,
                                                            (c * coords_c[1] + d * coords_d[1]) / q,
                                                            (c * coords_c[2] + d * coords_d[2]) / q};

                                std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0], xyz_p[1] - xyz_q[1],
                                                             xyz_p[2] - xyz_q[2]};

                                double alpha = p * q / (p + q);
                                double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};
                                double x = alpha * (dx * dx + dy * dy + dz * dz);
                                boys_f.calcFnx(x, &fnx[0]);

                                double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                                calcRInts<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                                int pos_ecoeffs_cd = icd * n_ecoeffs_cd;

                                // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab,
                                //             n_sph_cd, n_hermite_cd, 1.0, &rints[0], n_hermite_cd,
                                //             &ecoeffs_cd_tsp[pos_ecoeffs_cd], n_sph_cd, 1.0,
                                //             &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_cd);
                            }
                    }

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
                        int pos_ecoeffs_ab = iab * n_ecoeffs_ab;

                        // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd,
                        //             n_hermite_ab, 1.0, &ecoeffs_ab[pos_ecoeffs_ab], n_hermite_ab,
                        //             &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_cd, 1.0,
                        //             &eri4_batch[0], n_sph_cd);
                    }
            }
        }
    }
}
