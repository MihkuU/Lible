#pragma once

#include <lible/ints/boys_function.hpp>
#include <lible/ints/utils.hpp>

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
        // R-ints forward-declaration.
        template <int la, int lb>
        void calcRInts(const double alpha, const double fac, const double *fnx,
                       const double *xyz_ab, double *rints_out);

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
                std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);

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

                                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab,
                                            n_sph_cd, n_hermite_cd, 1.0, &rints[0], n_hermite_cd,
                                            &ecoeffs_cd_tsp[pos_ecoeffs_cd], n_sph_cd, 1.0,
                                            &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_cd);
                            }
                    }

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
                        int pos_ecoeffs_ab = iab * n_ecoeffs_ab;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd,
                                    n_hermite_ab, 1.0, &ecoeffs_ab[pos_ecoeffs_ab], n_hermite_ab,
                                    &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_cd, 1.0,
                                    &eri4_batch[0], n_sph_cd);
                    }
            }

            template <int la, int lb, int lc>
            void eri3Kernel(const int cdepth_a, const int cdepth_b, const int cdepth_c,
                            const double *exps_a, const double *exps_b, const double *exps_c,
                            const double *coords_a, const double *coords_b, const double *coords_c,
                            const double *ecoeffs_ab, const double *ecoeffs_c, double *eri3_batch)
            {
                constexpr int lab = la + lb;
                constexpr int labc = lab + lc;

                constexpr int n_sph_a = numSphericalsC(la);
                constexpr int n_sph_b = numSphericalsC(lb);
                constexpr int n_sph_c = numSphericalsC(lc);
                constexpr int n_hermite_ab = numHermitesC(lab);
                constexpr int n_hermite_c = numHermitesC(lc);
                constexpr int n_sph_ab = n_sph_a * n_sph_b;
                constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
                constexpr int n_ecoeffs_c = n_sph_c * n_hermite_c;

                std::fill(eri3_batch, eri3_batch + n_sph_ab * n_sph_c, 0);

                std::array<double, labc + 1> fnx;
                BoysF2<labc> boys_f;

                constexpr int n_hermites_abc = numHermitesC(lab) * numHermitesC(lc);
                std::array<double, n_hermites_abc> rints;

                constexpr int n_rints_x_ecoeffs = n_sph_c * n_hermite_ab;
                std::vector<double> rints_x_ecoeffs(cdepth_a * cdepth_b * n_rints_x_ecoeffs, 0);

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
                        for (int ic = 0; ic < cdepth_c; ic++)
                        {
                            double a = exps_a[ia];
                            double b = exps_b[ib];
                            double c = exps_c[ic];

                            double p = a + b;
                            double alpha = p * c / (p + c);

                            std::array<double, 3> xyz_p{(a * coords_a[0] + b * coords_b[0]) / p,
                                                        (a * coords_a[1] + b * coords_b[1]) / p,
                                                        (a * coords_a[2] + b * coords_b[2]) / p};

                            std::array<double, 3> xyz_pc{xyz_p[0] - coords_c[0],
                                                         xyz_p[1] - coords_c[1],
                                                         xyz_p[2] - coords_c[2]};

                            double dx{xyz_pc[0]}, dy{xyz_pc[1]}, dz{xyz_pc[2]};
                            double x = alpha * (dx * dx + dy * dy + dz * dz);
                            boys_f.calcFnx(x, &fnx[0]);

                            double fac = (2.0 * std::pow(M_PI, 2.5) / (p * c * std::sqrt(p + c)));
                            calcRInts<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);

                            int pos_ecoeffs_c = ic * n_ecoeffs_c;

                            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_hermite_ab,
                                        n_sph_c, n_hermite_c, 1.0, &rints[0], n_hermite_c,
                                        &ecoeffs_c[pos_ecoeffs_c], n_hermite_c, 1.0,
                                        &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_c);
                        }
                    }

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int pos_rints_x_ecoeffs = iab * n_rints_x_ecoeffs;
                        int pos_ecoeffs_ab = iab * n_ecoeffs_ab;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                                    n_hermite_ab, 1.0, &ecoeffs_ab[pos_ecoeffs_ab], n_hermite_ab,
                                    &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_c, 1.0,
                                    &eri3_batch[0], n_sph_c);
                    }
            }

            template <int la, int lb>
            void eri2Kernel(const int cdepth_a, const int cdepth_b,
                            const double *exps_a, const double *exps_b,
                            const double *coords_a, const double *coords_b,
                            const double *ecoeffs_a, const double *ecoeffs_b_tsp,
                            double *eri2_batch)
            {
                constexpr int lab = la + lb;
                constexpr int n_sph_a = numSphericalsC(la);
                constexpr int n_sph_b = numSphericalsC(lb);
                constexpr int n_hermite_a = numHermitesC(la);
                constexpr int n_hermite_b = numHermitesC(lb);
                constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;
                constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;

                std::fill(eri2_batch, eri2_batch + n_sph_a * n_sph_b, 0);

                std::array<double, lab + 1> fnx;
                BoysF2<lab> boys_f;

                constexpr int n_hermites_ab = numHermitesC(la) * numHermitesC(lb);
                std::array<double, n_hermites_ab> rints;

                constexpr int n_rints_x_ecoeffs = n_hermite_a * n_sph_b;
                std::vector<double> rints_x_ecoeffs(cdepth_a * n_rints_x_ecoeffs, 0);

                std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                                             coords_a[1] - coords_b[1],
                                             coords_a[2] - coords_b[2]};
                double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
                double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

                for (int ia = 0; ia < cdepth_a; ia++)
                {
                    int pos_rints_x_ecoeffs = ia * n_rints_x_ecoeffs;
                    for (int ib = 0; ib < cdepth_b; ib++)
                    {
                        double a = exps_a[ia];
                        double b = exps_b[ib];

                        double alpha = a * b / (a + b);
                        double x = alpha * xyz_ab_dot;
                        boys_f.calcFnx(x, &fnx[0]);

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
                        calcRInts<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);

                        int pos_ecoeffs_b = ib * n_ecoeffs_b;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_a,
                                    n_sph_b, n_hermite_b, 1.0, &rints[0], n_hermite_b,
                                    &ecoeffs_b_tsp[pos_ecoeffs_b], n_sph_b, 1.0,
                                    &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_b);
                    }
                }

                for (int ia = 0; ia < cdepth_a; ia++)
                {
                    int pos_rints_x_ecoeffs = ia * n_rints_x_ecoeffs;
                    int pos_ecoeffs_a = ia * n_ecoeffs_a;

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b,
                                n_hermite_a, 1.0, &ecoeffs_a[pos_ecoeffs_a], n_hermite_a,
                                &rints_x_ecoeffs[pos_rints_x_ecoeffs], n_sph_b, 1.0,
                                &eri2_batch[0], n_sph_b);
                }
            }
        }
    }
}