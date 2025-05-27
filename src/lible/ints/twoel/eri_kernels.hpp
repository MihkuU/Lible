#pragma once

#include <lible/ints/boys_function.hpp>
#include <lible/ints/shell_pair_data.hpp>
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
        // Forward declaration.
        template <int la, int lb>
        void calcRInts_ERI(const double alpha, const double fac, const double *fnx,
                           const double *xyz_ab, double *rints_out);

        // Forward declaration.
        template <int la, int lb>
        void calcRInts_ERI2_deriv1(const double alpha, const double fac, const double *fnx,
                                   const double *xyz_pq, double *rints_out);

        namespace two
        {
            template <int la, int lb, int lc, int ld>
            vec4d eri4Kernel(const int ipair_ab, const int ipair_cd,
                             const std::vector<double> &ecoeffs_ab,
                             const std::vector<double> &ecoeffs_cd_tsp,
                             const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd)
            {
                const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
                const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
                const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
                const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
                const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
                const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
                const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
                const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];

                const double *exps_a = &sp_data_ab.exps[cofs_a];
                const double *exps_b = &sp_data_ab.exps[cofs_b];
                const double *exps_c = &sp_data_cd.exps[cofs_c];
                const double *exps_d = &sp_data_cd.exps[cofs_d];
                const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
                const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
                const double *coords_c = &sp_data_cd.coords[6 * ipair_cd];
                const double *coords_d = &sp_data_cd.coords[6 * ipair_cd + 3];
                const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
                const double *pecoeffs_cd_tsp = &ecoeffs_cd_tsp[sp_data_cd.offsets_ecoeffs[ipair_cd]];

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

                std::array<double, labcd + 1> fnx;
                BoysF2<labcd> boys_f;

                constexpr int n_hermites_abcd = numHermitesC(lab) * numHermitesC(lcd);
                std::array<double, n_hermites_abcd> rints;

                constexpr int n_R_x_E = n_sph_cd * n_hermite_ab;
                std::vector<double> R_x_E(cdepth_a * cdepth_b * n_R_x_E, 0);

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int ofs_R_x_E = iab * n_R_x_E;
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
                                calcRInts_ERI<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                                int ofs_ecoeffs_cd = icd * n_ecoeffs_cd;

                                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_ab,
                                            n_sph_cd, n_hermite_cd, 1.0, &rints[0], n_hermite_cd,
                                            &pecoeffs_cd_tsp[ofs_ecoeffs_cd], n_sph_cd, 1.0,
                                            &R_x_E[ofs_R_x_E], n_sph_cd);
                            }
                    }

                vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int ofs_R_x_E = iab * n_R_x_E;
                        int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_cd,
                                    n_hermite_ab, 1.0, &pecoeffs_ab[ofs_ecoeffs_ab], n_hermite_ab,
                                    &R_x_E[ofs_R_x_E], n_sph_cd, 1.0, eri4_batch.memptr(),
                                    n_sph_cd);
                    }

                int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
                int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
                int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd];
                int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int nu = 0; nu < n_sph_b; nu++)
                        for (int ka = 0; ka < n_sph_c; ka++)
                            for (int ta = 0; ta < n_sph_d; ta++)
                            {
                                double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                                double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                                double norm_c = sp_data_cd.norms[ofs_norm_c + ka];
                                double norm_d = sp_data_cd.norms[ofs_norm_d + ta];

                                eri4_batch(mu, nu, ka, ta) *= norm_a * norm_b * norm_c * norm_d;
                            }

                return eri4_batch;
            }

            template <int la, int lb, int lc>
            vec3d eri3Kernel(const int ipair_ab, const int ishell_c,
                             const std::vector<double> &ecoeffs_ab,
                             const std::vector<double> &ecoeffs_c,
                             const ShellPairData &sp_data_ab, const ShellData &sh_data_c)
            {
                const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
                const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
                const int cdepth_c = sh_data_c.cdepths[ishell_c];
                const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
                const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
                const int cofs_c = sh_data_c.coffsets[ishell_c];

                const double *exps_a = &sp_data_ab.exps[cofs_a];
                const double *exps_b = &sp_data_ab.exps[cofs_b];
                const double *exps_c = &sh_data_c.exps[cofs_c];
                const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
                const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
                const double *coords_c = &sh_data_c.coords[3 * ishell_c];
                const double *pecoeffs_ab = &ecoeffs_ab[sp_data_ab.offsets_ecoeffs[ipair_ab]];
                const double *pecoeffs_c = &ecoeffs_c[sh_data_c.offsets_ecoeffs[ishell_c]];

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

                std::array<double, labc + 1> fnx;
                BoysF2<labc> boys_f;

                constexpr int n_hermites_abc = numHermitesC(lab) * numHermitesC(lc);
                std::array<double, n_hermites_abc> rints;

                constexpr int n_R_x_E = n_sph_c * n_hermite_ab;
                std::vector<double> R_x_E(cdepth_a * cdepth_b * n_R_x_E, 0);

                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int ofs_R_x_E = iab * n_R_x_E;
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
                            calcRInts_ERI<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);

                            int ofs_ecoeffs_c = ic * n_ecoeffs_c;

                            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n_hermite_ab,
                                        n_sph_c, n_hermite_c, 1.0, &rints[0], n_hermite_c,
                                        &pecoeffs_c[ofs_ecoeffs_c], n_hermite_c, 1.0,
                                        &R_x_E[ofs_R_x_E], n_sph_c);
                        }
                    }

                vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);
                for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                    for (int ib = 0; ib < cdepth_b; ib++, iab++)
                    {
                        int ofs_R_x_E = iab * n_R_x_E;
                        int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_ab, n_sph_c,
                                    n_hermite_ab, 1.0, &pecoeffs_ab[ofs_ecoeffs_ab], n_hermite_ab,
                                    &R_x_E[ofs_R_x_E], n_sph_c, 1.0, eri3_batch.memptr(), n_sph_c);
                    }

                int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab];
                int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
                int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int nu = 0; nu < n_sph_b; nu++)
                        for (int ka = 0; ka < n_sph_c; ka++)
                        {
                            double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                            double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                            double norm_c = sh_data_c.norms[ofs_norm_c + ka];

                            eri3_batch(mu, nu, ka) *= norm_a * norm_b * norm_c;
                        }

                return eri3_batch;
            }

            template <int la, int lb>
            vec2d eri2Kernel(const int ishell_a, const int ishell_b,
                             const std::vector<double> &ecoeffs_a,
                             const std::vector<double> &ecoeffs_b_tsp,
                             const ShellData &sh_data_a, const ShellData &sh_data_b)
            {
                const int cdepth_a = sh_data_a.cdepths[ishell_a];
                const int cdepth_b = sh_data_b.cdepths[ishell_b];
                const int cofs_a = sh_data_a.coffsets[ishell_a];
                const int cofs_b = sh_data_b.coffsets[ishell_b];

                const double *exps_a = &sh_data_a.exps[cofs_a];
                const double *exps_b = &sh_data_b.exps[cofs_b];
                const double *coords_a = &sh_data_a.coords[3 * ishell_a];
                const double *coords_b = &sh_data_b.coords[3 * ishell_b];
                const double *pecoeffs_a = &ecoeffs_a[sh_data_a.offsets_ecoeffs[ishell_a]];
                const double *pecoeffs_b_tsp = &ecoeffs_b_tsp[sh_data_b.offsets_ecoeffs[ishell_b]];

                constexpr int lab = la + lb;
                constexpr int n_sph_a = numSphericalsC(la);
                constexpr int n_sph_b = numSphericalsC(lb);
                constexpr int n_hermite_a = numHermitesC(la);
                constexpr int n_hermite_b = numHermitesC(lb);
                constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;
                constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;

                std::array<double, lab + 1> fnx;
                BoysF2<lab> boys_f;

                constexpr int n_hermites_ab = numHermitesC(la) * numHermitesC(lb);
                std::array<double, n_hermites_ab> rints;

                constexpr int n_R_x_E = n_hermite_a * n_sph_b;
                std::vector<double> R_x_E(cdepth_a * n_R_x_E, 0);

                std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                                             coords_a[1] - coords_b[1],
                                             coords_a[2] - coords_b[2]};
                double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
                double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

                for (int ia = 0; ia < cdepth_a; ia++)
                {
                    int ofs_R_x_E = ia * n_R_x_E;
                    for (int ib = 0; ib < cdepth_b; ib++)
                    {
                        double a = exps_a[ia];
                        double b = exps_b[ib];

                        double alpha = a * b / (a + b);
                        double x = alpha * xyz_ab_dot;
                        boys_f.calcFnx(x, &fnx[0]);

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
                        calcRInts_ERI<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);

                        int ofs_ecoeffs_b = ib * n_ecoeffs_b;

                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_hermite_a,
                                    n_sph_b, n_hermite_b, 1.0, &rints[0], n_hermite_b,
                                    &pecoeffs_b_tsp[ofs_ecoeffs_b], n_sph_b, 1.0,
                                    &R_x_E[ofs_R_x_E], n_sph_b);
                    }
                }

                vec2d eri2_batch(Fill(0), n_sph_a, n_sph_b);
                for (int ia = 0; ia < cdepth_a; ia++)
                {
                    int ofs_R_x_E = ia * n_R_x_E;
                    int pos_ecoeffs_a = ia * n_ecoeffs_a;

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b,
                                n_hermite_a, 1.0, &pecoeffs_a[pos_ecoeffs_a], n_hermite_a,
                                &R_x_E[ofs_R_x_E], n_sph_b, 1.0, eri2_batch.memptr(), n_sph_b);
                }

                int ofs_norm_a = sh_data_a.offsets_norms[ishell_a];
                int ofs_norm_b = sh_data_b.offsets_norms[ishell_b];
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int nu = 0; nu < n_sph_b; nu++)
                    {
                        double norm_a = sh_data_a.norms[ofs_norm_a + mu];
                        double norm_b = sh_data_b.norms[ofs_norm_b + nu];

                        eri2_batch(mu, nu) *= norm_a * norm_b;
                    }

                return eri2_batch;
            }

            template <int la, int lb>
            std::array<vec2d, 6> eri2d1Kernel(const int ishell_a, const int ishell_b,
                                              const std::vector<double> &ecoeffs_a,
                                              const std::vector<double> &ecoeffs_b_tsp,
                                              const ShellData &sh_data_a, const ShellData &sh_data_b)
            {
                const int cdepth_a = sh_data_a.cdepths[ishell_a];
                const int cdepth_b = sh_data_b.cdepths[ishell_b];
                const int cofs_a = sh_data_a.coffsets[ishell_a];
                const int cofs_b = sh_data_b.coffsets[ishell_b];

                const double *exps_a = &sh_data_a.exps[cofs_a];
                const double *exps_b = &sh_data_b.exps[cofs_b];
                const double *coords_a = &sh_data_a.coords[3 * ishell_a];
                const double *coords_b = &sh_data_b.coords[3 * ishell_b];
                const double *pecoeffs_a = &ecoeffs_a[sh_data_a.offsets_ecoeffs[ishell_a]];
                const double *pecoeffs_b_tsp = &ecoeffs_b_tsp[sh_data_b.offsets_ecoeffs[ishell_b]];                

                // vector<array<int, 3>> idxs_tuv_a = returnHermiteGaussianIdxs(la);
                // vector<array<int, 3>> idxs_tuv_b = returnHermiteGaussianIdxs(lb);

                constexpr int lab = la + lb;
                constexpr int n_sph_a = numSphericalsC(la);
                constexpr int n_sph_b = numSphericalsC(lb);
                constexpr int n_hermite_a = numHermitesC(la);
                constexpr int n_hermite_b = numHermitesC(lb);
                constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;
                constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;

                std::array<double, lab + 2> fnx;
                BoysF2<lab + 1> boys_f;

                constexpr int n_hermites_ab = numHermitesC(la) * numHermitesC(lb);
                std::array<double, 6 * n_hermites_ab> rints;

                std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                                             coords_a[1] - coords_b[1],
                                             coords_a[2] - coords_b[2]};
                double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
                double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

                constexpr int n_R_x_E = n_hermite_a * n_sph_b;
                std::vector<double> R_x_E(6 * cdepth_a * n_R_x_E, 0);
                for (int ia = 0; ia < cdepth_a; ia++)
                {
                    int ofs_R_x_E = 6 * ia * n_R_x_E;
                    for (int ib = 0; ib < cdepth_b; ib++)
                    {
                        double a = exps_a[ia];
                        double b = exps_b[ib];

                        double alpha = a * b / (a + b);
                        double x = alpha * xyz_ab_dot;
                        boys_f.calcFnx(x, &fnx[0]);                                                

                        double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
                        calcRInts_ERI2_deriv1<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);

                        int ofs_ecoeffs_b = ib * n_ecoeffs_b;
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6 * n_hermite_a,
                                    n_sph_b, n_hermite_b, 1.0, &rints[0], n_hermite_b,
                                    &pecoeffs_b_tsp[ofs_ecoeffs_b], n_sph_b, 1.0,
                                    &R_x_E[ofs_R_x_E], n_sph_b);
                    }
                }

                std::array<vec2d, 6> eri2_batch;
                for (int ideriv = 0; ideriv < 6; ideriv++)
                    eri2_batch[ideriv] = vec2d(Fill(0), n_sph_a, n_sph_b);

                for (int ia = 0; ia < cdepth_a; ia++)
                {
                    int start = ia * 6 * n_R_x_E;
                    int ofs0 = start + 0 * n_R_x_E;
                    int ofs1 = start + 1 * n_R_x_E;
                    int ofs2 = start + 2 * n_R_x_E;
                    int ofs3 = start + 3 * n_R_x_E;
                    int ofs4 = start + 4 * n_R_x_E;
                    int ofs5 = start + 5 * n_R_x_E;

                    int ofs_ecoeffs = ia * n_ecoeffs_a;

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                                1.0, &pecoeffs_a[ofs_ecoeffs], n_hermite_a, &R_x_E[ofs0], n_sph_b, 1.0,
                                eri2_batch[0].memptr(), n_sph_b);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                                1.0, &pecoeffs_a[ofs_ecoeffs], n_hermite_a, &R_x_E[ofs1], n_sph_b, 1.0,
                                eri2_batch[1].memptr(), n_sph_b);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                                1.0, &pecoeffs_a[ofs_ecoeffs], n_hermite_a, &R_x_E[ofs2], n_sph_b, 1.0,
                                eri2_batch[2].memptr(), n_sph_b);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                                1.0, &pecoeffs_a[ofs_ecoeffs], n_hermite_a, &R_x_E[ofs3], n_sph_b, 1.0,
                                eri2_batch[3].memptr(), n_sph_b);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                                1.0, &pecoeffs_a[ofs_ecoeffs], n_hermite_a, &R_x_E[ofs4], n_sph_b, 1.0,
                                eri2_batch[4].memptr(), n_sph_b);

                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_sph_a, n_sph_b, n_hermite_a,
                                1.0, &pecoeffs_a[ofs_ecoeffs], n_hermite_a, &R_x_E[ofs5], n_sph_b, 1.0,
                                eri2_batch[5].memptr(), n_sph_b);
                }

                int ofs_norm_a = sh_data_a.offsets_norms[ishell_a];
                int ofs_norm_b = sh_data_b.offsets_norms[ishell_b];

                const double *norms_a = &sh_data_a.norms[ofs_norm_a];
                const double *norms_b = &sh_data_b.norms[ofs_norm_b];

                for (int ideriv = 0; ideriv < 6; ideriv++)
                    for (int mu = 0; mu < n_sph_a; mu++)
                        for (int nu = 0; nu < n_sph_b; nu++)
                        {
                            double norm_a = norms_a[mu];
                            double norm_b = norms_b[nu];

                            eri2_batch[ideriv](mu, nu) *= norm_a * norm_b;
                        }

                return eri2_batch;
            }
        }
    }
}