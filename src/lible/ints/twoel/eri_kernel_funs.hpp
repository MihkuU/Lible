#pragma once

#include <lible/ints/boys_function.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>
#include <lible/ints/twoel/shark_mm_kernels.hpp>

#include <array>
#include <cstring> 
#include <functional>
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
        template <int la, int lb>
        void calcRInts_ERI(const double alpha, const double fac, const double *fnx,
                           const double *xyz_ab, double *rints_out);

        template <int la, int lb>
        void calcRInts_ERI2D1(const double alpha, const double fac, const double *fnx,
                              const double *xyz_ab, double *rints);

        template <int lab, int lc>
        void calcRInts_ERI3D1(const double alpha, const double fac, const double *fnx,
                              const double *xyz_pc, double *rints);

        template <int lab, int lcd>
        void calcRInts_ERI4D1(const double alpha, const double fac, const double *fnx,
                              const double *xyz_pq, const int n_rints, const int ofs_row,
                              const int ofs_col, const int n_cols, const int n_rows,
                              double *rints);

        // Forward decls.
        struct ERI4Kernel;
        struct ERI3Kernel;
        struct ERI2Kernel;

        template <int la, int lb, int lc, int ld>
        vec4d eri4KernelFun(const int ipair_ab, const int ipair_cd,
                            const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                            const ERI4Kernel *eri4_kernel)
        {
            // Compile-time data
            constexpr int lab = la + lb;
            constexpr int lcd = lc + ld;
            constexpr int labcd = lab + lcd;

            constexpr int n_sph_a = numSphericalsC(la);
            constexpr int n_sph_b = numSphericalsC(lb);
            constexpr int n_sph_c = numSphericalsC(lc);
            constexpr int n_sph_d = numSphericalsC(ld);
            constexpr int n_sph_ab = n_sph_a * n_sph_b;
            constexpr int n_sph_cd = n_sph_c * n_sph_d;
            constexpr int n_hermite_ab = numHermitesC(lab);
            constexpr int n_hermite_cd = numHermitesC(lcd);
            constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
            constexpr int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;

            // Read-in data
            const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
            const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
            const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd];
            const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
            const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
            const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
            const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd];
            const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
            const int ofs_E_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
            const int ofs_E_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];

            const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab];
            const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
            const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd];
            const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
            const double *exps_a = &sp_data_ab.exps[cofs_a];
            const double *exps_b = &sp_data_ab.exps[cofs_b];
            const double *exps_c = &sp_data_cd.exps[cofs_c];
            const double *exps_d = &sp_data_cd.exps[cofs_d];
            const double *ecoeffs_ab = &eri4_kernel->ecoeffs_bra[ofs_E_ab];
            const double *ecoeffs_cd = &eri4_kernel->ecoeffs_ket[ofs_E_cd];

            // SHARK integrals
            std::array<double, labcd + 1> fnx;
            BoysF2<labcd> boys_f;

            vec4d eri4_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);
            std::array<double, n_hermite_ab * n_hermite_cd> rints;
            for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (int ib = 0; ib < cdepth_b; ib++, iab++)
                {
                    std::array<double, n_hermite_ab * n_sph_cd> R_x_E{};
                    for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                        for (int id = 0; id < cdepth_d; id++, icd++)
                        {
                            double a = exps_a[ia];
                            double b = exps_b[ib];
                            double c = exps_c[ic];
                            double d = exps_d[id];

                            double p = a + b;
                            double q = c + d;
                            double alpha = p * q / (p + q);

                            std::array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                                        (a * xyz_a[1] + b * xyz_b[1]) / p,
                                                        (a * xyz_a[2] + b * xyz_b[2]) / p};

                            std::array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                                        (c * xyz_c[1] + d * xyz_d[1]) / q,
                                                        (c * xyz_c[2] + d * xyz_d[2]) / q};

                            std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                                         xyz_p[1] - xyz_q[1],
                                                         xyz_p[2] - xyz_q[2]};

                            double dx{xyz_pq[0]}, dy{xyz_pq[1]}, dz{xyz_pq[2]};
                            double x = alpha * (dx * dx + dy * dy + dz * dz);
                            boys_f.calcFnx(x, &fnx[0]);

                            double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                            calcRInts_ERI<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                            int ofs_ecoeffs_cd = icd * n_ecoeffs_cd;
                            shark_mm_ket2<lab, lc, ld>(&rints[0], &ecoeffs_cd[ofs_ecoeffs_cd], &R_x_E[0]);
                        }
                    int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs_ab[ofs_ecoeffs_ab], &R_x_E[0], &eri4_batch[0]);
                }

            // Norms
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

        vec4d eri4KernelFun(const int ipair_ab, const int ipair_cd,
                            const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                            const ERI4Kernel *eri4_kernel);

        template <int la, int lb, int lc>
        vec3d eri3KernelFun(const int ipair_ab, const int ishell_c,
                            const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                            const ERI3Kernel *eri3_kernel)
        {
            // Compile-time data
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

            // Read-in data
            const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
            const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
            const int cdepth_c = sh_data_c.cdepths[ishell_c];
            const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
            const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
            const int cofs_c = sh_data_c.coffsets[ishell_c];
            const int ofs_E_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
            const int ofs_E_c = sh_data_c.offsets_ecoeffs[ishell_c];

            const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
            const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
            const double *coords_c = &sh_data_c.coords[3 * ishell_c];
            const double *exps_a = &sp_data_ab.exps[cofs_a];
            const double *exps_b = &sp_data_ab.exps[cofs_b];
            const double *exps_c = &sh_data_c.exps[cofs_c];
            const double *ecoeffs_ab = &eri3_kernel->ecoeffs_bra[ofs_E_ab];
            const double *ecoeffs_c = &eri3_kernel->ecoeffs_ket[ofs_E_c];

            // SHARK integrals
            std::array<double, labc + 1> fnx;
            BoysF2<labc> boys_f;

            vec3d eri3_batch(Fill(0), n_sph_a, n_sph_b, n_sph_c);
            std::array<double, n_hermite_ab * n_hermite_c> rints;
            for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (int ib = 0; ib < cdepth_b; ib++, iab++)
                {
                    std::array<double, n_hermite_ab * n_sph_c> R_x_E{};
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
                        shark_mm_ket1<lab, lc>(&rints[0], &ecoeffs_c[ofs_ecoeffs_c], &R_x_E[0]);
                    }
                    int ofs_ecoeffs_ab = iab * n_ecoeffs_ab;
                    shark_mm_bra2<la, lb, lc>(&ecoeffs_ab[ofs_ecoeffs_ab], &R_x_E[0], &eri3_batch[0]);
                }

            // Norms
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

        vec3d eri3KernelFun(const int ipair_ab, const int ishell_c,
                            const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                            const ERI3Kernel *eri3_kernel);

        template <int la, int lb>
        vec2d eri2KernelFun(const int ishell_a, const int ishell_b,
                            const ShellData &sh_data_a, const ShellData &sh_data_b,
                            const ERI2Kernel *eri2_kernel)
        {
            // Compile-time data
            constexpr int lab = la + lb;

            constexpr int n_sph_a = numSphericalsC(la);
            constexpr int n_sph_b = numSphericalsC(lb);
            constexpr int n_hermite_a = numHermitesC(la);
            constexpr int n_hermite_b = numHermitesC(lb);
            constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;
            constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;

            // Read-in data
            const int cdepth_a = sh_data_a.cdepths[ishell_a];
            const int cdepth_b = sh_data_b.cdepths[ishell_b];
            const int cofs_a = sh_data_a.coffsets[ishell_a];
            const int cofs_b = sh_data_b.coffsets[ishell_b];
            const int ofs_E_a = sh_data_a.offsets_ecoeffs[ishell_a];
            const int ofs_E_b = sh_data_b.offsets_ecoeffs[ishell_b];

            const double *coords_a = &sh_data_a.coords[3 * ishell_a];
            const double *coords_b = &sh_data_b.coords[3 * ishell_b];
            const double *exps_a = &sh_data_a.exps[cofs_a];
            const double *exps_b = &sh_data_b.exps[cofs_b];
            const double *ecoeffs_a = &eri2_kernel->ecoeffs_bra[ofs_E_a];
            const double *ecoeffs_b = &eri2_kernel->ecoeffs_ket[ofs_E_b];

            // SHARK integrals
            std::array<double, lab + 1> fnx;
            BoysF2<lab> boys_f;

            std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                                         coords_a[1] - coords_b[1],
                                         coords_a[2] - coords_b[2]};
            double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
            double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

            vec2d eri2_batch(Fill(0), n_sph_a, n_sph_b);
            std::array<double, n_hermite_a * n_hermite_b> rints;
            for (int ia = 0; ia < cdepth_a; ia++)
            {
                std::array<double, n_hermite_a * n_sph_b> R_x_E{};
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
                    shark_mm_ket1<la, lb>(&rints[0], &ecoeffs_b[ofs_ecoeffs_b], &R_x_E[0]);
                }
                int ofs_ecoeffs_a = ia * n_ecoeffs_a;
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[0], &eri2_batch[0]);
            }

            // Norms
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

        vec2d eri2KernelFun(const int ishell_a, const int ishell_b,
                            const ShellData &sh_data_a, const ShellData &sh_data_b,
                            const ERI2Kernel *eri2_kernel);

        template <int la, int lb>
        std::array<vec2d, 6> eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                             const ShellData &sh_data_a,
                                             const ShellData &sh_data_b,
                                             const ERI2D1Kernel *eri2d1_kernel)
        {
            // Compile-time data
            constexpr int lab = la + lb;

            constexpr int n_sph_a = numSphericalsC(la);
            constexpr int n_sph_b = numSphericalsC(lb);
            constexpr int n_hermite_a = numHermitesC(la);
            constexpr int n_hermite_b = numHermitesC(lb);
            constexpr int n_rints = n_hermite_a * n_hermite_b;
            constexpr int n_ecoeffs_a = n_sph_a * n_hermite_a;
            constexpr int n_ecoeffs_b = n_sph_b * n_hermite_b;
            constexpr int n_R_x_E = n_hermite_a * n_sph_b;

            // Read-in data
            const int cdepth_a = sh_data_a.cdepths[ishell_a];
            const int cdepth_b = sh_data_b.cdepths[ishell_b];
            const int cofs_a = sh_data_a.coffsets[ishell_a];
            const int cofs_b = sh_data_b.coffsets[ishell_b];
            const int ofs_E_a = sh_data_a.offsets_ecoeffs[ishell_a];
            const int ofs_E_b = sh_data_b.offsets_ecoeffs[ishell_b];

            const double *exps_a = &sh_data_a.exps[cofs_a];
            const double *exps_b = &sh_data_b.exps[cofs_b];
            const double *coords_a = &sh_data_a.coords[3 * ishell_a];
            const double *coords_b = &sh_data_b.coords[3 * ishell_b];
            const double *ecoeffs_a = &eri2d1_kernel->ecoeffs_bra[ofs_E_a];
            const double *ecoeffs_b = &eri2d1_kernel->ecoeffs_ket[ofs_E_b];

            // SHARK integrals
            std::array<double, lab + 2> fnx;
            BoysF2<lab + 1> boys_f;

            std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                                         coords_a[1] - coords_b[1],
                                         coords_a[2] - coords_b[2]};
            double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
            double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

            std::array<vec2d, 6> eri2_batch;
            for (int ideriv = 0; ideriv < 6; ideriv++)
                eri2_batch[ideriv] = vec2d(Fill(0), n_sph_a, n_sph_b);

            std::array<double, 3 * n_rints> rints{};
            for (int ia = 0; ia < cdepth_a; ia++)
            {                
                std::array<double, 6 * n_R_x_E> R_x_E{};
                for (int ib = 0; ib < cdepth_b; ib++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];

                    double alpha = a * b / (a + b);
                    double x = alpha * xyz_ab_dot;
                    boys_f.calcFnx(x, &fnx[0]);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
                    calcRInts_ERI2D1<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);  

                    std::array<double, 3 * n_R_x_E> I{};

                    int ofs_ecoeffs_b = ib * n_ecoeffs_b;
                    shark_mm_ket1<la, lb>(&rints[0 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &I[0 * n_R_x_E]);
                    shark_mm_ket1<la, lb>(&rints[1 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &I[1 * n_R_x_E]);
                    shark_mm_ket1<la, lb>(&rints[2 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &I[2 * n_R_x_E]);

                    cblas_daxpy(3 * n_R_x_E, 1.0, &I[0], 1, &R_x_E[0 * n_R_x_E], 1);
                    cblas_daxpy(3 * n_R_x_E, -1.0, &I[0], 1, &R_x_E[3 * n_R_x_E], 1);
                }
                int ofs_ecoeffs_a = ia * n_ecoeffs_a;
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[0 * n_R_x_E], &eri2_batch[0][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[1 * n_R_x_E], &eri2_batch[1][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[2 * n_R_x_E], &eri2_batch[2][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[3 * n_R_x_E], &eri2_batch[3][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[4 * n_R_x_E], &eri2_batch[4][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[5 * n_R_x_E], &eri2_batch[5][0]);
            }

            int ofs_norm_a = sh_data_a.offsets_norms[ishell_a];
            int ofs_norm_b = sh_data_b.offsets_norms[ishell_b];
            for (int ideriv = 0; ideriv < 6; ideriv++)
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int nu = 0; nu < n_sph_b; nu++)
                    {
                        double norm_a = sh_data_a.norms[ofs_norm_a + mu];
                        double norm_b = sh_data_b.norms[ofs_norm_b + nu];

                        eri2_batch[ideriv](mu, nu) *= norm_a * norm_b;
                    }

            return eri2_batch;
        }

        std::array<vec2d, 6> eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                             const ShellData &sh_data_a,
                                             const ShellData &sh_data_b,
                                             const ERI2D1Kernel *eri2d1_kernel);

        template <int la, int lb, int lc>
        std::array<vec3d, 9> eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                             const ShellPairData &sp_data_ab,
                                             const ShellData &sh_data_c,
                                             const ERI3D1Kernel *eri3d1_kernel)
        {
            // Compile-time data
            constexpr int lab = la + lb;
            constexpr int labc = lab + lc;

            constexpr int n_sph_a = numSphericalsC(la);
            constexpr int n_sph_b = numSphericalsC(lb);
            constexpr int n_sph_c = numSphericalsC(lc);
            constexpr int n_hermite_ab = numHermitesC(lab);
            constexpr int n_hermite_c = numHermitesC(lc);
            constexpr int n_sph_ab = n_sph_a * n_sph_b;
            constexpr int n_sph_abc = n_sph_a * n_sph_b * n_sph_c;
            constexpr int n_rints = n_hermite_ab * n_hermite_c;
            constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
            constexpr int n_ecoeffs_c = n_sph_c * n_hermite_c;
            constexpr int n_R_x_E = n_hermite_ab * n_sph_c;

            // Read-in data
            const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab];
            const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
            const int cdepth_c = sh_data_c.cdepths[ishell_c];
            const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab];
            const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
            const int cofs_c = sh_data_c.coffsets[ishell_c];
            const int ofs_E0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
            const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
            const int ofs_E0_c = sh_data_c.offsets_ecoeffs[ishell_c];

            const double *exps_a = &sp_data_ab.exps[cofs_a];
            const double *exps_b = &sp_data_ab.exps[cofs_b];
            const double *exps_c = &sh_data_c.exps[cofs_c];
            const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
            const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
            const double *coords_c = &sh_data_c.coords[3 * ishell_c];
            const double *ecoeffs0_ab = &eri3d1_kernel->ecoeffs0_bra[ofs_E0_ab];
            const double *ecoeffs1_ab = &eri3d1_kernel->ecoeffs1_bra[ofs_E1_ab];
            const double *ecoeffs0_c = &eri3d1_kernel->ecoeffs0_ket[ofs_E0_c];

            // SHARK integrals
            std::array<double, labc + 2> fnx;
            BoysF2<labc + 1> boys_f;

            std::array<vec3d, 9> eri3_batch;
            for (int ideriv = 0; ideriv < 9; ideriv++)
                eri3_batch[ideriv] = vec3d(Fill(0), n_sph_a, n_sph_b, n_sph_c);

            std::array<double, 4 * n_rints> rints;
            for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (int ib = 0; ib < cdepth_b; ib++, iab++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double p = a + b;

                    std::array<double, 7 * n_R_x_E> R_x_E{};
                    for (int ic = 0; ic < cdepth_c; ic++)
                    {
                        double c = exps_c[ic];

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
                        calcRInts_ERI3D1<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);

                        std::array<double, 3 * n_R_x_E> I{};

                        int ofs_ecoeffs_c = ic * n_ecoeffs_c;
                        shark_mm_ket1<lab, lc>(&rints[0 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &I[0 * n_R_x_E]);
                        shark_mm_ket1<lab, lc>(&rints[1 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &I[1 * n_R_x_E]);
                        shark_mm_ket1<lab, lc>(&rints[2 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &I[2 * n_R_x_E]);

                        shark_mm_ket1<lab, lc>(&rints[3 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &R_x_E[3 * n_R_x_E]);

                        cblas_daxpy(3 * n_R_x_E, 1.0, &I[0], 1, &R_x_E[0 * n_R_x_E], 1);
                        cblas_daxpy(3 * n_R_x_E, -1.0, &I[0], 1, &R_x_E[4 * n_R_x_E], 1);
                    }

                    int ofs_ecoeffs0_ab = iab * n_ecoeffs_ab;
                    int ofs_ecoeffs1_ab = 3 * iab * n_ecoeffs_ab;

                    // P & R
                    std::array<double, 3 * n_sph_abc> P{};
                    std::array<double, 3 * n_sph_abc> R{};

                    shark_mm_bra2<la, lb, lc>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[0 * n_R_x_E], &P[0 * n_sph_abc]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[1 * n_R_x_E], &P[1 * n_sph_abc]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[2 * n_R_x_E], &P[2 * n_sph_abc]);                    

                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_ecoeffs1_ab + 0 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[0 * n_sph_abc]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_ecoeffs1_ab + 1 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[1 * n_sph_abc]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_ecoeffs1_ab + 2 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[2 * n_sph_abc]);

                    // A
                    cblas_daxpy(n_sph_abc, (a / p), &P[0 * n_sph_abc], 1, &eri3_batch[0][0], 1);
                    cblas_daxpy(n_sph_abc, (a / p), &P[1 * n_sph_abc], 1, &eri3_batch[1][0], 1);
                    cblas_daxpy(n_sph_abc, (a / p), &P[2 * n_sph_abc], 1, &eri3_batch[2][0], 1);

                    cblas_daxpy(n_sph_abc, 1.0, &R[0 * n_sph_abc], 1, &eri3_batch[0][0], 1);
                    cblas_daxpy(n_sph_abc, 1.0, &R[1 * n_sph_abc], 1, &eri3_batch[1][0], 1);
                    cblas_daxpy(n_sph_abc, 1.0, &R[2 * n_sph_abc], 1, &eri3_batch[2][0], 1);

                    // B
                    cblas_daxpy(n_sph_abc, (b / p), &P[0 * n_sph_abc], 1, &eri3_batch[3][0], 1);
                    cblas_daxpy(n_sph_abc, (b / p), &P[1 * n_sph_abc], 1, &eri3_batch[4][0], 1);
                    cblas_daxpy(n_sph_abc, (b / p), &P[2 * n_sph_abc], 1, &eri3_batch[5][0], 1);

                    cblas_daxpy(n_sph_abc, -1.0, &R[0 * n_sph_abc], 1, &eri3_batch[3][0], 1);
                    cblas_daxpy(n_sph_abc, -1.0, &R[1 * n_sph_abc], 1, &eri3_batch[4][0], 1);
                    cblas_daxpy(n_sph_abc, -1.0, &R[2 * n_sph_abc], 1, &eri3_batch[5][0], 1);

                    // C
                    shark_mm_bra2<la, lb, lc>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[4 * n_R_x_E], &eri3_batch[6][0]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[5 * n_R_x_E], &eri3_batch[7][0]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[6 * n_R_x_E], &eri3_batch[8][0]);
                }

            // Norms
            int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab + 0];
            int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
            int ofs_norm_c = sh_data_c.offsets_norms[ishell_c];
            for (int ideriv = 0; ideriv < 9; ideriv++)
                for (int a = 0; a < n_sph_a; a++)
                    for (int b = 0; b < n_sph_b; b++)
                        for (int c = 0; c < n_sph_c; c++)
                        {
                            double norm_a = sp_data_ab.norms[ofs_norm_a + a];
                            double norm_b = sp_data_ab.norms[ofs_norm_b + b];
                            double norm_c = sh_data_c.norms[ofs_norm_c + c];

                            eri3_batch[ideriv](a, b, c) *= norm_a * norm_b * norm_c;
                        }

            return eri3_batch;
        }

        std::array<vec3d, 9> eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                             const ShellPairData &sp_data_ab,
                                             const ShellData &sh_data_c,
                                             const ERI3D1Kernel *eri3d1_kernel);

        template <int la, int lb, int lc, int ld>
        std::array<vec4d, 12> eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                              const ShellPairData &sp_data_ab,
                                              const ShellPairData &sp_data_cd,
                                              const ERI4D1Kernel *eri4d1_kernel)
        {
            // Compile-time data
            constexpr int lab = la + lb;
            constexpr int lcd = lc + ld;
            constexpr int labcd = lab + lcd;

            constexpr int n_sph_a = numSphericalsC(la);
            constexpr int n_sph_b = numSphericalsC(lb);
            constexpr int n_sph_c = numSphericalsC(lc);
            constexpr int n_sph_d = numSphericalsC(ld);
            constexpr int n_sph_ab = n_sph_a * n_sph_b;
            constexpr int n_sph_cd = n_sph_c * n_sph_d;
            constexpr int n_sph_abcd = n_sph_ab * n_sph_cd;
            constexpr int n_hermite_ab = numHermitesC(lab);
            constexpr int n_hermite_cd = numHermitesC(lcd);

            // Read-in data
            const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
            const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
            const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
            const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
            const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
            const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
            const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
            const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
            const int ofs_E0_bra = sp_data_ab.offsets_ecoeffs[ipair_ab];
            const int ofs_E1_bra = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
            const int ofs_E0_ket = sp_data_cd.offsets_ecoeffs[ipair_cd];
            const int ofs_E1_ket = sp_data_cd.offsets_ecoeffs_deriv1[ipair_cd];

            const double *exps_a = &sp_data_ab.exps[cofs_a];
            const double *exps_b = &sp_data_ab.exps[cofs_b];
            const double *exps_c = &sp_data_cd.exps[cofs_c];
            const double *exps_d = &sp_data_cd.exps[cofs_d];
            const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
            const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
            const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
            const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
            const double *ecoeffs0_bra = &eri4d1_kernel->ecoeffs0_bra[ofs_E0_bra];
            const double *ecoeffs1_bra = &eri4d1_kernel->ecoeffs1_bra[ofs_E1_bra];
            const double *ecoeffs0_ket = &eri4d1_kernel->ecoeffs0_ket[ofs_E0_ket];
            const double *ecoeffs1_ket = &eri4d1_kernel->ecoeffs1_ket[ofs_E1_ket];

            // R-integrals
            std::array<double, labcd + 2> fnx;
            BoysF2<labcd + 1> boys_f;

            int n_rrows = (cdepth_a * cdepth_b * n_hermite_ab);
            int n_rcols = (cdepth_c * cdepth_d * n_hermite_cd);
            int n_rints = n_rrows * n_rcols;
            std::vector<double> ecoeffs0_bra_ap(n_sph_ab * n_rrows);
            std::vector<double> ecoeffs0_ket_cq(n_sph_cd * n_rcols);
            std::vector<double> rints(8 * n_rints);
            for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (int ib = 0; ib < cdepth_b; ib++, iab++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double p = a + b;
                    for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                        for (int id = 0; id < cdepth_d; id++, icd++)
                        {
                            double c = exps_c[ic];
                            double d = exps_d[id];

                            double q = c + d;
                            double alpha = p * q / (p + q);

                            std::array<double, 3> xyz_p{(a * xyz_a[0] + b * xyz_b[0]) / p,
                                                        (a * xyz_a[1] + b * xyz_b[1]) / p,
                                                        (a * xyz_a[2] + b * xyz_b[2]) / p};

                            std::array<double, 3> xyz_q{(c * xyz_c[0] + d * xyz_d[0]) / q,
                                                        (c * xyz_c[1] + d * xyz_d[1]) / q,
                                                        (c * xyz_c[2] + d * xyz_d[2]) / q};

                            std::array<double, 3> xyz_pq{xyz_p[0] - xyz_q[0],
                                                         xyz_p[1] - xyz_q[1],
                                                         xyz_p[2] - xyz_q[2]};

                            double xx{xyz_pq[0]}, xy{xyz_pq[1]}, xz{xyz_pq[2]};
                            double x = alpha * (xx * xx + xy * xy + xz * xz);

                            boys_f.calcFnx(x, &fnx[0]);

                            double fac = (2.0 * std::pow(M_PI, 2.5) / (p * q * std::sqrt(p + q)));
                            int ofs_row = iab * n_hermite_ab;
                            int ofs_col = icd * n_hermite_cd;

                            calcRInts_ERI4D1<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], n_rints, ofs_row,
                                                       ofs_col, n_rcols, n_rrows, &rints[0]);
                        }

                    for (int munu = 0; munu < n_sph_ab; munu++)
                    {
                        int ofs = munu * n_rrows + iab * n_hermite_ab;
                        for (int tuv = 0; tuv < n_hermite_ab; tuv++)
                            ecoeffs0_bra_ap[ofs + tuv] = (a / p) * ecoeffs0_bra[ofs + tuv];
                    }
                }

            for (int ic = 0, icd = 0; ic < cdepth_c; ic++)
                for (int id = 0; id < cdepth_d; id++, icd++)
                {
                    double c = exps_c[ic];
                    double d = exps_d[id];
                    double q = c + d;

                    for (int kata = 0; kata < n_sph_cd; kata++)
                    {
                        int ofs = kata * n_rcols + icd * n_hermite_cd;
                        for (int tuv = 0; tuv < n_hermite_cd; tuv++)
                            ecoeffs0_ket_cq[ofs + tuv] = (c / q) * ecoeffs0_ket[ofs + tuv];
                    }
                }

            // SHARK integrals

            int n_R_x_E = n_rrows * n_sph_cd;
            int n_E_x_R = n_sph_ab * n_rcols;
            std::vector<double> R_x_E(4 * n_R_x_E, 0);
            std::vector<double> E_x_R(4 * n_E_x_R, 0);

            int m = 4 * n_rrows;
            int n = n_sph_cd;
            int k = n_rcols;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0, &rints[0 * n_rints], k,
                        &ecoeffs0_ket[0], k, 1.0, &R_x_E[0], n);

            m = 4 * n_rcols;
            n = n_sph_ab;
            k = n_rrows;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0, &rints[4 * n_rints], k,
                        &ecoeffs0_bra[0], k, 1.0, &E_x_R[0], n);

            std::vector<double> eri4_batch_raw(12 * n_sph_abcd, 0);

            // bra P
            m = n_sph_ab;
            n = n_sph_cd;
            k = n_rrows;
            std::vector<double> eri4_batch_P(3 * n_sph_abcd, 0);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra[0], k,
                        &R_x_E[0 * n_R_x_E], n, 1.0, &eri4_batch_P[0 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra[0], k,
                        &R_x_E[1 * n_R_x_E], n, 1.0, &eri4_batch_P[1 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra[0], k,
                        &R_x_E[2 * n_R_x_E], n, 1.0, &eri4_batch_P[2 * n_sph_abcd], n);

            // bra A from P and R
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra_ap[0], k,
                        &R_x_E[0 * n_R_x_E], n, 1.0, &eri4_batch_raw[0 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra_ap[0], k,
                        &R_x_E[1 * n_R_x_E], n, 1.0, &eri4_batch_raw[1 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_bra_ap[0], k,
                        &R_x_E[2 * n_R_x_E], n, 1.0, &eri4_batch_raw[2 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3 * m, n, k, 1.0, &ecoeffs1_bra[0], k,
                        &R_x_E[3 * n_R_x_E], n, 1.0, &eri4_batch_raw[0 * n_sph_abcd], n);

            // B from P and A
            cblas_daxpy(3 * n_sph_abcd, 1.0, &eri4_batch_P[0], 1, &eri4_batch_raw[3 * n_sph_abcd], 1);                 // P
            cblas_daxpy(3 * n_sph_abcd, -1.0, &eri4_batch_raw[0 * n_sph_abcd], 1, &eri4_batch_raw[3 * n_sph_abcd], 1); // A

            // ket Q
            m = n_sph_cd;
            n = n_sph_ab;
            k = n_rcols;
            std::vector<double> eri4_batch_Q(3 * n_sph_abcd, 0);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket[0], k,
                        &E_x_R[0 * n_E_x_R], n, 1.0, &eri4_batch_Q[0 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket[0], k,
                        &E_x_R[1 * n_E_x_R], n, 1.0, &eri4_batch_Q[1 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket[0], k,
                        &E_x_R[2 * n_E_x_R], n, 1.0, &eri4_batch_Q[2 * n_sph_abcd], n);

            // ket C from Q and S
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket_cq[0], k,
                        &E_x_R[0 * n_E_x_R], n, 1.0, &eri4_batch_raw[6 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket_cq[0], k,
                        &E_x_R[1 * n_E_x_R], n, 1.0, &eri4_batch_raw[7 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, &ecoeffs0_ket_cq[0], k,
                        &E_x_R[2 * n_E_x_R], n, 1.0, &eri4_batch_raw[8 * n_sph_abcd], n);

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3 * m, n, k, 1.0, &ecoeffs1_ket[0], k,
                        &E_x_R[3 * n_E_x_R], n, 1.0, &eri4_batch_raw[6 * n_sph_abcd], n);

            // D from Q and C
            cblas_daxpy(3 * n_sph_abcd, 1.0, &eri4_batch_Q[0], 1, &eri4_batch_raw[9 * n_sph_abcd], 1);                 // Q
            cblas_daxpy(3 * n_sph_abcd, -1.0, &eri4_batch_raw[6 * n_sph_abcd], 1, &eri4_batch_raw[9 * n_sph_abcd], 1); // C

            std::array<vec4d, 12> eri4_batch;
            for (int ideriv = 0; ideriv < 12; ideriv++)
                eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

            int ofs_norm_a = sp_data_ab.offsets_norms[2 * ipair_ab + 0];
            int ofs_norm_b = sp_data_ab.offsets_norms[2 * ipair_ab + 1];
            int ofs_norm_c = sp_data_cd.offsets_norms[2 * ipair_cd + 0];
            int ofs_norm_d = sp_data_cd.offsets_norms[2 * ipair_cd + 1];
            for (int ideriv = 0; ideriv < 12; ideriv++)
                for (int mu = 0; mu < n_sph_a; mu++)
                    for (int nu = 0; nu < n_sph_b; nu++)
                        for (int ka = 0; ka < n_sph_c; ka++)
                            for (int ta = 0; ta < n_sph_d; ta++)
                            {
                                double norm_a = sp_data_ab.norms[ofs_norm_a + mu];
                                double norm_b = sp_data_ab.norms[ofs_norm_b + nu];
                                double norm_c = sp_data_cd.norms[ofs_norm_c + ka];
                                double norm_d = sp_data_cd.norms[ofs_norm_d + ta];

                                int idx;
                                if (ideriv < 6)
                                    idx = ideriv * n_sph_abcd + mu * (n_sph_b * n_sph_c * n_sph_d) +
                                          nu * (n_sph_c * n_sph_d) + ka * n_sph_d + ta;
                                else
                                    idx = ideriv * n_sph_abcd + ka * (n_sph_d * n_sph_a * n_sph_b) +
                                          ta * (n_sph_a * n_sph_b) + mu * n_sph_b + nu;

                                eri4_batch[ideriv](mu, nu, ka, ta) = norm_a * norm_b * norm_c * norm_d *
                                                                     eri4_batch_raw[idx];
                            }

            return eri4_batch;
        }

        std::array<lible::vec4d, 12> eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                                     const ShellPairData &sp_data_ab,
                                                     const ShellPairData &sp_data_cd,
                                                     const ERI4D1Kernel *eri4d1_kernel);
    }
}