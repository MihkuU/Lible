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

        template <int la, int lb>
        void calcRInts_ERI2D2(const double alpha, const double fac, const double *fnx,
                              const double *xyz_ab, double *rints);

        template <int lab, int lc>
        void calcRInts_ERI3D1(const double alpha, const double fac, const double *fnx,
                              const double *xyz_pc, double *rints);

        template <int lbra, int lket>
        void calcRInts_ERISOC(const double alpha, const double fac, const double *fnx,
                              const double *xyz_pq, double *rints);

        // ERI kernels for arbitrary L

        vec4d eri4KernelFun(const int ipair_ab, const int ipair_cd,
                            const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                            const ERI4Kernel *eri4_kernel);

        vec3d eri3KernelFun(const int ipair_ab, const int ishell_c,
                            const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                            const ERI3Kernel *eri3_kernel);

        vec2d eri2KernelFun(const int ishell_a, const int ishell_b,
                            const ShellData &sh_data_a, const ShellData &sh_data_b,
                            const ERI2Kernel *eri2_kernel);

        std::array<lible::vec4d, 12> eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                                     const ShellPairData &sp_data_ab,
                                                     const ShellPairData &sp_data_cd,
                                                     const ERI4D1Kernel *eri4d1_kernel);

        std::array<vec3d, 9> eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                             const ShellPairData &sp_data_ab,
                                             const ShellData &sh_data_c,
                                             const ERI3D1Kernel *eri3d1_kernel);

        std::array<vec2d, 6> eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                             const ShellData &sh_data_a,
                                             const ShellData &sh_data_b,
                                             const ERI2D1Kernel *eri2d1_kernel);

        // arr2d<vec4d, 12, 12> eri4d2KernelFun(const int ipair_ab, const int ipair_cd,
        //                                      const ShellPairData &sp_data_ab,
        //                                      const ShellPairData &sp_data_cd,
        //                                      const ERI4D2Kernel *eri4d2_kernel);

        // arr2d<vec3d, 9, 9> eri3d2KernelFun(const int ipair_ab, const int ishell_c,
        //                                    const ShellPairData &sp_data_ab,
        //                                    const ShellData &sh_data_c,
        //                                    const ERI3D2Kernel *eri3d2_kernel);

        arr2d<vec2d, 6, 6> eri2d2KernelFun(const int ishell_a, const int ishell_b,
                                           const ShellData &sh_data_a,
                                           const ShellData &sh_data_b,
                                           const ERI2D2Kernel *eri2d2_kernel);

        std::array<vec4d, 3> eri4socKernelFun(const int ipair_ab, const int ipair_cd,
                                              const ShellPairData &sp_data_ab,
                                              const ShellPairData &sp_data_cd,
                                              const ERI4SOCKernel *eri4soc_kernel);

        std::array<vec3d, 3> eri3socKernelFun(const int ipair_ab, const int ishell_c,
                                              const ShellPairData &sp_data_ab,
                                              const ShellData &sh_data_c,
                                              const ERI3SOCKernel *eri3soc_kernel);

        // Templated ERI kernels for limited L

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

            return eri4_batch;
        }

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

            return eri3_batch;
        }

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

            return eri2_batch;
        }

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
            constexpr int n_rints = n_hermite_ab * n_hermite_cd;
            constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
            constexpr int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;
            constexpr int n_R_x_E = n_hermite_ab * n_sph_cd;

            // Read-in data
            const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
            const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
            const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
            const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
            const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
            const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
            const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
            const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
            const int ofs_E0_ab = sp_data_ab.offsets_ecoeffs[ipair_ab];
            const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
            const int ofs_E0_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];
            const int ofs_E1_cd = sp_data_cd.offsets_ecoeffs_deriv1[ipair_cd];

            const double *exps_a = &sp_data_ab.exps[cofs_a];
            const double *exps_b = &sp_data_ab.exps[cofs_b];
            const double *exps_c = &sp_data_cd.exps[cofs_c];
            const double *exps_d = &sp_data_cd.exps[cofs_d];
            const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
            const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
            const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
            const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
            const double *ecoeffs0_ab = &eri4d1_kernel->ecoeffs0_bra[ofs_E0_ab];
            const double *ecoeffs1_ab = &eri4d1_kernel->ecoeffs1_bra[ofs_E1_ab];
            const double *ecoeffs0_cd = &eri4d1_kernel->ecoeffs0_ket[ofs_E0_cd];
            const double *ecoeffs1_cd = &eri4d1_kernel->ecoeffs1_ket[ofs_E1_cd];

            // SHARK integrals
            std::array<double, labcd + 2> fnx;
            BoysF2<labcd + 1> boys_f;

            std::array<vec4d, 12> eri4_batch;
            for (int ideriv = 0; ideriv < 12; ideriv++)
                eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

            std::array<double, 4 * n_rints> rints;
            for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (int ib = 0; ib < cdepth_b; ib++, iab++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double p = a + b;

                    std::array<double, 7 * n_R_x_E> R_x_E{};
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
                            calcRInts_ERI3D1<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                            std::array<double, 3 * n_R_x_E> I{};

                            int ofs_ecoeffs0_cd = icd * n_ecoeffs_cd;
                            int ofs_ecoeffs1_cd = 3 * icd * n_ecoeffs_cd;
                            shark_mm_ket2<lab, lc, ld>(&rints[0 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &I[0 * n_R_x_E]);
                            shark_mm_ket2<lab, lc, ld>(&rints[1 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &I[1 * n_R_x_E]);
                            shark_mm_ket2<lab, lc, ld>(&rints[2 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &I[2 * n_R_x_E]);

                            shark_mm_ket2<lab, lc, ld>(&rints[3 * n_rints], &ecoeffs0_cd[ofs_ecoeffs0_cd], &R_x_E[3 * n_R_x_E]);

                            shark_mm_ket2<lab, lc, ld>(&rints[3 * n_rints], &ecoeffs1_cd[ofs_ecoeffs1_cd + 0 * n_ecoeffs_cd], &R_x_E[4 * n_R_x_E]);
                            shark_mm_ket2<lab, lc, ld>(&rints[3 * n_rints], &ecoeffs1_cd[ofs_ecoeffs1_cd + 1 * n_ecoeffs_cd], &R_x_E[5 * n_R_x_E]);
                            shark_mm_ket2<lab, lc, ld>(&rints[3 * n_rints], &ecoeffs1_cd[ofs_ecoeffs1_cd + 2 * n_ecoeffs_cd], &R_x_E[6 * n_R_x_E]);

                            cblas_daxpy(3 * n_R_x_E, 1.0, &I[0], 1, &R_x_E[0 * n_R_x_E], 1);
                            cblas_daxpy(3 * n_R_x_E, -(c / q), &I[0], 1, &R_x_E[4 * n_R_x_E], 1);
                        }

                    int ofs_ecoeffs0_ab = iab * n_ecoeffs_ab;
                    int ofs_ecoeffs1_ab = 3 * iab * n_ecoeffs_ab;

                    // P & R
                    std::array<double, 3 * n_sph_abcd> P{};
                    std::array<double, 3 * n_sph_abcd> R{};

                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[0 * n_R_x_E], &P[0 * n_sph_abcd]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[1 * n_R_x_E], &P[1 * n_sph_abcd]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[2 * n_R_x_E], &P[2 * n_sph_abcd]);

                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_ecoeffs1_ab + 0 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[0 * n_sph_abcd]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_ecoeffs1_ab + 1 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[1 * n_sph_abcd]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_ecoeffs1_ab + 2 * n_ecoeffs_ab], &R_x_E[3 * n_R_x_E], &R[2 * n_sph_abcd]);

                    // A
                    cblas_daxpy(n_sph_abcd, (a / p), &P[0 * n_sph_abcd], 1, &eri4_batch[0][0], 1);
                    cblas_daxpy(n_sph_abcd, (a / p), &P[1 * n_sph_abcd], 1, &eri4_batch[1][0], 1);
                    cblas_daxpy(n_sph_abcd, (a / p), &P[2 * n_sph_abcd], 1, &eri4_batch[2][0], 1);

                    cblas_daxpy(n_sph_abcd, 1.0, &R[0 * n_sph_abcd], 1, &eri4_batch[0][0], 1);
                    cblas_daxpy(n_sph_abcd, 1.0, &R[1 * n_sph_abcd], 1, &eri4_batch[1][0], 1);
                    cblas_daxpy(n_sph_abcd, 1.0, &R[2 * n_sph_abcd], 1, &eri4_batch[2][0], 1);

                    // B
                    cblas_daxpy(n_sph_abcd, (b / p), &P[0 * n_sph_abcd], 1, &eri4_batch[3][0], 1);
                    cblas_daxpy(n_sph_abcd, (b / p), &P[1 * n_sph_abcd], 1, &eri4_batch[4][0], 1);
                    cblas_daxpy(n_sph_abcd, (b / p), &P[2 * n_sph_abcd], 1, &eri4_batch[5][0], 1);

                    cblas_daxpy(n_sph_abcd, -1.0, &R[0 * n_sph_abcd], 1, &eri4_batch[3][0], 1);
                    cblas_daxpy(n_sph_abcd, -1.0, &R[1 * n_sph_abcd], 1, &eri4_batch[4][0], 1);
                    cblas_daxpy(n_sph_abcd, -1.0, &R[2 * n_sph_abcd], 1, &eri4_batch[5][0], 1);

                    // C
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[4 * n_R_x_E], &eri4_batch[6][0]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[5 * n_R_x_E], &eri4_batch[7][0]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs0_ab[ofs_ecoeffs0_ab], &R_x_E[6 * n_R_x_E], &eri4_batch[8][0]);
                }

            // D
            for (int ideriv = 9; ideriv < 12; ideriv++)
                eri4_batch[ideriv] = -1 * (eri4_batch[ideriv - 9] + eri4_batch[ideriv - 6] + eri4_batch[ideriv - 3]);

            return eri4_batch;
        }

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

                    std::array<double, 4 * n_R_x_E> R_x_E{};
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

                        int ofs_ecoeffs_c = ic * n_ecoeffs_c;
                        shark_mm_ket1<lab, lc>(&rints[0 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &R_x_E[0 * n_R_x_E]);
                        shark_mm_ket1<lab, lc>(&rints[1 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &R_x_E[1 * n_R_x_E]);
                        shark_mm_ket1<lab, lc>(&rints[2 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &R_x_E[2 * n_R_x_E]);                        

                        shark_mm_ket1<lab, lc>(&rints[3 * n_rints], &ecoeffs0_c[ofs_ecoeffs_c], &R_x_E[3 * n_R_x_E]);
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

                    // C
                    cblas_daxpy(n_sph_abc, -1.0, &P[0 * n_sph_abc], 1, &eri3_batch[6][0], 1);
                    cblas_daxpy(n_sph_abc, -1.0, &P[1 * n_sph_abc], 1, &eri3_batch[7][0], 1);
                    cblas_daxpy(n_sph_abc, -1.0, &P[2 * n_sph_abc], 1, &eri3_batch[8][0], 1);
                }

            // B
            for (int ideriv = 3; ideriv < 6; ideriv++)
                eri3_batch[ideriv] = -1 * (eri3_batch[ideriv - 3] + eri3_batch[ideriv + 3]);

            return eri3_batch;
        }

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
                std::array<double, 3 * n_R_x_E> R_x_E{};                
                for (int ib = 0; ib < cdepth_b; ib++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];

                    double alpha = a * b / (a + b);
                    double x = alpha * xyz_ab_dot;
                    boys_f.calcFnx(x, &fnx[0]);

                    double fac = (2.0 * std::pow(M_PI, 2.5) / (a * b * std::sqrt(a + b)));
                    calcRInts_ERI2D1<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);             

                    int ofs_ecoeffs_b = ib * n_ecoeffs_b;
                    shark_mm_ket1<la, lb>(&rints[0 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &R_x_E[0 * n_R_x_E]);
                    shark_mm_ket1<la, lb>(&rints[1 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &R_x_E[1 * n_R_x_E]);
                    shark_mm_ket1<la, lb>(&rints[2 * n_rints], &ecoeffs_b[ofs_ecoeffs_b], &R_x_E[2 * n_R_x_E]);

                }

                // A
                int ofs_ecoeffs_a = ia * n_ecoeffs_a;
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[0 * n_R_x_E], &eri2_batch[0][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[1 * n_R_x_E], &eri2_batch[1][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_ecoeffs_a], &R_x_E[2 * n_R_x_E], &eri2_batch[2][0]);
            }

            // B
            for (int ideriv = 3; ideriv < 6; ideriv++)
                eri2_batch[ideriv] = -1 * eri2_batch[ideriv - 3];

            return eri2_batch;
        }

        // template <int la, int lb, int lc>
        // arr2d<vec3d, 9, 9> eri3d2KernelFun(const int ipair_ab, const int ishell_c,
        //                                       const ShellPairData &sp_data_ab,
        //                                       const ShellData &sh_data_c,
        //                                       const ERI3D2Kernel *eri3d2_kernel)
        // {
        //     arr2d<vec3d, 9, 9> eri3_batch;

        //     return eri3_batch;
        // }

        template <int la, int lb>
        arr2d<vec2d, 6, 6> eri2d2KernelFun(const int ishell_a, const int ishell_b,
                                           const ShellData &sh_data_a,
                                           const ShellData &sh_data_b,
                                           const ERI2D2Kernel *eri2d2_kernel)
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
            const double *ecoeffs_a = &eri2d2_kernel->ecoeffs_bra[ofs_E_a];
            const double *ecoeffs_b = &eri2d2_kernel->ecoeffs_ket[ofs_E_b];

            // SHARK integrals
            std::array<double, lab + 3> fnx;
            BoysF2<lab + 2> boys_f;

            std::array<double, 3> xyz_ab{coords_a[0] - coords_b[0],
                                         coords_a[1] - coords_b[1],
                                         coords_a[2] - coords_b[2]};
            double dx{xyz_ab[0]}, dy{xyz_ab[1]}, dz{xyz_ab[2]};
            double xyz_ab_dot = dx * dx + dy * dy + dz * dz;

            arr2d<vec2d, 6, 6> eri2_batch;
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    eri2_batch[i][j] = vec2d(Fill(0), n_sph_a, n_sph_b);

            std::array<double, 6 * n_rints> rints{};
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
                    calcRInts_ERI2D2<la, lb>(alpha, fac, &fnx[0], &xyz_ab[0], &rints[0]);

                    int ofs_e0_b = ib * n_ecoeffs_b;
                    shark_mm_ket1<la, lb>(&rints[0 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[0 * n_R_x_E]); // 200
                    shark_mm_ket1<la, lb>(&rints[1 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[1 * n_R_x_E]); // 110
                    shark_mm_ket1<la, lb>(&rints[2 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[2 * n_R_x_E]); // 101
                    shark_mm_ket1<la, lb>(&rints[3 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[3 * n_R_x_E]); // 020
                    shark_mm_ket1<la, lb>(&rints[4 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[4 * n_R_x_E]); // 011
                    shark_mm_ket1<la, lb>(&rints[5 * n_rints], &ecoeffs_b[ofs_e0_b], &R_x_E[5 * n_R_x_E]); // 002
                }

                int ofs_e0_a = ia * n_ecoeffs_a;

                // AA upper triangle
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_e0_a], &R_x_E[0 * n_R_x_E], &eri2_batch[0][0][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_e0_a], &R_x_E[1 * n_R_x_E], &eri2_batch[0][1][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_e0_a], &R_x_E[2 * n_R_x_E], &eri2_batch[0][2][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_e0_a], &R_x_E[3 * n_R_x_E], &eri2_batch[1][1][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_e0_a], &R_x_E[4 * n_R_x_E], &eri2_batch[1][2][0]);
                shark_mm_bra1<la, lb>(&ecoeffs_a[ofs_e0_a], &R_x_E[5 * n_R_x_E], &eri2_batch[2][2][0]);
            }

            // Complete AA
            for (int ideriv = 0; ideriv < 3; ideriv++)
                for (int jderiv = ideriv + 1; jderiv < 3; jderiv++)
                    eri2_batch[jderiv][ideriv] = eri2_batch[ideriv][jderiv];

            // AB
            for (int ideriv = 0; ideriv < 3; ideriv++)
                for (int jderiv = 0; jderiv < 3; jderiv++)
                {
                    eri2_batch[jderiv][3 + ideriv] = -1 * eri2_batch[ideriv][jderiv];
                    eri2_batch[3 + ideriv][jderiv] = -1 * eri2_batch[ideriv][jderiv];
                }

            // BB
            for (int ideriv = 3; ideriv < 6; ideriv++)
                for (int jderiv = 3; jderiv < 6; jderiv++)
                    eri2_batch[ideriv][jderiv] = -1 * eri2_batch[ideriv - 3][jderiv];

            return eri2_batch;
        }

        template <int la, int lb, int lc, int ld>
        std::array<vec4d, 3> eri4socKernelFun(const int ipair_ab, const int ipair_cd,
                                              const ShellPairData &sp_data_ab,
                                              const ShellPairData &sp_data_cd,
                                              const ERI4SOCKernel *eri4soc_kernel)
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
            constexpr int n_rints = n_hermite_ab * n_hermite_cd;
            constexpr int n_ecoeffs_ab = n_sph_ab * n_hermite_ab;
            constexpr int n_ecoeffs_cd = n_sph_cd * n_hermite_cd;
            constexpr int n_R_x_E = n_hermite_ab * n_sph_cd;

            // Read-in data
            const int cdepth_a = sp_data_ab.cdepths[2 * ipair_ab + 0];
            const int cdepth_b = sp_data_ab.cdepths[2 * ipair_ab + 1];
            const int cdepth_c = sp_data_cd.cdepths[2 * ipair_cd + 0];
            const int cdepth_d = sp_data_cd.cdepths[2 * ipair_cd + 1];
            const int cofs_a = sp_data_ab.coffsets[2 * ipair_ab + 0];
            const int cofs_b = sp_data_ab.coffsets[2 * ipair_ab + 1];
            const int cofs_c = sp_data_cd.coffsets[2 * ipair_cd + 0];
            const int cofs_d = sp_data_cd.coffsets[2 * ipair_cd + 1];
            const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
            const int ofs_E0_cd = sp_data_cd.offsets_ecoeffs[ipair_cd];

            const double *exps_a = &sp_data_ab.exps[cofs_a];
            const double *exps_b = &sp_data_ab.exps[cofs_b];
            const double *exps_c = &sp_data_cd.exps[cofs_c];
            const double *exps_d = &sp_data_cd.exps[cofs_d];
            const double *xyz_a = &sp_data_ab.coords[6 * ipair_ab + 0];
            const double *xyz_b = &sp_data_ab.coords[6 * ipair_ab + 3];
            const double *xyz_c = &sp_data_cd.coords[6 * ipair_cd + 0];
            const double *xyz_d = &sp_data_cd.coords[6 * ipair_cd + 3];
            const double *ecoeffs1_ab = &eri4soc_kernel->ecoeffs1_bra[ofs_E1_ab];
            const double *ecoeffs0_cd = &eri4soc_kernel->ecoeffs0_ket[ofs_E0_cd];

            // SHARK integrals
            std::array<double, labcd + 2> fnx;
            BoysF2<labcd + 1> boys_f;

            std::array<vec4d, 3> eri4_batch;
            for (int ideriv = 0; ideriv < 3; ideriv++)
                eri4_batch[ideriv] = vec4d(Fill(0), n_sph_a, n_sph_b, n_sph_c, n_sph_d);

            std::array<double, 3 * n_rints> rints;
            for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (int ib = 0; ib < cdepth_b; ib++, iab++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double p = a + b;

                    std::array<double, 3 * n_R_x_E> R_x_E{};
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
                            calcRInts_ERISOC<lab, lcd>(alpha, fac, &fnx[0], &xyz_pq[0], &rints[0]);

                            int ofs_e0_cd = icd * n_ecoeffs_cd;
                            shark_mm_ket2<lab, lc, ld>(&rints[0 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &R_x_E[0 * n_R_x_E]);
                            shark_mm_ket2<lab, lc, ld>(&rints[1 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &R_x_E[1 * n_R_x_E]);
                            shark_mm_ket2<lab, lc, ld>(&rints[2 * n_rints], &ecoeffs0_cd[ofs_e0_cd], &R_x_E[2 * n_R_x_E]);
                        }

                    int ofs_e1_ab = 3 * iab * n_ecoeffs_ab;
                    int ofs_e1_0 = ofs_e1_ab + 0 * n_ecoeffs_ab;
                    int ofs_e1_1 = ofs_e1_ab + 1 * n_ecoeffs_ab;
                    int ofs_e1_2 = ofs_e1_ab + 2 * n_ecoeffs_ab;

                    // P & R
                    std::array<double, n_sph_abcd> PR_yz{};
                    std::array<double, n_sph_abcd> PR_zy{};
                    std::array<double, n_sph_abcd> PR_zx{};
                    std::array<double, n_sph_abcd> PR_xz{};
                    std::array<double, n_sph_abcd> PR_xy{};
                    std::array<double, n_sph_abcd> PR_yx{};

                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_e1_2], &R_x_E[1 * n_R_x_E], &PR_yz[0]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_e1_1], &R_x_E[2 * n_R_x_E], &PR_zy[0]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_e1_0], &R_x_E[2 * n_R_x_E], &PR_zx[0]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_e1_2], &R_x_E[0 * n_R_x_E], &PR_xz[0]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_e1_1], &R_x_E[0 * n_R_x_E], &PR_xy[0]);
                    shark_mm_bra2<la, lb, lc, ld>(&ecoeffs1_ab[ofs_e1_0], &R_x_E[1 * n_R_x_E], &PR_yx[0]);

                    // A & B
                    cblas_daxpy(n_sph_abcd, 1.0, &PR_zy[0], 1, &eri4_batch[0][0], 1);
                    cblas_daxpy(n_sph_abcd, 1.0, &PR_xz[0], 1, &eri4_batch[1][0], 1);
                    cblas_daxpy(n_sph_abcd, 1.0, &PR_yx[0], 1, &eri4_batch[2][0], 1);

                    cblas_daxpy(n_sph_abcd, -1.0, &PR_yz[0], 1, &eri4_batch[0][0], 1);
                    cblas_daxpy(n_sph_abcd, -1.0, &PR_zx[0], 1, &eri4_batch[1][0], 1);
                    cblas_daxpy(n_sph_abcd, -1.0, &PR_xy[0], 1, &eri4_batch[2][0], 1);
                }

            return eri4_batch;
        }

        template <int la, int lb, int lc>
        std::array<vec3d, 3> eri3socKernelFun(const int ipair_ab, const int ishell_c,
                                              const ShellPairData &sp_data_ab,
                                              const ShellData &sh_data_c,
                                              const ERI3SOCKernel *eri3soc_kernel)
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
            const int ofs_E1_ab = sp_data_ab.offsets_ecoeffs_deriv1[ipair_ab];
            const int ofs_E0_c = sh_data_c.offsets_ecoeffs[ishell_c];

            const double *exps_a = &sp_data_ab.exps[cofs_a];
            const double *exps_b = &sp_data_ab.exps[cofs_b];
            const double *exps_c = &sh_data_c.exps[cofs_c];
            const double *coords_a = &sp_data_ab.coords[6 * ipair_ab];
            const double *coords_b = &sp_data_ab.coords[6 * ipair_ab + 3];
            const double *coords_c = &sh_data_c.coords[3 * ishell_c];
            const double *ecoeffs1_ab = &eri3soc_kernel->ecoeffs1_bra[ofs_E1_ab];
            const double *ecoeffs0_c = &eri3soc_kernel->ecoeffs0_ket[ofs_E0_c];

            // SHARK integrals
            std::array<double, labc + 2> fnx;
            BoysF2<labc + 1> boys_f;

            std::array<vec3d, 3> eri3_batch;
            for (int ideriv = 0; ideriv < 3; ideriv++)
                eri3_batch[ideriv] = vec3d(Fill(0), n_sph_a, n_sph_b, n_sph_c);

            std::array<double, 3 * n_rints> rints;
            for (int ia = 0, iab = 0; ia < cdepth_a; ia++)
                for (int ib = 0; ib < cdepth_b; ib++, iab++)
                {
                    double a = exps_a[ia];
                    double b = exps_b[ib];
                    double p = a + b;

                    std::array<double, 3 * n_R_x_E> R_x_E{};
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
                        calcRInts_ERISOC<lab, lc>(alpha, fac, &fnx[0], &xyz_pc[0], &rints[0]);

                        int ofs_e0_c = ic * n_ecoeffs_c;
                        shark_mm_ket1<lab, lc>(&rints[0 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[0 * n_R_x_E]);
                        shark_mm_ket1<lab, lc>(&rints[1 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[1 * n_R_x_E]);
                        shark_mm_ket1<lab, lc>(&rints[2 * n_rints], &ecoeffs0_c[ofs_e0_c], &R_x_E[2 * n_R_x_E]);
                    }

                    int ofs_e1_ab = 3 * iab * n_ecoeffs_ab;
                    int ofs_e1_0 = ofs_e1_ab + 0 * n_ecoeffs_ab;
                    int ofs_e1_1 = ofs_e1_ab + 1 * n_ecoeffs_ab;
                    int ofs_e1_2 = ofs_e1_ab + 2 * n_ecoeffs_ab;

                    // P & R
                    std::array<double, n_sph_abc> PR_yz{};
                    std::array<double, n_sph_abc> PR_zy{};
                    std::array<double, n_sph_abc> PR_zx{};
                    std::array<double, n_sph_abc> PR_xz{};
                    std::array<double, n_sph_abc> PR_xy{};
                    std::array<double, n_sph_abc> PR_yx{};

                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_e1_2], &R_x_E[1 * n_R_x_E], &PR_yz[0]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_e1_1], &R_x_E[2 * n_R_x_E], &PR_zy[0]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_e1_0], &R_x_E[2 * n_R_x_E], &PR_zx[0]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_e1_2], &R_x_E[0 * n_R_x_E], &PR_xz[0]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_e1_1], &R_x_E[0 * n_R_x_E], &PR_xy[0]);
                    shark_mm_bra2<la, lb, lc>(&ecoeffs1_ab[ofs_e1_0], &R_x_E[1 * n_R_x_E], &PR_yx[0]);

                    // A & B
                    cblas_daxpy(n_sph_abc, 1.0, &PR_zy[0], 1, &eri3_batch[0][0], 1);
                    cblas_daxpy(n_sph_abc, 1.0, &PR_xz[0], 1, &eri3_batch[1][0], 1);
                    cblas_daxpy(n_sph_abc, 1.0, &PR_yx[0], 1, &eri3_batch[2][0], 1);

                    cblas_daxpy(n_sph_abc, -1.0, &PR_yz[0], 1, &eri3_batch[0][0], 1);
                    cblas_daxpy(n_sph_abc, -1.0, &PR_zx[0], 1, &eri3_batch[1][0], 1);
                    cblas_daxpy(n_sph_abc, -1.0, &PR_xy[0], 1, &eri3_batch[2][0], 1);
                }

            return eri3_batch;
        }
    }
}