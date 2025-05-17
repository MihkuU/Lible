#pragma once

#include <lible/types.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <functional>

namespace lible
{
    namespace ints
    {
        namespace two // TODO: remove two-namespace
        {
            std::vector<double> calcERI2Diagonal(const Structure &structure);

            vec2d calcERI2(const Structure &structure);

            vec2d calcERI4Diagonal(const Structure &structure);

            vec3d calcERI3(const Structure &structure);

            vec4d calcERI4(const Structure &structure);

            // vec4d calcERI4New(const Structure &structure);

            void calcERI4Benchmark(const Structure &structure);

            void calcERI4BenchmarkNew(const Structure &structure);

            // TODO: remove E-coeffs arg?
            template <int la, int lb>
            vec2d eri2Kernel(const int ishell_a, const int ishell_b,
                             const std::vector<double> &ecoeffs_a,
                             const std::vector<double> &ecoeffs_b_tsp,
                             const ShellData &sh_data_a, const ShellData &sh_data_b);

            // TODO: remove E-coeffs arg?
            template <int la, int lb, int lc>
            vec3d eri3Kernel(const int ipair_ab, const int ishell_c,
                             const std::vector<double> &ecoeffs_ab,
                             const std::vector<double> &ecoeffs_c,
                             const ShellPairData &sp_data_ab, const ShellData &sh_data_c);

            // TODO: remove E-coeffs arg?
            template <int la, int lb, int lc, int ld>
            vec4d eri4Kernel(const int ipair_ab, const int ipair_cd,
                             const std::vector<double> &ecoeffs_ab,
                             const std::vector<double> &ecoeffs_cd_tsp,
                             const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd);

            using kernel_eri2_t = std::function<vec2d(const int ishell_a, const int ishell_b,
                                                      const std::vector<double> &ecoeffs_a,
                                                      const std::vector<double> &ecoeffs_b_tsp,
                                                      const ShellData &sh_data_a,
                                                      const ShellData &sh_data_b)>;

            using kernel_eri3_t = std::function<vec3d(const int ipair_ab, const int ishell_c,
                                                      const std::vector<double> &ecoeffs_ab,
                                                      const std::vector<double> &ecoeffs_c,
                                                      const ShellPairData &sp_data_ab,
                                                      const ShellData &sh_data_c)>;

            using kernel_eri4_t = std::function<vec4d(const int ipair_ab, const int ipair_cd,
                                                      const std::vector<double> &ecoeffs_ab,
                                                      const std::vector<double> &ecoeffs_cd_tsp,
                                                      const ShellPairData &sp_data_ab,
                                                      const ShellPairData &sp_data_cd)>;

            kernel_eri2_t deployERI2Kernel(const int la, const int lb);

            kernel_eri3_t deployERI3Kernel(const int la, const int lb, const int lc);

            kernel_eri4_t deployERI4Kernel(const int la, const int lb, const int lc, const int ld);

            std::array<vec2d, 6> kernelERI2Deriv1(const int ishell_a, const int ishell_b,
                                                  const std::vector<double> &ecoeffs_a,
                                                  const std::vector<double> &ecoeffs_b_tsp,
                                                  const BoysGrid &boys_grid,
                                                  const ShellData &sh_data_a,
                                                  const ShellData &sh_data_b);

            std::array<vec3d, 9> kernelERI3Deriv1(const int ipair_ab, const int ishell_c,
                                                  const std::vector<double> &ecoeffs_ab,
                                                  const std::vector<double> &ecoeffs1_ab,
                                                  const std::vector<double> &ecoeffs_c,
                                                  const BoysGrid &boys_grid,
                                                  const ShellPairData &sp_data_ab,
                                                  const ShellData &sh_data_c);

            std::array<vec4d, 12> kernelERI4Deriv1(const int ipair_ab, const int ipair_cd,
                                                   const std::vector<double> &ecoeffs_ab,
                                                   const std::vector<double> &ecoeffs1_ab,
                                                   const std::vector<double> &ecoeffs_cd_tsp,
                                                   const std::vector<double> &ecoeffs1_cd_tsp,
                                                   const BoysGrid &boys_grid,
                                                   const ShellPairData &sp_data_ab,
                                                   const ShellPairData &sp_data_cd);

            // template <int la, int lb>
            // void eri2Kernel(const int cdepth_a, const int cdepth_b,
            //                 const double *exps_a, const double *exps_b,
            //                 const double *coords_a, const double *coords_b,
            //                 const double *ecoeffs_a, const double *ecoeffs_b_tsp,
            //                 double *eri2_batch);

            // template <int la, int lb, int lc>
            // void eri3Kernel(const int cdepth_a, const int cdepth_b,
            //                 const int cdepth_c, const double *exps_a,
            //                 const double *exps_b, const double *exps_c,
            //                 const double *coords_a, const double *coords_b,
            //                 const double *coords_c, const double *ecoeffs_ab,
            //                 const double *ecoeffs_c, double *eri3_batch);

            // template <int la, int lb, int lc, int ld>
            // void eri4Kernel(const int cdepth_a, const int cdepth_b,
            //                 const int cdepth_c, const int cdepth_d,
            //                 const double *exps_a, const double *exps_b,
            //                 const double *exps_c, const double *exps_d,
            //                 const double *coords_a, const double *coords_b,
            //                 const double *coords_c, const double *coords_d,
            //                 const double *ecoeffs_ab,
            //                 const double *ecoeffs_cd_tsp,
            //                 double *eri4_batch);

            // using kernel_eri2_t = std::function<void(const int cdepth_a, const int cdepth_b,
            //                                          const double *exps_a, const double *exps_b,
            //                                          const double *coords_a, const double *coords_b,
            //                                          const double *ecoeffs_a,
            //                                          const double *ecoeffs_b_tsp,
            //                                          double *eri2_batch)>;

            // using kernel_eri3_t = std::function<void(const int cdepth_a, const int cdepth_b,
            //                                          const int cdepth_c, const double *exps_a,
            //                                          const double *exps_b, const double *exps_c,
            //                                          const double *coords_a, const double *coords_b,
            //                                          const double *coords_c, const double *ecoeffs_ab,
            //                                          const double *ecoeffs_c, double *eri3_batch)>;

            // using kernel_eri4_t = std::function<void(const int cdepth_a, const int cdepth_b,
            //                                          const int cdepth_c, const int cdepth_d,
            //                                          const double *exps_a, const double *exps_b,
            //                                          const double *exps_c, const double *exps_d,
            //                                          const double *coords_a, const double *coords_b,
            //                                          const double *coords_c, const double *coords_d,
            //                                          const double *ecoeffs_ab,
            //                                          const double *ecoeffs_cd_tsp,
            //                                          double *eri4_batch)>;

            // void kernelERI2Deriv1(const int la, const int lb, const int cdepth_a, const int cdepth_b,
            //                       const double *exps_a, const double *exps_b, const double *coords_a,
            //                       const double *coords_b, const double *ecoeffs_a,
            //                       const double *ecoeffs_b_tsp, const double *norms_a,
            //                       const double *norms_b, const BoysGrid &boys_grid,
            //                       double *eri2_batch);

            // void kernelERI3Deriv1(const int la, const int lb, const int lc,
            //                       const int cdepth_a, const int cdepth_b, const int cdepth_c,
            //                       const double *exps_a, const double *exps_b, const double *exps_c,
            //                       const double *coords_a, const double *coords_b,
            //                       const double *coords_c, const double *ecoeffs_ab,
            //                       const double *ecoeffs_deriv1_ab, const double *ecoeffs_c,
            //                       const double *norms_a, const double *norms_b,
            //                       const double *norms_c, const BoysGrid &boys_grid,
            //                       double *eri3_batch);

            // void kernelERI4Deriv1(const int la, const int lb, const int lc, const int ld,
            //                       const int cdepth_a, const int cdepth_b, const int cdepth_c,
            //                       const int cdepth_d, const double *exps_a, const double *exps_b,
            //                       const double *exps_c, const double *exps_d,
            //                       const double *xyz_a, const double *xyz_b,
            //                       const double *xyz_c, const double *xyz_d,
            //                       const double *ecoeffs_ab, const double *ecoeffs_deriv1_ab,
            //                       const double *ecoeffs_cd_tsp, const double *ecoeffs_deriv1_cd_tsp,
            //                       const double *norms_a, const double *norms_b,
            //                       const double *norms_c, const double *norms_d,
            //                       const BoysGrid &boys_grid, double *eri4_batch);
        }
    }
}