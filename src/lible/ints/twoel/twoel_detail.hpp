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
        namespace two
        {
            std::vector<double> calcERI2Diagonal(const Structure &structure);

            vec2d calcERI2(const Structure &structure);

            vec2d calcERI4Diagonal(const Structure &structure);

            vec3d calcERI3(const Structure &structure);

            vec4d calcERI4(const Structure &structure);

            vec4d calcERI4New(const Structure &structure);

            void calcERI4Benchmark(const Structure &structure);

            void calcERI4BenchmarkNew(const Structure &structure);

            template <int la, int lb>
            void eri2Kernel(const int cdepth_a, const int cdepth_b,
                            const double *exps_a, const double *exps_b,
                            const double *coords_a, const double *coords_b,
                            const double *ecoeffs_a, const double *ecoeffs_b_tsp,
                            double *eri2_batch);

            template <int la, int lb, int lc>
            void eri3Kernel(const int cdepth_a, const int cdepth_b,
                            const int cdepth_c, const double *exps_a,
                            const double *exps_b, const double *exps_c,
                            const double *coords_a, const double *coords_b,
                            const double *coords_c, const double *ecoeffs_ab,
                            const double *ecoeffs_c, double *eri3_batch);

            template <int la, int lb, int lc, int ld>
            void eri4Kernel(const int cdepth_a, const int cdepth_b,
                            const int cdepth_c, const int cdepth_d,
                            const double *exps_a, const double *exps_b,
                            const double *exps_c, const double *exps_d,
                            const double *coords_a, const double *coords_b,
                            const double *coords_c, const double *coords_d,
                            const double *ecoeffs_ab,
                            const double *ecoeffs_cd_tsp,
                            double *eri4_batch);

            using kernel_eri2_t = std::function<void(const int cdepth_a, const int cdepth_b,
                                                     const double *exps_a, const double *exps_b,
                                                     const double *coords_a, const double *coords_b,
                                                     const double *ecoeffs_a,
                                                     const double *ecoeffs_b_tsp,
                                                     double *eri2_batch)>;

            using kernel_eri3_t = std::function<void(const int cdepth_a, const int cdepth_b,
                                                     const int cdepth_c, const double *exps_a,
                                                     const double *exps_b, const double *exps_c,
                                                     const double *coords_a, const double *coords_b,
                                                     const double *coords_c, const double *ecoeffs_ab,
                                                     const double *ecoeffs_c, double *eri3_batch)>;

            using kernel_eri4_t = std::function<void(const int cdepth_a, const int cdepth_b,
                                                     const int cdepth_c, const int cdepth_d,
                                                     const double *exps_a, const double *exps_b,
                                                     const double *exps_c, const double *exps_d,
                                                     const double *coords_a, const double *coords_b,
                                                     const double *coords_c, const double *coords_d,
                                                     const double *ecoeffs_ab,
                                                     const double *ecoeffs_cd_tsp,
                                                     double *eri4_batch)>;

            kernel_eri2_t deployERI2Kernel(const int la, const int lb);

            kernel_eri3_t deployERI3Kernel(const int la, const int lb, const int lc);

            kernel_eri4_t deployERI4Kernel(const int la, const int lb, const int lc, const int ld);

            void kernelERI2Deriv1(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                                  const double *exps_a, const double *exps_b, const double *coords_a,
                                  const double *coords_b, const double *ecoeffs_a,
                                  const double *ecoeffs_b_tsp, const double *norms_a,
                                  const double *norms_b, const BoysGrid &boys_grid,
                                  double *eri2_batch);

            void kernelERI3Deriv1(const int la, const int lb, const int lc,
                                  const int cdepth_a, const int cdepth_b, const int cdepth_c,
                                  const double *exps_a, const double *exps_b, const double *exps_c,
                                  const double *coords_a, const double *coords_b, 
                                  const double *coords_c, const double *ecoeffs_ab, 
                                  const double *ecoeffs_deriv1_ab, const double *ecoeffs_c, 
                                  const double *norms_a, const double *norms_b, 
                                  const double *norms_c, const BoysGrid &boys_grid, 
                                  double *eri3_batch);
        }
    }
}