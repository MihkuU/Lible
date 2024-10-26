#pragma once

#include <lible/types.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <functional>

#include <armadillo>

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

            using kernel_eri4_t = std::function<void(const int cdepth_a, const int cdepth_b,
                                                     const int cdepth_c, const int cdepth_d,
                                                     const double *exps_a, const double *exps_b,
                                                     const double *exps_c, const double *exps_d,
                                                     const double *coords_a, const double *coords_b,
                                                     const double *coords_c, const double *coords_d,
                                                     const double *ecoeffs_ab,
                                                     const double *ecoeffs_cd_tsp,
                                                     double *eri4_batch)>;

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

            kernel_eri4_t deployERI4Kernel(const int la, const int lb, const int lc, const int ld);
        }
    }
}