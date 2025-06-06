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
        namespace two // TODO: delete this file.
        {
            std::vector<double> calcERI2Diagonal(const Structure &structure);

            vec2d calcERI2(const Structure &structure);

            vec2d calcERI4Diagonal(const Structure &structure);

            vec3d calcERI3(const Structure &structure);

            vec4d calcERI4(const Structure &structure);

            void calcERI4Benchmark(const Structure &structure);

            void calcERI4BenchmarkTest(const Structure &structure);


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

            std::array<vec3d, 9> kernelERI3Deriv1Test(const int ipair_ab, const int ishell_c,
                                                      const std::vector<double> &ecoeffs0_bra,
                                                      const std::vector<double> &ecoeffs1_bra,
                                                      const std::vector<double> &ecoeffs0_ket,
                                                      const BoysGrid &boys_grid,
                                                      const ShellPairData &sp_data_ab,
                                                      const ShellData &sh_data_c);

            std::array<vec4d, 12> kernelERI4Deriv1Test(const int ipair_ab, const int ipair_cd,
                                                       const std::vector<double> &ecoeffs0_bra,
                                                       const std::vector<double> &ecoeffs1_bra,
                                                       const std::vector<double> &ecoeffs0_ket,
                                                       const std::vector<double> &ecoeffs1_ket,
                                                       const BoysGrid &boys_grid,
                                                       const ShellPairData &sp_data_ab,
                                                       const ShellPairData &sp_data_cd);

            std::array<vec4d, 12> kernelERI4Deriv1(const int ipair_ab, const int ipair_cd,
                                                   const std::vector<double> &ecoeffs_ab,
                                                   const std::vector<double> &ecoeffs1_ab,
                                                   const std::vector<double> &ecoeffs_cd_tsp,
                                                   const std::vector<double> &ecoeffs1_cd_tsp,
                                                   const BoysGrid &boys_grid,
                                                   const ShellPairData &sp_data_ab,
                                                   const ShellPairData &sp_data_cd);
        }  
    }
}