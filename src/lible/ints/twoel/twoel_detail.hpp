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

            void calcERI4BenchmarkTest(const Structure &structure);

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

            using kernel_eri2d1_t = std::function<std::array<vec2d, 6>(
                const int ishell_a, const int ishell_b,
                const std::vector<double> &ecoeffs_a,
                const std::vector<double> &ecoeffs_b_tsp,
                const ShellData &sh_data_a, const ShellData &sh_data_b)>;

            kernel_eri2d1_t deployERI2Deriv1Kernel(const int la, const int lb);

            template <int la, int lb>
            std::array<vec2d, 6> eri2d1Kernel(const int ishell_a, const int ishell_b,
                                              const std::vector<double> &ecoeffs_a,
                                              const std::vector<double> &ecoeffs_b_tsp,
                                              const ShellData &sh_data_a, const ShellData &sh_data_b);

        //     template<int la, int lb>            
        //     struct ERI4Kernel
        //     {
        //         // std::vector<>
        //         ERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd)
        //         {
        //             // max cdepth ab _x_ max cdepth cd
        //         }

        //         std::vector<double> ecoeffs_bra;
        //         std::vector<double> ecoeffs_ket;
        //         std::vector<double> rints;
        //         std::vector<double> R_x_E;

        //         void operator()() const
        //         {
        //         }
        //     };
        }

        ///////////////////////////////////////////// NEW SHIT

        // struct ERI4Kernel;

        // template <int la, int lb, int lc, int ld>
        // vec4d eri4KernelFun(const int ipair_ab, const int ipair_cd,
        //                     const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
        //                     const ERI4Kernel *eri4_kernel);        
    }
}