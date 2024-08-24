#pragma once

#include <lible/types.hpp>
#include <lible/ints/boys_function.hpp>
#include <lible/ints/cart_exps.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <armadillo>

namespace lible
{
    namespace ints
    {
        namespace two
        {
            vec2d calcERI2(const Structure &structure);

            vec3d calcERI3(const Structure &structure);

            vec4d calcERI4(const Structure &structure);

            // TODO: fix the names of the following functions once everything has matured enough

            vec4d calcERI4Shark(const Structure &structure);

            vec4d calcERI4SharkFlat(const Structure &structure);            

            void calcERI4Benchmark(const Structure &structure);

            void calcERI4BenchmarkShark(const Structure &structure);

            void calcERI4BenchmarkSharkFlat(const Structure &structure);

            // TODO: perhaps declare the functions in the x.cpp file.
            void kernelERI4(const int lab, const int lcd,
                            const size_t ipair_ab, const size_t ipair_cd,
                            const std::vector<std::vector<vec4d>> &ecoeffs_lalb,
                            const std::vector<std::vector<vec4d>> &ecoeffs_lcld,
                            const std::vector<CartExps> &cart_exps_a,
                            const std::vector<CartExps> &cart_exps_b,
                            const std::vector<CartExps> &cart_exps_c,
                            const std::vector<CartExps> &cart_exps_d,
                            const ShellPairData &shell_pair_data_ab,
                            const ShellPairData &shell_pair_data_cd,
                            const BoysF &boys_f, vec4d &eri4_shells_cart);

            void kernelERI4Shark(const int lab, const int lcd,
                                 const size_t ipair_ab, const size_t ipair_cd,
                                 const std::vector<std::vector<arma::dmat>> &ecoeffs_lalb,
                                 const std::vector<std::vector<arma::dmat>> &ecoeffs_lcld,
                                 const std::vector<IdxsTUV> &idxs_tuv_ab,
                                 const std::vector<IdxsTUV> &idxs_tuv_cd,
                                 const ShellPairData &shell_pair_data_ab,
                                 const ShellPairData &shell_pair_data_cd,
                                 const BoysF &boys_f, arma::dmat &eri4_shells_sph);

            // void kernelERI4SharkFlat(const int lab, const int lcd,
            //                          const size_t ipair_ab, const size_t ipair_cd,
            //                          const std::vector<double> &ecoeffs_lalb,
            //                          const std::vector<double> &ecoeffs_lcld,
            //                          const std::vector<MD::IdxsTUV> &idxs_tuv_ab,
            //                          const std::vector<MD::IdxsTUV> &idxs_tuv_cd,
            //                          const ShellPairData &shell_pair_data_ab,
            //                          const ShellPairData &shell_pair_data_cd,
            //                          const BoysF &boys_f, std::vector<double> &eri4_shells_sph,
            //                          std::vector<double> &rints, std::vector<double> &fnx,
            //                          vec4d &rints_tmp);
        }
    }
}