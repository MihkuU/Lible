#pragma once

#include <lible/boys_function.hpp>
#include <lible/cart_exps.hpp>
#include <lible/mcmurchie_davidson.hpp>
#include <lible/shell_pair_data.hpp>
#include <lible/types.hpp>

namespace lible
{
    namespace ints
    {
        namespace two
        {
            vec4d calcERI4(const Structure &structure);

            vec4d calcERI4Shark(const Structure &structure);

            void calcERI4Benchmark(const Structure &structure);

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
                                 const std::vector<MD::IdxsTUV> &idxs_tuv_ab,
                                 const std::vector<MD::IdxsTUV> &idxs_tuv_cd,
                                 const ShellPairData &shell_pair_data_ab,
                                 const ShellPairData &shell_pair_data_cd,
                                 const BoysF &boys_f, arma::dmat &eri4_shells_sph);
        }
    }
}