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

            void calcERI4Benchmark(const Structure &structure);
            
            using eri4_kernel_t = std::function<void(void)>;

            eri4_kernel_t deployERI4Kernel(const int lab, const int lcd);            
        }
    }
}