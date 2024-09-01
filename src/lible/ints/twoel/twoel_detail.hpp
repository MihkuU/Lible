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

            vec4d calcERI4_new(const Structure &structure);

            void calcERI4Benchmark(const Structure &structure);

            void calcERI4Benchmark_new(const Structure &structure);
        }
    }
}