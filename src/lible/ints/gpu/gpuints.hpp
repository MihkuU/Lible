#pragma once

#include <lible/types.hpp>
#include <lible/ints/structure.hpp>

namespace lible
{
    namespace ints
    {
        namespace gpu
        {
            enum class Option
            {
                overlap
            };

            template <Option opt>
            vec2d calculate(const Structure &structure);

            vec2d calculateS_L0(const Structure &structure);

            vec2d calculateS(const Structure &structure);
        }
    }
}