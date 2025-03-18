#pragma once

#include <lible/types.hpp>
#include <lible/utils.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <chrono>
#include <stdexcept>
#include <string>

#include <armadillo>
#include <fmt/core.h>

namespace lible
{
    namespace ints
    {
        namespace one
        {
            /** */
            enum class Option
            {
                dipole_moment,
                kinetic_energy,
                nuclear_attraction,
                overlap
            };

            template <Option opt>
            void kernel(const int la, const int lb, const ShellPairData &sp_data,
                        vec2d &ints_out);

            template <Option opt>
            vec2d calculate(const Structure &structure)
            {                
                int l_max = structure.getMaxL();
                size_t dim_ao = structure.getDimAO();

                vec2d ints(dim_ao, dim_ao, 0);
                for (int la = l_max; la >= 0; la--)
                {
                    ShellPairData sp_data = constructShellPairData(la, la, structure);

                    kernel<opt>(la, la, sp_data, ints);
                }

                for (int la = l_max; la >= 0; la--)
                    for (int lb = la - 1; lb >= 0; lb--)
                    {
                        ShellPairData sp_data = constructShellPairData(la, lb, structure);

                        kernel<opt>(la, lb, sp_data, ints);
                    }
                
                return ints;
            }
        }
    }
}