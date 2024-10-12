#pragma once

#include <lible/log.hpp>
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

            static std::string returnPreamble(const Option &opt)
            {
                switch (opt)
                {
                case Option::dipole_moment:
                    return "Dipole moment integrals...";
                case Option::kinetic_energy:
                    return "Kinetic energy integrals...";
                case Option::nuclear_attraction:
                    return "Nuclear attraction integrals...";
                case Option::overlap:
                    return "Overlap integrals...";
                default:
                    throw std::runtime_error("Inappropriate one-electron integral option.\n");
                }
            }

            template <Option opt>
            void kernel(const int la, const int lb, const ShellPairData &sp_data,
                        vec2d &ints_out);

            template <Option opt>
            vec2d calculate(const Structure &structure)
            {
                auto start{std::chrono::steady_clock::now()};

                std::string msg = returnPreamble(opt);            
                log::logger << fmt::format("Lible::{:<40}", msg);

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

                auto end{std::chrono::steady_clock::now()};
                std::chrono::duration<double> duration{end - start};                
                
                log::logger << fmt::format(" {:.2e} s\n", duration.count());

                return ints;
            }
        }
    }
}