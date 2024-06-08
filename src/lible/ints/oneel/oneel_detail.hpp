#pragma once

#include <lible/ints/mcmurchie_davidson.hpp>
#include <lible/ints/shell_pair_data.hpp>
#include <lible/types.hpp>
#include <lible/util.hpp>

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

            // template <Option opt>
            // void kernel(const ShellPairData &shell_pair_data,
            //             const std::vector<std::vector<arma::dmat>> &ecoeffs, vec2d &ints_out);

            // template <Option opt>
            // std::vector<std::vector<arma::dmat>> calcECoeffs(const int la, const int lb,
            //                                                  const ShellPairData &shell_pair_data);

            template <Option opt>
            void kernel(const int la, const int lb, const ShellPairData &shell_pair_data,
                        vec2d &ints_out);

            template <Option opt>
            vec2d calculate(const Structure &structure)
            {
                auto start{std::chrono::steady_clock::now()};

                std::string msg = returnPreamble(opt);
                palPrint(fmt::format("Lible::{:<40}", msg));

                int l_max = structure.getMaxL();
                size_t dim_ao = structure.getDimAO();

                vec2d ints(dim_ao, dim_ao, 0);
                for (int la = l_max; la >= 0; la--)
                {
                    ShellPairData shell_pair_data = ShellPairData(la, la, structure);

                    kernel<opt>(la, la, shell_pair_data, ints);
                }

                for (int la = l_max; la >= 0; la--)
                    for (int lb = la - 1; lb >= 0; lb--)
                    {
                        ShellPairData shell_pair_data = ShellPairData(la, lb, structure);

                        kernel<opt>(la, lb, shell_pair_data, ints);
                    }

                auto end{std::chrono::steady_clock::now()};
                std::chrono::duration<double> duration{end - start};
                palPrint(fmt::format(" {:.2e} s\n", duration.count()));

                return ints;
            }
        }
    }
}