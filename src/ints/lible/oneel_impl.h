#pragma once

#include <lible/ints.h>
#include <lible/types.h>
#include <lible/util.h>

#include <chrono>
#include <stdexcept>
#include <string>

#include <fmt/core.h>

namespace lible
{
    namespace ints
    {
        namespace one
        {
            enum class Option
            {
                DIPOLE,
                KINETIC,
                NUCLEAR,
                OVERLAP
            };

            static std::string returnPreamble(const Option &opt)
            {
                switch (opt)
                {
                case Option::DIPOLE:
                    return "Dipole moment integrals...";
                case Option::KINETIC:
                    return "Kinetic energy integrals...";
                case Option::NUCLEAR:
                    return "Nuclear attraction integrals...";
                case Option::OVERLAP:
                    return "Overlap integrals...";
                default:
                    throw std::runtime_error("Inappropriate one-electron integral option.\n");
                }
            }

            template <Option opt>
            vec2d calc(const Structure &structure)
            {
                auto start{std::chrono::steady_clock::now()};

                std::string msg = returnPreamble(opt);
                palPrint(fmt::format("Lible::{:<40}", msg));

                int l_max = structure.max_angular_momentum;
                size_t n_ao = structure.n_atomic_orbitals;

                vec2d ints(n_ao, n_ao, 0);
                for (int la = l_max; la >= 0; la--)
                {
                    /*
                     * Stage the shell-pairs/calculate Hermite exp. coeffs
                     */
                    continue;
                }

                for (int la = l_max; la >= 0; la--)
                    for (int lb = la - 1; lb >= 0; lb--)
                    {
                        /*
                         * Stage the shell-pairs/calculate Hermite exp. coeffs
                         */
                        continue;
                    }

                auto end(std::chrono::steady_clock::now());
                std::chrono::duration<double> duration{end - start};                
                palPrint(fmt::format(" {:.2e} s\n", duration.count()));

                return ints;
            }

            template <Option opt>
            void kernel()
            {
            }
        }
    }
}