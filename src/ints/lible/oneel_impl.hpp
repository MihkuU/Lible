#pragma once

#include <lible/ints.hpp>
#include <lible/mcmurchie_davidson.hpp>
#include <lible/shell_pair_data.hpp>
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
            void kernel(const ShellPairData &shell_pair_data,
                        const std::vector<std::vector<arma::dmat>> &h_coeffs,
                        vec2d &ints_out);

            template <Option opt>
            vec2d calc(const Structure &structure)
            {
                auto start{std::chrono::steady_clock::now()};

                // TODO: log instead of print, fcks sake
                std::string msg = returnPreamble(opt);
                palPrint(fmt::format("Lible::{:<40}", msg));

                int l_max = structure.max_angular_momentum;
                size_t n_ao = structure.n_atomic_orbitals;

                vec2d ints(n_ao, n_ao, 0);
                for (int la = l_max; la >= 0; la--)
                {
                    auto shell_pair_data = ShellPairData(la, la, structure);

                    std::vector<std::vector<arma::dmat>> h_coeffs;
                    MD::calcHCoeffs(la, la, shell_pair_data, h_coeffs);
                                        
                    kernel<opt>(shell_pair_data, h_coeffs, ints);
                }

                for (int la = l_max; la >= 0; la--)
                    for (int lb = la - 1; lb >= 0; lb--)
                    {
                        auto shell_pair_data = ShellPairData(la, lb, structure);

                        std::vector<std::vector<arma::dmat>> h_coeffs;
                        MD::calcHCoeffs(la, lb, shell_pair_data, h_coeffs);

                        kernel<opt>(shell_pair_data, h_coeffs, ints);                        
                    }

                auto end(std::chrono::steady_clock::now());
                std::chrono::duration<double> duration{end - start};                
                palPrint(fmt::format(" {:.2e} s\n", duration.count()));

                return ints;
            }
        }
    }
}