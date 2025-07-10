#pragma once

#include <lible/types.hpp>
#include <lible/utils.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <chrono>
#include <stdexcept>
#include <string>

#include <armadillo>

namespace lible
{
    namespace ints
    {
        // TODO: this stuff should prolly be removed, its not well generalizeable and 
        // doesnt give much benefit.

        namespace one // TODO: remove one-namespace
        {
            /** */
            enum class Option
            {
                dipole_moment,
                kinetic_energy,
                nuclear_attraction,
                nuclear_attraction_erf,
                overlap
            };

            /** TODO: rename?  */
            template <Option opt>
            void kernel(const int la, const int lb, const ShellPairData &sp_data,
                        vec2d &ints_out); // TODO: remove la, lb

            /** for erf attenuated ints */
            template <Option opt>
            void kernel(const double omega, const int la, const int lb, const ShellPairData &sp_data,
                        vec2d &ints_out); // TODO: remove la, lb

            /** TODO: rename? */
            template <Option opt, typename T>
            void kernel(const int la, const int lb, const ShellPairData &sp_data,
                        const T& arg, std::array<vec2d, 3> &ints_out); // TODO: remove la, lb

            /** For various one-electron integrals. */
            template <Option opt>
            vec2d calculate(const Structure &structure)
            {
                int l_max = structure.getMaxL();
                size_t dim_ao = structure.getDimAO();

                vec2d ints(Fill(0), dim_ao, dim_ao);
                for (int la = l_max; la >= 0; la--)
                {
                    ShellPairData sp_data = shellPairDataSymm(la, la, structure);

                    kernel<opt>(la, la, sp_data, ints);
                }

                for (int la = l_max; la >= 0; la--)
                    for (int lb = la - 1; lb >= 0; lb--)
                    {
                        ShellPairData sp_data = shellPairDataSymm(la, lb, structure);

                        kernel<opt>(la, lb, sp_data, ints);
                    }

                return ints;
            }

            /** For the erf-attenuated nuclear attraction integrals. */
            template <Option opt>
            vec2d calculate(const Structure &structure, const double omega)
            {
                int l_max = structure.getMaxL();
                size_t dim_ao = structure.getDimAO();

                vec2d ints(Fill(0), dim_ao, dim_ao);
                for (int la = l_max; la >= 0; la--)
                {
                    ShellPairData sp_data = shellPairDataSymm(la, la, structure);

                    kernel<opt>(omega, la, la, sp_data, ints);
                }

                for (int la = l_max; la >= 0; la--)
                    for (int lb = la - 1; lb >= 0; lb--)
                    {
                        ShellPairData sp_data = shellPairDataSymm(la, lb, structure);

                        kernel<opt>(omega, la, lb, sp_data, ints);
                    }

                return ints;
            }

            /** For dipole mom, linear mom and angmom. */
            template <Option opt, typename T>
            std::array<vec2d, 3> calculate3D(const Structure &structure, const T& arg)
            {
                int l_max = structure.getMaxL();
                size_t dim_ao = structure.getDimAO();

                std::array<vec2d, 3> ints{vec2d(Fill(0), dim_ao, dim_ao),
                                          vec2d(Fill(0), dim_ao, dim_ao),
                                          vec2d(Fill(0), dim_ao, dim_ao)};

                for (int la = l_max; la >= 0; la--)
                {
                    ShellPairData sp_data = shellPairDataSymm(la, la, structure);

                    kernel<opt, T>(la, la, sp_data, arg, ints);
                }

                for (int la = l_max; la >= 0; la--)
                    for (int lb = la - 1; lb >= 0; lb--)
                    {
                        ShellPairData sp_data = shellPairDataSymm(la, lb, structure);

                        kernel<opt, T>(la, lb, sp_data, arg, ints);
                    }

                return ints;
            }
        }
    }
}