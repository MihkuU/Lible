#pragma once

#include "ints.h"
#include "shell.h"

namespace Lible
{
    namespace Kernels1El
    {
        template <Ints::Option1El option>
        void oneElIntKernel(const Shells::ShellPair &shell_pair, std::vector<double> &one_el_ints);

        namespace ObaraSaika
        {
            void overlapX(const int &angmom_a, const int &angmom_b, const double &one_over_2p,
                          const double &x_pa, const double &x_pb, const double &overlap_00,
                          std::vector<double> &overlap_x);
        }
    }
}