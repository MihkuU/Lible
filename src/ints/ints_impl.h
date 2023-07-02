#pragma once

#include "kernel_1el_ints.h"
#include "spherical_trafo.h"
#include "structure.h"
#include "ints.h"

namespace Lible
{
    template <Ints::Option1El option>
    std::vector<double> Ints::calcOneElInts()
    {
        int max_l = structure->max_angular_momentum;

        std::vector<double> one_el_ints;
        for (int la = 0; la < max_l; la++)
        {
            for (int lb = 0; lb <= la; lb++)
            {                
                std::size_t dim_cart_a = Shells::calcShellDimCartesian(la);
                std::size_t dim_cart_b = Shells::calcShellDimCartesian(lb);
                std::vector<double> one_el_ints_cart(dim_cart_a * dim_cart_b, 0);
                // This design already starts to break down, looks like the following should all go inside the kernel-function
                for (const Shells::ShellPair &shell_pair : structure->shell_pairs[std::make_pair(la, lb)])
                {
                    Kernels1El::oneElIntKernel<option>(shell_pair, one_el_ints_cart);
                    SphericalTrafo::transformAndTransferInts(one_el_ints_cart, shell_pair, one_el_ints);
                }                
            }
        }

        return one_el_ints;
    }
}