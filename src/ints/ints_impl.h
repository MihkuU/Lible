#pragma once

#include "kernels_1el_ints.h"
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
                Kernels1El::oneElIntKernel<option>(la, lb, structure->shell_pairs.at(std::make_pair(la, lb)),
                                                   one_el_ints);
        }

        return one_el_ints;
    }
}