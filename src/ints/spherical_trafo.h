#pragma once

#include "shell.h"

namespace Lible
{
    namespace SphericalTrafo
    {
        void transformAndTransferInts(const std::vector<double> &one_el_ints_cart, const Shells::ShellPair &shell_pair,
                                      std::vector<double> &one_el_ints);
    }
}