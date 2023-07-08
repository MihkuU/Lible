#pragma once

#include "geomopt.h"

namespace Lible
{
    template <GeomOpt::Option option>
    std::vector<double> GeomOpt::optimize(std::function<void(double &energy, std::vector<double> &coords, std::vector<double> &gradient)> singlePointCalculation)
    {
        // std::vector<double> coords_opt = coords;
        // for (std::size_t iter = 0; iter < max_iter; iter++)
        // {
        //     double energy;
        //     std::vector<double> geometry, gradient;
        //     singlePointCalculation(energy, geometry, gradient);

        //     coords_opt = update<option>(coords_opt);
        // }

        // return coords_opt;
    }
}