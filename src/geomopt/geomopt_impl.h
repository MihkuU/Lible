#pragma once 

#include "geomopt.h"

namespace Lible
{
    template <GeomOpt::Option option>
    std::vector<double> GeomOpt::optimize(std::function<
                                 void(double &energy, std::vector<double> &geometry, std::vector<double> &gradient)>
                                     singlePointCalculation)
    {
        std::vector<double> opt_geometry = geometry;
        for (std::size_t iter = 0; iter < max_iter; iter++)
        {
            double energy;
            std::vector<double> geometry, gradient;
            singlePointCalculation(energy, geometry, gradient);

            opt_geometry = update<option>(geometry);
        }

        return opt_geometry;
    }
}