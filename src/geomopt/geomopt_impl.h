#pragma once

#include "geomopt.h"
#include "geometry.h"

namespace Lible
{
    template <GeomOpt::Option option>
    std::vector<double> GeomOpt::optimize(std::function<void(const std::vector<double> &coords_cart, double &energy, std::vector<double> &gradient_cart)> singlePointCalculation)
    {
        std::vector<double> coords_cart = geometry->getCoordsCart();
        for (std::size_t iter = 0; iter < max_iter; iter++)
        {
            double energy;
            std::vector<double> grad_cart;
            singlePointCalculation(coords_cart, energy, grad_cart);

            std::vector<double> coords_redint;
            std::vector<double> grad_redint;

            coords_cart = update<option>(coords_redint, grad_redint);
        }

        return coords_cart;
    }
}