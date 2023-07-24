#pragma once

#include "geomopt.h"
#include "geometry.h"

namespace Lible
{
    template <GeomOpt::Option option>
    std::vector<double> GeomOpt::optimize(const double &tol_grad_norm, const std::size_t &max_iter,
                                          const std::vector<double> &coords,  const std::vector<std::string> &atoms,
                                          single_point_calc_t singlePointCalc)
    {
        Geometry geometry(coords, atoms);

        std::vector<double> coords_cart = geometry->getCoordsCart();
        for (std::size_t iter = 0; iter < max_iter; iter++)
        {
            auto [energy, grad_cart] = singlePointCalc(coords_cart);

            std::vector<double> coords_redint; // get from coords_cart
            std::vector<double> grad_redint;   // get from

            coords_cart = update<option>(coords_redint, grad_redint);
        }

        return coords_cart;
    }
}