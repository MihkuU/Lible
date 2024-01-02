#pragma once

#include "geomopt.h"
#include "geometry.h"
#include "geomopt_utils.h"

namespace lible
{
    namespace geomopt
    {
        template <>
        struct Optimizer<BFGS>
        {
            std::vector<double> update(const std::vector<double> &coords_redint_in,
                                       const std::vector<double> &grad_redint);

        private:
        };

        template <>
        struct Optimizer<GDESCENT>
        {
            std::vector<double> update(const std::vector<double> &coords,
                                       const std::vector<double> &grad,
                                       const size_t &iter);

        private:
            std::vector<double> grad_prev;
            std::vector<double> step_prev;

            enum Method
            {
                FletcherReeves,
                PolakRibiere
            };

            inline static int method = FletcherReeves;
        };

        template <Option option>
        std::vector<double> optimize(const double &tol_grad_norm, const std::size_t &max_iter,
                                     const std::vector<double> &coords_bohr,
                                     const std::vector<std::string> &atomic_symbols,
                                     single_point_calc_t singlePointCalc)
        {
            Geometry geometry(coords_bohr, atomic_symbols);
            optimizePrintPreamble(geometry);

            Optimizer<option> optimizer;

            bool converged = false;
            std::vector<double> coords_cart = coords_bohr;
            std::vector<double> coords_redint = geometry.getCoordsRedint();
            for (std::size_t iter = 0; iter < max_iter; iter++)
            {
                std::vector<double> grad_cart = singlePointCalc(coords_cart);
                std::vector<double> grad_redint = geometry.transformGradCartToRedint(grad_cart);

                double grad_norm = calcGradNorm(grad_redint);
                if (grad_norm < tol_grad_norm)
                {
                    converged = true;
                    break;
                }

                std::vector<double> step_redint = optimizer.update(coords_redint, grad_redint, iter);
                std::vector<double> step_cart = geometry.transformStepRedIntToCart(coords_cart, coords_redint,
                                                                                   step_redint);

                std::transform(coords_redint.begin(), coords_redint.end(),
                               step_redint.begin(), coords_redint.begin(),
                               std::plus<double>());

                std::transform(coords_cart.begin(), coords_cart.end(),
                               step_cart.begin(), coords_cart.begin(),
                               std::plus<double>());

                optimizePrintIter<option>(iter);
            }
            optimizePrintEpilogue(converged);

            return coords_cart;
        }
    }
}