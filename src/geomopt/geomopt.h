#pragma once

#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace Lible
{
    namespace GeomOpt
    {
        enum Option
        {
            BFGS,
            GDESCENT,
            GDIIS,
            KRIGING
        };

        enum OptionConicalIntersection
        {
            IGNACIO,
            MOROKUMA,
            YARKONY
        };

        typedef std::function<std::vector<double>(const std::vector<double> &coords_bohr)> single_point_calc_t;

        template <Option option>
        std::vector<double> optimize(const double &tol_grad_norm, const std::size_t &max_iter,
                                     const std::vector<double> &coords_bohr, 
                                     const std::vector<std::string> &atomic_symbols,
                                     single_point_calc_t singlePointCalc);

        // template <Option option>
        // std::vector<double> optimizeConstrained();

        // template <Option option> 
        // std::vector<double> optimizeTransitionState();

        // template <Option option>
        // std::vector<std::vector<double>> scanRelaxed();

        // template <Option option>
        // std::vector<std::vector<double>> scanUnelaxed();                

        template <OptionConicalIntersection option_ci>
        std::vector<double> optimizeConicalIntersection();

        template <Option option>
        struct Optimizer;
    }
}