#pragma once

#include <functional>
#include <string>
#include <vector>

namespace Lible
{
    /*
     * A class for geometry optimization methods
     *
     * Abbreviations:
     *   BFGS     -
     *   DIIS     -
     *   GDESCENT -
     *   GDIIS    - generalized DIIS
     *   GEOM     - geometry
     *   Opt      - optimization
     */
    class GeomOpt
    {
    public:
        enum class Option
        {
            BFGS,
            GDESCENT,
            GDIIS,
            KRIGING
        };

        enum class OptionForCIOpt
        {
            IGNACIO,
            YARKONY
        };

        GeomOpt(const std::vector<double> geometry_) : geometry{geometry_}
        {
        }

        template <Option option>
        std::vector<double> optimize(std::function<
                                     void(double &energy, std::vector<double> &geometry, std::vector<double> &gradient)>
                                         singlePointCalculation);
        // {
        //     std::vector<double> opt_geometry = geometry;
        //     for (std::size_t iter = 0; iter < max_iter; iter++)
        //     {
        //         double energy;
        //         std::vector<double> geometry, gradient;
        //         singlePointCalculation(energy, geometry, gradient);

        //         opt_geometry = update<option>(geometry);
        //     }

        //     return opt_geometry;
        // }

        template <OptionForCIOpt option_for_ciopt>
        std::vector<double> optimizeConicalIntersection();

    private:
        template <Option option>
        std::vector<double> update(const std::vector<double> &previous_geometry);

        double tol_grad_norm{1e-3};
        double tol_energy_diff{1e-6};
        std::size_t max_iter{1};
        std::vector<double> geometry{};
    };
}