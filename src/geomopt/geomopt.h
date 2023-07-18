#pragma once

#include <functional>
#include <memory>
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

        // TODO: consider removing the 'class'
        enum class Option
        {
            BFGS,
            GDESCENT,
            GDIIS,
            KRIGING
        };

        // TODO: consider removing the 'class'
        enum class OptionForCIOpt
        {
            IGNACIO,
            MOROKUMA,
            YARKONY           
        };

        GeomOpt();

        template <Option option>
        std::vector<double> optimize(std::function<void(const std::vector<double> &coords_cart, double &energy, std::vector<double> &gradient_cart)> singlePointCalculation);

        template <OptionForCIOpt option_for_ciopt>
        std::vector<double> optimizeConicalIntersection();

    private:
        struct Geometry;
        std::unique_ptr<Geometry> geometry;

        template <Option option>
        std::vector<double> update(const std::vector<double> &coords_redint, const std::vector<double> &grad_redint);

        double tol_grad_norm{1e-3};
        double tol_energy_diff{1e-6};
        std::size_t max_iter{1};
    };
}