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

        GeomOpt();

        template <Option option>
        std::vector<double> optimize(std::function<void(double &energy, std::vector<double> &coords, std::vector<double> &gradient)> singlePointCalculation);

        template <OptionForCIOpt option_for_ciopt>
        std::vector<double> optimizeConicalIntersection();

    private:
        struct Geometry;
        std::unique_ptr<Geometry> geometry;

        template <Option option>
        std::vector<double> update(const std::vector<double> &coords_previous);

        double tol_grad_norm{1e-3};
        double tol_energy_diff{1e-6};
        std::size_t max_iter{1};
    };
}