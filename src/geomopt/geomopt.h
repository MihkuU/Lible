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

        typedef std::function<std::tuple<double, std::vector<double>>(const std::vector<double> &coords)> single_point_calc_t;

        template <Option option>
        std::vector<double> optimize(const double &tol_grad_norm, const std::size_t &max_iter,
                                     const std::vector<double> &coords, const std::vector<std::string> &atoms,
                                     single_point_calc_t singlePointCalc);

        template <OptionConicalIntersection option_ci>
        std::vector<double> optimizeConicalIntersection();

        template <Option option>
        std::vector<double> update(const std::vector<double> &coords_redint, const std::vector<double> &grad_redint);
    }

    // /*
    //  * A class for geometry optimization methods
    //  *
    //  * Abbreviations:
    //  *   BFGS     -
    //  *   DIIS     -
    //  *   GDESCENT -
    //  *   GDIIS    - generalized DIIS
    //  *   GEOM     - geometry
    //  *   Opt      - optimization
    //  */

    // class GeomOpt
    // {
    // public:
    //     GeomOpt(const double &tol_grad_norm_, const std::size_t &max_iter_);

    //     enum Option
    //     {
    //         BFGS,
    //         GDESCENT,
    //         GDIIS,
    //         KRIGING
    //     };

    //     enum OptionCIOpt
    //     {
    //         IGNACIO,
    //         MOROKUMA,
    //         YARKONY
    //     };

    //     // typedef std::function<void(const std::vector<double> &coords, double &energy, std::vector<double> &grad)> single_point_calc_t;
    //     typedef std::function<std::tuple<double, std::vector<double>>(const std::vector<double> &coords)> single_point_calc_t;

    //     template <Option option>
    //     std::vector<double> optimize(single_point_calc_t singlePointCalc);

    //     template <OptionCIOpt option_ciopt>
    //     std::vector<double> optimizeConicalIntersection();

    //     double getTolGradNorm() const
    //     {
    //         return tol_grad_norm;
    //     }

    //     std::size_t getMaxIter() const
    //     {
    //         return max_iter;
    //     }

    // private:
    //     struct Geometry;
    //     std::unique_ptr<Geometry> geometry;

    //     template <Option option>
    //     std::vector<double> update(const std::vector<double> &coords_redint, const std::vector<double> &grad_redint);

    //     double tol_grad_norm{1e-3};
    //     std::size_t max_iter{1};
    // };
}