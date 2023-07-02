#include "geomopt.h"

using namespace Lible;
using std::vector;

// template <GeomOpt::Option option>
// vector<double> GeomOpt::optimize(std::function<void(double &energy, vector<double> &geometry, vector<double> &gradient)> singlePointCalculation)
// {
//     vector<double> opt_geometry = geometry;
//     for (std::size_t iter = 0; iter < max_iter; iter++)
//     {
//         double energy;
//         vector<double> geometry, gradient;
//         singlePointCalculation(energy, geometry, gradient);

//         opt_geometry = update<option>(geometry);
//     }

//     return opt_geometry;
// }
