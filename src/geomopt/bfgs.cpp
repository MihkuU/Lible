#include "geomopt.h"

using namespace Lible; 

template<>
std::vector<double> GeomOpt::update<GeomOpt::Option::BFGS>(const std::vector<double> &coords_redint, const std::vector<double> &grad_redint)
{
    std::vector<double> coords_new = coords_redint;
    /* Do stuff */
    return coords_new;
}