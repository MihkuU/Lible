#include "geomopt.h"
using namespace Lible;

template<>
std::vector<double> GeomOpt::update<GeomOpt::Option::KRIGING>(const std::vector<double> &previous_geometry)
{
    std::vector<double> new_geometry = previous_geometry;
    /* Do stuff */
    return new_geometry;
}