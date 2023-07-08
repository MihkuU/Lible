#include "geomopt.h"

using namespace Lible; 

template<>
std::vector<double> GeomOpt::update<GeomOpt::Option::GDESCENT>(const std::vector<double> &coords_previous)
{
    std::vector<double> coords_new = coords_previous;
    /* Do stuff */
    return coords_new;
}