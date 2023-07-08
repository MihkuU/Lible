#include "geomopt.h"

using namespace Lible;

template<>
std::vector<double> GeomOpt::update<GeomOpt::Option::GDIIS>(const std::vector<double> &coords_prev)
{
    std::vector<double> coords_new = coords_prev;
    /* Do stuff */
    return coords_new;    
}