#include "defs_geomopt.h"
#include "geomopt_utils.h"

#include <cmath>
#include <fmt/core.h>
#include <iostream>

namespace LG = Lible::GeomOpt;
using namespace LG;

using std::string;
using std::vector;
using fmt::format;

double LG::calcGradNorm(const vector<double> &grad)
{
    double norm = 0;
    for (const double &val : grad)
        norm = pow(val, 2);
    norm = sqrt(norm);

    return norm;
}

void LG::optimizePrintPreamble(const Geometry &geometry)
{
    std::cout << string(Lible::GeomOptDefs::n_print_chars, '-') << std::endl;
    std::cout << string(3, ' ') << "Lible-GeomOpt optimize called" << std::endl;
    std::cout << string(3, ' ') << "#-redundant internal coordinates:" << std::endl; //TODO: replace # with the actual number of coords.
    std::cout << string(6, ' ') << format("Bonds: {}", geometry.getNumBonds()) << std::endl; //TODO: make it like n-bonds if # > 0.
    std::cout << string(6, ' ') << format("Angles: {}", geometry.getNumAngles()) << std::endl; //TODO: make it like m-angles if # > 0.
    std::cout << string(6, ' ') << format("Dihedrals: {}", geometry.getNumDihedrals()) << std::endl; //TODO: make it like o-dihedrals if # > 0.

    // std::string 
    // std::cout << fmt::format() << std::endl;

    std::cout << std::string(Lible::GeomOptDefs::n_print_chars, '-') << std::endl;
}

void LG::optimizePrintEpilogue(const bool &converged)
{

}

template<>
void LG::optimizePrintIter<GDESCENT>(const size_t &iter)
{
    std::cout << "Hi, I am optimizing your geometry using gradient descent." << std::endl;
    std::cout << "Hi, I am optimizing your geometry using gradient descent." << std::endl;
    std::cout << "Hi, I am optimizing your geometry using gradient descent." << std::endl;
    std::cout << "Hi, I am optimizing your geometry using gradient descent." << std::endl;
}

template<>
void LG::optimizePrintIter<BFGS>(const size_t &iter)
{
    
}