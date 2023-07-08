#include "structure.h"
#include "ints.h"

#include "spherical_trafo.h"

using namespace Lible;

Ints::Ints(const std::string &basis_set, const std::vector<double> &coordinates, const std::vector<std::string> &elements)
{
    assert(coordinates.size() / 3 == elements.size());    
    
    structure.reset(new Structure(basis_set, coordinates, elements));
}

template <Ints::Option2El option>
std::vector<double> Ints::calcTwoElInts(const std::vector<double> &density)
{
}