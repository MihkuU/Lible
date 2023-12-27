#include "structure.h"
#include "ints.h"
#include "experimental/structure_gpu.h"

#include "spherical_trafo.h"

using namespace Lible;

// Ints::Ints(const std::string &basis_set, const std::vector<double> &coordinates_ang,
//            const std::vector<std::string> &elements)
// {
//     assert(coordinates_ang.size() / 3 == elements.size());    
    
//     structure.reset(new Structure(basis_set, coordinates_ang, elements));

//     structure_gpu.reset(new StructureGPU(structure));
// }

// template <Ints::Option2El option>
// std::vector<double> Ints::calcTwoElInts(const std::vector<double> &density)
// {
// }