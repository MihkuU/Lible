#include "ints.h"

#include "structure.h"
#include "experimental/structure_gpu.h"

#include "spherical_trafo.h"

using namespace lible;

// ints::ints(const std::string &basis_set, const std::vector<double> &coordinates_ang,
//            const std::vector<std::string> &elements)
// {
//     assert(coordinates_ang.size() / 3 == elements.size());    
    
//     structure.reset(new Structure(basis_set, coordinates_ang, elements));

//     structure_gpu.reset(new StructureGPU(structure));
// }

// template <ints::Option2El option>
// std::vector<double> ints::calcTwoElInts(const std::vector<double> &density)
// {
// }