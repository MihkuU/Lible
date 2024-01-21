#include <lible/ints.h>
#include <lible/oneel_impl.h>
#include <lible/structure.h>

namespace LI = lible::ints;
namespace LIO = lible::ints::one;

using namespace lible;

vec2d LI::overlap(const Structure &structure)
{
    return LIO::calc<LIO::Option::OVERLAP>(structure);
}

vec2d LI::kineticEnergy(const Structure &structure)
{
}

vec2d LI::coulombAttraction(const Structure &structure)
{
}

vec2d LI::dipoleMoment(const Structure &structure)
{
}

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