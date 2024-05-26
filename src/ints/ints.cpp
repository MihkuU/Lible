#include <lible/ints.hpp>
#include <lible/oneel_detail.hpp>
#include <lible/twoel_detail.hpp>

namespace LI = lible::ints;
namespace LIO = lible::ints::one;
namespace LIT = lible::ints::two;

using namespace lible;

vec2d LI::overlap(const Structure &structure)
{
    return LIO::calculate<LIO::Option::overlap>(structure);
}

vec2d LI::kineticEnergy(const Structure &structure)
{
    return LIO::calculate<LIO::Option::kinetic_energy>(structure);
}

vec2d LI::nuclearAttraction(const Structure &structure)
{
    return LIO::calculate<LIO::Option::nuclear_attraction>(structure);
}

vec4d LI::eri4(const Structure &structure)
{
    return LIT::calcERI4(structure);
}