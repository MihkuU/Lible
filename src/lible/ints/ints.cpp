#include <lible/ints/ints.hpp>
#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/twoel/twoel_detail.hpp>

#ifdef _LIBLE_USE_HIP_
#include <lible/ints/gpu/gpuints.hpp>
#endif

namespace LI = lible::ints;
namespace LIO = lible::ints::one;
namespace LIT = lible::ints::two;

lible::vec2d LI::overlap(const Structure &structure)
{
    return LIO::calculate<LIO::Option::overlap>(structure);
}

lible::vec2d LI::overlap_new(const Structure &structure)
{
    return LIO::calculate_new<LIO::Option::overlap>(structure);
}

lible::vec2d LI::kineticEnergy(const Structure &structure)
{
    return LIO::calculate<LIO::Option::kinetic_energy>(structure);
}

lible::vec2d LI::kineticEnergy_new(const Structure &structure)
{
    return LIO::calculate_new<LIO::Option::kinetic_energy>(structure);
}

lible::vec2d LI::nuclearAttraction(const Structure &structure)
{
    return LIO::calculate<LIO::Option::nuclear_attraction>(structure);
}

lible::vec2d LI::nuclearAttraction_new(const Structure &structure)
{
    return LIO::calculate_new<LIO::Option::nuclear_attraction>(structure);
}

lible::vec4d LI::eri4(const Structure &structure)
{
    return LIT::calcERI4(structure);
}

void LI::eri4Benchmark(const Structure &structure)
{
    LIT::calcERI4Benchmark(structure);
}

#ifdef _LIBLE_USE_HIP_
namespace LIG = lible::ints::gpu;

lible::vec2d LIG::overlap0(const Structure &structure)
{
    return LIG::calculateS_L0(structure);
}

lible::vec2d LIG::overlap(const Structure &structure)
{
    return LIG::calculate<LIG::Option::overlap>(structure);
}
#endif