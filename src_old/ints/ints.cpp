#include <lible/ints.hpp>
#include <lible/oneel_detail.hpp>
#include <lible/twoel_detail.hpp>

#ifdef _LIBLE_USE_HIP_
#include <lible/oneel_detail_gpu.hpp>
#endif 

namespace LI = lible::ints;
namespace LIO = lible::ints::one;
namespace LIT = lible::ints::two;

lible::vec2d LI::overlap(const Structure &structure)
{
    return LIO::calculate<LIO::Option::overlap>(structure);
}

lible::vec2d LI::kineticEnergy(const Structure &structure)
{
    return LIO::calculate<LIO::Option::kinetic_energy>(structure);
}

lible::vec2d LI::nuclearAttraction(const Structure &structure)
{
    return LIO::calculate<LIO::Option::nuclear_attraction>(structure);
}

lible::vec4d LI::eri4(const Structure &structure)
{
    return LIT::calcERI4(structure);
}

lible::vec4d LI::eri4Shark(const Structure &structure)
{
    return LIT::calcERI4Shark(structure);
}

void LI::eri4Benchmark(const Structure &structure)
{
    LIT::calcERI4Benchmark(structure);
}

#ifdef _LIBLE_USE_HIP_
namespace LIG = lible::ints::gpu;
namespace LIOG = lible::ints::one_gpu;

lible::vec2d LIG::overlap(const Structure &structure)
{
    return LIOG::calculateS_L0(structure);
}
#endif