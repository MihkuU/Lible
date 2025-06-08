#include <lible/ints/ints.hpp>
#include <lible/ints/oneel/oneel_detail.hpp>

#ifdef _LIBLE_USE_HIP_
#include <lible/ints/gpu/gpuints.hpp>
#endif

namespace LI = lible::ints;
namespace LIO = lible::ints::one;

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

std::array<lible::vec2d, 3> LI::dipoleMoment(const std::array<double, 3> &origin,
                                             const Structure &structure)
{
    return LIO::calculate3D<LIO::Option::dipole_moment>(structure, origin);
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
