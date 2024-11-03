#include <lible/ints/ints.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/oneel/oneel_detail.hpp>
#include <lible/ints/twoel/twoel_detail.hpp>

#ifdef _LIBLE_USE_HIP_
#include <lible/ints/gpu/gpuints.hpp>
#endif

namespace LI = lible::ints;
namespace LIO = lible::ints::one;
namespace LIT = lible::ints::two;

std::vector<double> LI::eri2Diagonal(const Structure &structure)
{
    return LIT::calcERI2Diagonal(structure);
}

lible::vec2d LI::eri4Diagonal(const Structure &structure)
{
    return LIT::calcERI4Diagonal(structure);
}

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

lible::vec2d LI::eri2(const Structure &structure)
{
    return LIT::calcERI2(structure);
}

lible::vec3d LI::eri3(const Structure &structure)
{
    return LIT::calcERI3(structure);
}

lible::vec4d LI::eri4(const Structure &structure)
{
    return LIT::calcERI4(structure);
}

void LI::eri4Benchmark(const Structure &structure)
{
    LIT::calcERI4Benchmark(structure);
}

void LI::eri4BenchmarkNew(const Structure &structure)
{
    LIT::calcERI4BenchmarkNew(structure);
}

lible::vec4d LI::eri4New(const Structure &structure)
{
    return LIT::calcERI4New(structure);
}

LI::kernel_eri4_t LI::deployERI4Kernel(const int la, const int lb, const int lc, const int ld)
{
    return LIT::deployERI4Kernel(la, lb, lc, ld);
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