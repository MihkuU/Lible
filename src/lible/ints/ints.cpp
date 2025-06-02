#include <lible/ints/ints.hpp>
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

std::array<lible::vec2d, 3> LI::dipoleMoment(const std::array<double, 3> &origin,
                                             const Structure &structure)
{
    return LIO::calculate3D<LIO::Option::dipole_moment>(structure, origin);
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

void LI::eri4BenchmarkTest(const Structure &structure)
{
    LIT::calcERI4BenchmarkTest(structure);
}

LI::kernel_eri4_t LI::deployERI4Kernel(const int la, const int lb, const int lc, const int ld)
{
    return LIT::deployERI4Kernel(la, lb, lc, ld);
}

LI::kernel_eri3_t LI::deployERI3Kernel(const int la, const int lb, const int lc)
{
    return LIT::deployERI3Kernel(la, lb, lc);
}

LI::kernel_eri2_t LI::deployERI2Kernel(const int la, const int lb)
{
    return LIT::deployERI2Kernel(la, lb);
}

std::array<lible::vec2d, 6> LI::kernelERI2Deriv1(const int ishell_a, const int ishell_b,
                                                 const std::vector<double> &ecoeffs_a,
                                                 const std::vector<double> &ecoeffs_b_tsp,
                                                 const BoysGrid &boys_grid,
                                                 const ShellData &sh_data_a,
                                                 const ShellData &sh_data_b)
{
    return LIT::kernelERI2Deriv1(ishell_a, ishell_b, ecoeffs_a, ecoeffs_b_tsp, boys_grid,
                                 sh_data_a, sh_data_b);
}

std::array<lible::vec3d, 9> LI::kernelERI3Deriv1(const int ipair_ab, const int ishell_c,
                                                 const std::vector<double> &ecoeffs_ab,
                                                 const std::vector<double> &ecoeffs1_ab,
                                                 const std::vector<double> &ecoeffs_c,
                                                 const BoysGrid &boys_grid,
                                                 const ShellPairData &sp_data_ab,
                                                 const ShellData &sh_data_c)
{
    return LIT::kernelERI3Deriv1(ipair_ab, ishell_c, ecoeffs_ab, ecoeffs1_ab, ecoeffs_c, boys_grid,
                                 sp_data_ab, sh_data_c);
}

std::array<lible::vec3d, 9> LI::kernelERI3Deriv1Test(const int ipair_ab, const int ishell_c,
                                              const std::vector<double> &ecoeffs0_bra,
                                              const std::vector<double> &ecoeffs1_bra,
                                              const std::vector<double> &ecoeffs0_ket,
                                              const BoysGrid &boys_grid,
                                              const ShellPairData &sp_data_ab,
                                              const ShellData &sh_data_c)
{
    return LIT::kernelERI3Deriv1Test(ipair_ab, ishell_c, ecoeffs0_bra, ecoeffs1_bra, ecoeffs0_ket,
                                     boys_grid, sp_data_ab, sh_data_c);
}

std::array<lible::vec4d, 12> LI::kernelERI4Deriv1(const int ipair_ab, const int ipair_cd,
                                                  const std::vector<double> &ecoeffs_ab,
                                                  const std::vector<double> &ecoeffs1_ab,
                                                  const std::vector<double> &ecoeffs_cd_tsp,
                                                  const std::vector<double> &ecoeffs1_cd_tsp,
                                                  const BoysGrid &boys_grid,
                                                  const ShellPairData &sp_data_ab,
                                                  const ShellPairData &sp_data_cd)
{
    return LIT::kernelERI4Deriv1(ipair_ab, ipair_cd, ecoeffs_ab, ecoeffs1_ab, ecoeffs_cd_tsp,
                                 ecoeffs1_cd_tsp, boys_grid, sp_data_ab, sp_data_cd);
}

std::array<lible::vec4d, 12> LI::kernelERI4Deriv1Test(const int ipair_ab, const int ipair_cd,
                                                      const std::vector<double> &ecoeffs0_bra,
                                                      const std::vector<double> &ecoeffs1_bra,
                                                      const std::vector<double> &ecoeffs0_ket,
                                                      const std::vector<double> &ecoeffs1_ket,
                                                      const BoysGrid &boys_grid,
                                                      const ShellPairData &sp_data_ab,
                                                      const ShellPairData &sp_data_cd)
{
    return LIT::kernelERI4Deriv1Test(ipair_ab, ipair_cd, ecoeffs0_bra, ecoeffs1_bra, ecoeffs0_ket,
                                     ecoeffs1_ket, boys_grid, sp_data_ab, sp_data_cd);
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