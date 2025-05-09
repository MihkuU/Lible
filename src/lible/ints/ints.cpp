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

lible::vec4d LI::eri4New(const Structure &structure)
{
    return LIT::calcERI4New(structure);
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

void LI::kernelERI2Deriv1(const int la, const int lb, const int cdepth_a, const int cdepth_b,
                          const double *exps_a, const double *exps_b, const double *coords_a,
                          const double *coords_b, const double *ecoeffs_a,
                          const double *ecoeffs_b_tsp, const double *norms_a,
                          const double *norms_b, const BoysGrid &boys_grid,
                          double *eri2_batch)
{
    LIT::kernelERI2Deriv1(la, lb, cdepth_a, cdepth_b, exps_a, exps_b, coords_a, coords_b,
                          ecoeffs_a, ecoeffs_b_tsp, norms_a, norms_b, boys_grid, eri2_batch);
}

void LI::kernelERI3Deriv1(const int la, const int lb, const int lc,
                          const int cdepth_a, const int cdepth_b, const int cdepth_c,
                          const double *exps_a, const double *exps_b, const double *exps_c,
                          const double *coords_a, const double *coords_b, const double *coords_c,
                          const double *ecoeffs_ab, const double *ecoeffs_deriv1_ab,
                          const double *ecoeffs_c, const double *norms_a, const double *norms_b,
                          const double *norms_c, const BoysGrid &boys_grid, double *eri3_batch)
{
    LIT::kernelERI3Deriv1(la, lb, lc, cdepth_a, cdepth_b, cdepth_c, exps_a, exps_b, exps_c,
                          coords_a, coords_b, coords_c, ecoeffs_ab, ecoeffs_deriv1_ab,
                          ecoeffs_c, norms_a, norms_b, norms_c, boys_grid, eri3_batch);
}

void LI::kernelERI4Deriv1(const int la, const int lb, const int lc, const int ld,
                          const int cdepth_a, const int cdepth_b, const int cdepth_c,
                          const int cdepth_d, const double *exps_a, const double *exps_b,
                          const double *exps_c, const double *exps_d,
                          const double *xyz_a, const double *xyz_b,
                          const double *xyz_c, const double *xyz_d,
                          const double *ecoeffs_ab, const double *ecoeffs_deriv1_ab,
                          const double *ecoeffs_cd_tsp, const double *ecoeffs_deriv1_cd_tsp,
                          const double *norms_a, const double *norms_b,
                          const double *norms_c, const double *norms_d,
                          const BoysGrid &boys_grid, double *eri4_batch)
{
    LIT::kernelERI4Deriv1(la, lb, lc, ld, cdepth_a, cdepth_b, cdepth_c, cdepth_d, exps_a, exps_b,
                          exps_c, exps_d, xyz_a, xyz_b, xyz_c, xyz_d, ecoeffs_ab, 
                          ecoeffs_deriv1_ab, ecoeffs_cd_tsp, ecoeffs_deriv1_cd_tsp,
                          norms_a, norms_b, norms_c, norms_d, boys_grid, eri4_batch);
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