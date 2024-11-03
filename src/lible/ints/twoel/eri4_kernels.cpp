#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/rints_meta.hpp>
#include <lible/ints/utils.hpp>

#include <format>

#include <fmt/core.h>

#ifdef _LIBLE_USE_MKL_
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

namespace LIT = lible::ints::two;

using LIT::kernel_eri4_t;

using std::array, std::string, std::vector;

// Table of available kernels
namespace lible::ints::two
{
    const std::array<kernel_eri4_t, 120> eri4_kernels{
        eri4Kernel<0, 0, 0, 0>,
        eri4Kernel<1, 0, 0, 0>,
        eri4Kernel<1, 0, 1, 0>,
        eri4Kernel<1, 1, 0, 0>,
        eri4Kernel<1, 1, 1, 0>,
        eri4Kernel<1, 1, 1, 1>,
        eri4Kernel<2, 0, 0, 0>,
        eri4Kernel<2, 0, 1, 0>,
        eri4Kernel<2, 0, 1, 1>,
        eri4Kernel<2, 0, 2, 0>,
        eri4Kernel<2, 1, 0, 0>,
        eri4Kernel<2, 1, 1, 0>,
        eri4Kernel<2, 1, 1, 1>,
        eri4Kernel<2, 1, 2, 0>,
        eri4Kernel<2, 1, 2, 1>,
        eri4Kernel<2, 2, 0, 0>,
        eri4Kernel<2, 2, 1, 0>,
        eri4Kernel<2, 2, 1, 1>,
        eri4Kernel<2, 2, 2, 0>,
        eri4Kernel<2, 2, 2, 1>,
        eri4Kernel<2, 2, 2, 2>,
        eri4Kernel<3, 0, 0, 0>,
        eri4Kernel<3, 0, 1, 0>,
        eri4Kernel<3, 0, 1, 1>,
        eri4Kernel<3, 0, 2, 0>,
        eri4Kernel<3, 0, 2, 1>,
        eri4Kernel<3, 0, 2, 2>,
        eri4Kernel<3, 0, 3, 0>,
        eri4Kernel<3, 1, 0, 0>,
        eri4Kernel<3, 1, 1, 0>,
        eri4Kernel<3, 1, 1, 1>,
        eri4Kernel<3, 1, 2, 0>,
        eri4Kernel<3, 1, 2, 1>,
        eri4Kernel<3, 1, 2, 2>,
        eri4Kernel<3, 1, 3, 0>,
        eri4Kernel<3, 1, 3, 1>,
        eri4Kernel<3, 2, 0, 0>,
        eri4Kernel<3, 2, 1, 0>,
        eri4Kernel<3, 2, 1, 1>,
        eri4Kernel<3, 2, 2, 0>,
        eri4Kernel<3, 2, 2, 1>,
        eri4Kernel<3, 2, 2, 2>,
        eri4Kernel<3, 2, 3, 0>,
        eri4Kernel<3, 2, 3, 1>,
        eri4Kernel<3, 2, 3, 2>,
        eri4Kernel<3, 3, 0, 0>,
        eri4Kernel<3, 3, 1, 0>,
        eri4Kernel<3, 3, 1, 1>,
        eri4Kernel<3, 3, 2, 0>,
        eri4Kernel<3, 3, 2, 1>,
        eri4Kernel<3, 3, 2, 2>,
        eri4Kernel<3, 3, 3, 0>,
        eri4Kernel<3, 3, 3, 1>,
        eri4Kernel<3, 3, 3, 2>,
        eri4Kernel<3, 3, 3, 3>,
        eri4Kernel<4, 0, 0, 0>,
        eri4Kernel<4, 0, 1, 0>,
        eri4Kernel<4, 0, 1, 1>,
        eri4Kernel<4, 0, 2, 0>,
        eri4Kernel<4, 0, 2, 1>,
        eri4Kernel<4, 0, 2, 2>,
        eri4Kernel<4, 0, 3, 0>,
        eri4Kernel<4, 0, 3, 1>,
        eri4Kernel<4, 0, 3, 2>,
        eri4Kernel<4, 0, 3, 3>,
        eri4Kernel<4, 0, 4, 0>,
        eri4Kernel<4, 1, 0, 0>,
        eri4Kernel<4, 1, 1, 0>,
        eri4Kernel<4, 1, 1, 1>,
        eri4Kernel<4, 1, 2, 0>,
        eri4Kernel<4, 1, 2, 1>,
        eri4Kernel<4, 1, 2, 2>,
        eri4Kernel<4, 1, 3, 0>,
        eri4Kernel<4, 1, 3, 1>,
        eri4Kernel<4, 1, 3, 2>,
        eri4Kernel<4, 1, 3, 3>,
        eri4Kernel<4, 1, 4, 0>,
        eri4Kernel<4, 1, 4, 1>,
        eri4Kernel<4, 2, 0, 0>,
        eri4Kernel<4, 2, 1, 0>,
        eri4Kernel<4, 2, 1, 1>,
        eri4Kernel<4, 2, 2, 0>,
        eri4Kernel<4, 2, 2, 1>,
        eri4Kernel<4, 2, 2, 2>,
        eri4Kernel<4, 2, 3, 0>,
        eri4Kernel<4, 2, 3, 1>,
        eri4Kernel<4, 2, 3, 2>,
        eri4Kernel<4, 2, 3, 3>,
        eri4Kernel<4, 2, 4, 0>,
        eri4Kernel<4, 2, 4, 1>,
        eri4Kernel<4, 2, 4, 2>,
        eri4Kernel<4, 3, 0, 0>,
        eri4Kernel<4, 3, 1, 0>,
        eri4Kernel<4, 3, 1, 1>,
        eri4Kernel<4, 3, 2, 0>,
        eri4Kernel<4, 3, 2, 1>,
        eri4Kernel<4, 3, 2, 2>,
        eri4Kernel<4, 3, 3, 0>,
        eri4Kernel<4, 3, 3, 1>,
        eri4Kernel<4, 3, 3, 2>,
        eri4Kernel<4, 3, 3, 3>,
        eri4Kernel<4, 3, 4, 0>,
        eri4Kernel<4, 3, 4, 1>,
        eri4Kernel<4, 3, 4, 2>,
        eri4Kernel<4, 3, 4, 3>,
        eri4Kernel<4, 4, 0, 0>,
        eri4Kernel<4, 4, 1, 0>,
        eri4Kernel<4, 4, 1, 1>,
        eri4Kernel<4, 4, 2, 0>,
        eri4Kernel<4, 4, 2, 1>,
        eri4Kernel<4, 4, 2, 2>,
        eri4Kernel<4, 4, 3, 0>,
        eri4Kernel<4, 4, 3, 1>,
        eri4Kernel<4, 4, 3, 2>,
        eri4Kernel<4, 4, 3, 3>,
        eri4Kernel<4, 4, 4, 0>,
        eri4Kernel<4, 4, 4, 1>,
        eri4Kernel<4, 4, 4, 2>,
        eri4Kernel<4, 4, 4, 3>,
        eri4Kernel<4, 4, 4, 4>};
}

kernel_eri4_t LIT::deployERI4Kernel(const int la, const int lb, const int lc, const int ld)
{
    int idx_ab = la * (la + 1) / 2 + lb;
    int idx_cd = lc * (lc + 1) / 2 + ld;

    if (idx_ab < idx_cd)
        throw std::runtime_error(std::format("(la, lb) >= (lc, ld) condition must be satisfied, given was: {} vs {}!\n",
                                             idx_ab, idx_cd));

    int labcd = la + lb + lc + ld;
    int l_max = _eri_kernel_max_l_;
    if (labcd > _eri_kernel_max_l_)
        throw std::runtime_error(std::format("lab + lcd = {} is larger than the allowed max: {}!\n",
                                             labcd, l_max));

    int idx_abcd = idx_ab * (idx_ab + 1) / 2 + idx_cd;

    return eri4_kernels.at(idx_abcd);
}