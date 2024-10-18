#include <lible/ints/twoel/twoel_detail.hpp>
#include <lible/ints/defs.hpp>

#include <fmt/core.h>

namespace LIT = lible::ints::two;

using LIT::eri4_kernel_t;

template <int lab, int lcd>
void eri4Kernel()
{
}

const std::map<std::pair<int, int>, eri4_kernel_t> eri4_kernels{
    {{0, 0}, eri4Kernel<0, 0>},
    {{1, 0}, eri4Kernel<1, 0>},
    {{2, 0}, eri4Kernel<2, 0>},
    {{3, 0}, eri4Kernel<3, 0>}};

eri4_kernel_t LIT::deployERI4Kernel(const int lab, const int lcd)
{
    // if (lab < lcd)
    //     throw std::runtime_error("lab >= lcd condition must be satisfied!\n");

    // if (lab + lcd > _eri_kernel_max_l_)
    //     throw std::runtime_error(fmt::format("lab + lcd = {} is larger than the allowed max: {}\n!",
                                            //  lab + lcd, _eri_kernel_max_l_));

    return eri4_kernels.at(std::make_pair(lab, lcd));
}