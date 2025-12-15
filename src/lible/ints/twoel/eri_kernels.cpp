#include <lible/ints/defs.hpp>
#include <lible/ints/ints.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <tuple>

namespace lints = lible::ints;

namespace lible::ints
{
    // regular ERI kernels

    template <int la, int lb, int lc, int ld>
    vec4d eri4KernelFun(size_t ipair_ab, size_t ipair_cd, const ShellPairData &sp_data_ab,
                        const ShellPairData &sp_data_cd, const ERI4Kernel *eri4_kernel);

    vec4d eri4KernelFun(size_t ipair_ab, size_t ipair_cd, const ShellPairData &sp_data_ab,
                        const ShellPairData &sp_data_cd, const ERI4Kernel *eri4_kernel);

    template <int la, int lb, int lc>
    vec3d eri3KernelFun(size_t ipair_ab, size_t ishell_c,
                        const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                        const ERI3Kernel *eri3_kernel);

    vec3d eri3KernelFun(size_t ipair_ab, size_t ishell_c, const ShellPairData &sp_data_ab,
                        const ShellData &sh_data_c, const ERI3Kernel *eri3_kernel);

    template <int la, int lb>
    vec2d eri2KernelFun(size_t ishell_a, size_t ishell_b, const ShellData &sh_data_a,
                        const ShellData &sh_data_b, const ERI2Kernel *eri2_kernel);

    vec2d eri2KernelFun(size_t ishell_a, size_t ishell_b, const ShellData &sh_data_a,
                        const ShellData &sh_data_b, const ERI2Kernel *eri2_kernel);

    const std::map<std::tuple<int, int, int, int>, eri4_kernelfun_t> eri4_kernelfuns{
        {{0, 0, 0, 0}, eri4KernelFun<0, 0, 0, 0>},
        {{0, 0, 1, 0}, eri4KernelFun<0, 0, 1, 0>},
        {{0, 0, 0, 1}, eri4KernelFun<0, 0, 0, 1>},
        {{0, 0, 1, 1}, eri4KernelFun<0, 0, 1, 1>},
        {{0, 0, 2, 0}, eri4KernelFun<0, 0, 2, 0>},
        {{0, 0, 0, 2}, eri4KernelFun<0, 0, 0, 2>},
        {{0, 0, 2, 1}, eri4KernelFun<0, 0, 2, 1>},
        {{0, 0, 1, 2}, eri4KernelFun<0, 0, 1, 2>},
        {{0, 0, 2, 2}, eri4KernelFun<0, 0, 2, 2>},
        {{0, 0, 3, 0}, eri4KernelFun<0, 0, 3, 0>},
        {{0, 0, 0, 3}, eri4KernelFun<0, 0, 0, 3>},
        {{0, 0, 3, 1}, eri4KernelFun<0, 0, 3, 1>},
        {{0, 0, 1, 3}, eri4KernelFun<0, 0, 1, 3>},
        {{0, 0, 3, 2}, eri4KernelFun<0, 0, 3, 2>},
        {{0, 0, 2, 3}, eri4KernelFun<0, 0, 2, 3>},
        {{0, 0, 3, 3}, eri4KernelFun<0, 0, 3, 3>},
        {{0, 0, 4, 0}, eri4KernelFun<0, 0, 4, 0>},
        {{0, 0, 0, 4}, eri4KernelFun<0, 0, 0, 4>},
        {{0, 0, 4, 1}, eri4KernelFun<0, 0, 4, 1>},
        {{0, 0, 1, 4}, eri4KernelFun<0, 0, 1, 4>},
        {{0, 0, 4, 2}, eri4KernelFun<0, 0, 4, 2>},
        {{0, 0, 2, 4}, eri4KernelFun<0, 0, 2, 4>},
        {{0, 0, 5, 0}, eri4KernelFun<0, 0, 5, 0>},
        {{0, 0, 0, 5}, eri4KernelFun<0, 0, 0, 5>},
        {{0, 0, 5, 1}, eri4KernelFun<0, 0, 5, 1>},
        {{0, 0, 1, 5}, eri4KernelFun<0, 0, 1, 5>},
        {{0, 0, 6, 0}, eri4KernelFun<0, 0, 6, 0>},
        {{0, 0, 0, 6}, eri4KernelFun<0, 0, 0, 6>},
        {{1, 0, 0, 0}, eri4KernelFun<1, 0, 0, 0>},
        {{0, 1, 0, 0}, eri4KernelFun<0, 1, 0, 0>},
        {{1, 0, 1, 0}, eri4KernelFun<1, 0, 1, 0>},
        {{0, 1, 0, 1}, eri4KernelFun<0, 1, 0, 1>},
        {{1, 0, 0, 1}, eri4KernelFun<1, 0, 0, 1>},
        {{0, 1, 1, 0}, eri4KernelFun<0, 1, 1, 0>},
        {{1, 0, 1, 1}, eri4KernelFun<1, 0, 1, 1>},
        {{0, 1, 1, 1}, eri4KernelFun<0, 1, 1, 1>},
        {{1, 0, 2, 0}, eri4KernelFun<1, 0, 2, 0>},
        {{0, 1, 0, 2}, eri4KernelFun<0, 1, 0, 2>},
        {{1, 0, 0, 2}, eri4KernelFun<1, 0, 0, 2>},
        {{0, 1, 2, 0}, eri4KernelFun<0, 1, 2, 0>},
        {{1, 0, 2, 1}, eri4KernelFun<1, 0, 2, 1>},
        {{0, 1, 1, 2}, eri4KernelFun<0, 1, 1, 2>},
        {{1, 0, 1, 2}, eri4KernelFun<1, 0, 1, 2>},
        {{0, 1, 2, 1}, eri4KernelFun<0, 1, 2, 1>},
        {{1, 0, 2, 2}, eri4KernelFun<1, 0, 2, 2>},
        {{0, 1, 2, 2}, eri4KernelFun<0, 1, 2, 2>},
        {{1, 0, 3, 0}, eri4KernelFun<1, 0, 3, 0>},
        {{0, 1, 0, 3}, eri4KernelFun<0, 1, 0, 3>},
        {{1, 0, 0, 3}, eri4KernelFun<1, 0, 0, 3>},
        {{0, 1, 3, 0}, eri4KernelFun<0, 1, 3, 0>},
        {{1, 0, 3, 1}, eri4KernelFun<1, 0, 3, 1>},
        {{0, 1, 1, 3}, eri4KernelFun<0, 1, 1, 3>},
        {{1, 0, 1, 3}, eri4KernelFun<1, 0, 1, 3>},
        {{0, 1, 3, 1}, eri4KernelFun<0, 1, 3, 1>},
        {{1, 0, 3, 2}, eri4KernelFun<1, 0, 3, 2>},
        {{0, 1, 2, 3}, eri4KernelFun<0, 1, 2, 3>},
        {{1, 0, 2, 3}, eri4KernelFun<1, 0, 2, 3>},
        {{0, 1, 3, 2}, eri4KernelFun<0, 1, 3, 2>},
        {{1, 0, 4, 0}, eri4KernelFun<1, 0, 4, 0>},
        {{0, 1, 0, 4}, eri4KernelFun<0, 1, 0, 4>},
        {{1, 0, 0, 4}, eri4KernelFun<1, 0, 0, 4>},
        {{0, 1, 4, 0}, eri4KernelFun<0, 1, 4, 0>},
        {{1, 0, 4, 1}, eri4KernelFun<1, 0, 4, 1>},
        {{0, 1, 1, 4}, eri4KernelFun<0, 1, 1, 4>},
        {{1, 0, 1, 4}, eri4KernelFun<1, 0, 1, 4>},
        {{0, 1, 4, 1}, eri4KernelFun<0, 1, 4, 1>},
        {{1, 0, 5, 0}, eri4KernelFun<1, 0, 5, 0>},
        {{0, 1, 0, 5}, eri4KernelFun<0, 1, 0, 5>},
        {{1, 0, 0, 5}, eri4KernelFun<1, 0, 0, 5>},
        {{0, 1, 5, 0}, eri4KernelFun<0, 1, 5, 0>},
        {{1, 1, 0, 0}, eri4KernelFun<1, 1, 0, 0>},
        {{1, 1, 1, 0}, eri4KernelFun<1, 1, 1, 0>},
        {{1, 1, 0, 1}, eri4KernelFun<1, 1, 0, 1>},
        {{1, 1, 1, 1}, eri4KernelFun<1, 1, 1, 1>},
        {{1, 1, 2, 0}, eri4KernelFun<1, 1, 2, 0>},
        {{1, 1, 0, 2}, eri4KernelFun<1, 1, 0, 2>},
        {{1, 1, 2, 1}, eri4KernelFun<1, 1, 2, 1>},
        {{1, 1, 1, 2}, eri4KernelFun<1, 1, 1, 2>},
        {{1, 1, 2, 2}, eri4KernelFun<1, 1, 2, 2>},
        {{1, 1, 3, 0}, eri4KernelFun<1, 1, 3, 0>},
        {{1, 1, 0, 3}, eri4KernelFun<1, 1, 0, 3>},
        {{1, 1, 3, 1}, eri4KernelFun<1, 1, 3, 1>},
        {{1, 1, 1, 3}, eri4KernelFun<1, 1, 1, 3>},
        {{1, 1, 4, 0}, eri4KernelFun<1, 1, 4, 0>},
        {{1, 1, 0, 4}, eri4KernelFun<1, 1, 0, 4>},
        {{2, 0, 0, 0}, eri4KernelFun<2, 0, 0, 0>},
        {{0, 2, 0, 0}, eri4KernelFun<0, 2, 0, 0>},
        {{2, 0, 1, 0}, eri4KernelFun<2, 0, 1, 0>},
        {{0, 2, 0, 1}, eri4KernelFun<0, 2, 0, 1>},
        {{2, 0, 0, 1}, eri4KernelFun<2, 0, 0, 1>},
        {{0, 2, 1, 0}, eri4KernelFun<0, 2, 1, 0>},
        {{2, 0, 1, 1}, eri4KernelFun<2, 0, 1, 1>},
        {{0, 2, 1, 1}, eri4KernelFun<0, 2, 1, 1>},
        {{2, 0, 2, 0}, eri4KernelFun<2, 0, 2, 0>},
        {{0, 2, 0, 2}, eri4KernelFun<0, 2, 0, 2>},
        {{2, 0, 0, 2}, eri4KernelFun<2, 0, 0, 2>},
        {{0, 2, 2, 0}, eri4KernelFun<0, 2, 2, 0>},
        {{2, 0, 2, 1}, eri4KernelFun<2, 0, 2, 1>},
        {{0, 2, 1, 2}, eri4KernelFun<0, 2, 1, 2>},
        {{2, 0, 1, 2}, eri4KernelFun<2, 0, 1, 2>},
        {{0, 2, 2, 1}, eri4KernelFun<0, 2, 2, 1>},
        {{2, 0, 2, 2}, eri4KernelFun<2, 0, 2, 2>},
        {{0, 2, 2, 2}, eri4KernelFun<0, 2, 2, 2>},
        {{2, 0, 3, 0}, eri4KernelFun<2, 0, 3, 0>},
        {{0, 2, 0, 3}, eri4KernelFun<0, 2, 0, 3>},
        {{2, 0, 0, 3}, eri4KernelFun<2, 0, 0, 3>},
        {{0, 2, 3, 0}, eri4KernelFun<0, 2, 3, 0>},
        {{2, 0, 3, 1}, eri4KernelFun<2, 0, 3, 1>},
        {{0, 2, 1, 3}, eri4KernelFun<0, 2, 1, 3>},
        {{2, 0, 1, 3}, eri4KernelFun<2, 0, 1, 3>},
        {{0, 2, 3, 1}, eri4KernelFun<0, 2, 3, 1>},
        {{2, 0, 4, 0}, eri4KernelFun<2, 0, 4, 0>},
        {{0, 2, 0, 4}, eri4KernelFun<0, 2, 0, 4>},
        {{2, 0, 0, 4}, eri4KernelFun<2, 0, 0, 4>},
        {{0, 2, 4, 0}, eri4KernelFun<0, 2, 4, 0>},
        {{2, 1, 0, 0}, eri4KernelFun<2, 1, 0, 0>},
        {{1, 2, 0, 0}, eri4KernelFun<1, 2, 0, 0>},
        {{2, 1, 1, 0}, eri4KernelFun<2, 1, 1, 0>},
        {{1, 2, 0, 1}, eri4KernelFun<1, 2, 0, 1>},
        {{2, 1, 0, 1}, eri4KernelFun<2, 1, 0, 1>},
        {{1, 2, 1, 0}, eri4KernelFun<1, 2, 1, 0>},
        {{2, 1, 1, 1}, eri4KernelFun<2, 1, 1, 1>},
        {{1, 2, 1, 1}, eri4KernelFun<1, 2, 1, 1>},
        {{2, 1, 2, 0}, eri4KernelFun<2, 1, 2, 0>},
        {{1, 2, 0, 2}, eri4KernelFun<1, 2, 0, 2>},
        {{2, 1, 0, 2}, eri4KernelFun<2, 1, 0, 2>},
        {{1, 2, 2, 0}, eri4KernelFun<1, 2, 2, 0>},
        {{2, 1, 2, 1}, eri4KernelFun<2, 1, 2, 1>},
        {{1, 2, 1, 2}, eri4KernelFun<1, 2, 1, 2>},
        {{2, 1, 1, 2}, eri4KernelFun<2, 1, 1, 2>},
        {{1, 2, 2, 1}, eri4KernelFun<1, 2, 2, 1>},
        {{2, 1, 3, 0}, eri4KernelFun<2, 1, 3, 0>},
        {{1, 2, 0, 3}, eri4KernelFun<1, 2, 0, 3>},
        {{2, 1, 0, 3}, eri4KernelFun<2, 1, 0, 3>},
        {{1, 2, 3, 0}, eri4KernelFun<1, 2, 3, 0>},
        {{2, 2, 0, 0}, eri4KernelFun<2, 2, 0, 0>},
        {{2, 2, 1, 0}, eri4KernelFun<2, 2, 1, 0>},
        {{2, 2, 0, 1}, eri4KernelFun<2, 2, 0, 1>},
        {{2, 2, 1, 1}, eri4KernelFun<2, 2, 1, 1>},
        {{2, 2, 2, 0}, eri4KernelFun<2, 2, 2, 0>},
        {{2, 2, 0, 2}, eri4KernelFun<2, 2, 0, 2>},
        {{3, 0, 0, 0}, eri4KernelFun<3, 0, 0, 0>},
        {{0, 3, 0, 0}, eri4KernelFun<0, 3, 0, 0>},
        {{3, 0, 1, 0}, eri4KernelFun<3, 0, 1, 0>},
        {{0, 3, 0, 1}, eri4KernelFun<0, 3, 0, 1>},
        {{3, 0, 0, 1}, eri4KernelFun<3, 0, 0, 1>},
        {{0, 3, 1, 0}, eri4KernelFun<0, 3, 1, 0>},
        {{3, 0, 1, 1}, eri4KernelFun<3, 0, 1, 1>},
        {{0, 3, 1, 1}, eri4KernelFun<0, 3, 1, 1>},
        {{3, 0, 2, 0}, eri4KernelFun<3, 0, 2, 0>},
        {{0, 3, 0, 2}, eri4KernelFun<0, 3, 0, 2>},
        {{3, 0, 0, 2}, eri4KernelFun<3, 0, 0, 2>},
        {{0, 3, 2, 0}, eri4KernelFun<0, 3, 2, 0>},
        {{3, 0, 2, 1}, eri4KernelFun<3, 0, 2, 1>},
        {{0, 3, 1, 2}, eri4KernelFun<0, 3, 1, 2>},
        {{3, 0, 1, 2}, eri4KernelFun<3, 0, 1, 2>},
        {{0, 3, 2, 1}, eri4KernelFun<0, 3, 2, 1>},
        {{3, 0, 3, 0}, eri4KernelFun<3, 0, 3, 0>},
        {{0, 3, 0, 3}, eri4KernelFun<0, 3, 0, 3>},
        {{3, 0, 0, 3}, eri4KernelFun<3, 0, 0, 3>},
        {{0, 3, 3, 0}, eri4KernelFun<0, 3, 3, 0>},
        {{3, 1, 0, 0}, eri4KernelFun<3, 1, 0, 0>},
        {{1, 3, 0, 0}, eri4KernelFun<1, 3, 0, 0>},
        {{3, 1, 1, 0}, eri4KernelFun<3, 1, 1, 0>},
        {{1, 3, 0, 1}, eri4KernelFun<1, 3, 0, 1>},
        {{3, 1, 0, 1}, eri4KernelFun<3, 1, 0, 1>},
        {{1, 3, 1, 0}, eri4KernelFun<1, 3, 1, 0>},
        {{3, 1, 1, 1}, eri4KernelFun<3, 1, 1, 1>},
        {{1, 3, 1, 1}, eri4KernelFun<1, 3, 1, 1>},
        {{3, 1, 2, 0}, eri4KernelFun<3, 1, 2, 0>},
        {{1, 3, 0, 2}, eri4KernelFun<1, 3, 0, 2>},
        {{3, 1, 0, 2}, eri4KernelFun<3, 1, 0, 2>},
        {{1, 3, 2, 0}, eri4KernelFun<1, 3, 2, 0>},
        {{3, 2, 0, 0}, eri4KernelFun<3, 2, 0, 0>},
        {{2, 3, 0, 0}, eri4KernelFun<2, 3, 0, 0>},
        {{3, 2, 1, 0}, eri4KernelFun<3, 2, 1, 0>},
        {{2, 3, 0, 1}, eri4KernelFun<2, 3, 0, 1>},
        {{3, 2, 0, 1}, eri4KernelFun<3, 2, 0, 1>},
        {{2, 3, 1, 0}, eri4KernelFun<2, 3, 1, 0>},
        {{3, 3, 0, 0}, eri4KernelFun<3, 3, 0, 0>},
        {{4, 0, 0, 0}, eri4KernelFun<4, 0, 0, 0>},
        {{0, 4, 0, 0}, eri4KernelFun<0, 4, 0, 0>},
        {{4, 0, 1, 0}, eri4KernelFun<4, 0, 1, 0>},
        {{0, 4, 0, 1}, eri4KernelFun<0, 4, 0, 1>},
        {{4, 0, 0, 1}, eri4KernelFun<4, 0, 0, 1>},
        {{0, 4, 1, 0}, eri4KernelFun<0, 4, 1, 0>},
        {{4, 0, 1, 1}, eri4KernelFun<4, 0, 1, 1>},
        {{0, 4, 1, 1}, eri4KernelFun<0, 4, 1, 1>},
        {{4, 0, 2, 0}, eri4KernelFun<4, 0, 2, 0>},
        {{0, 4, 0, 2}, eri4KernelFun<0, 4, 0, 2>},
        {{4, 0, 0, 2}, eri4KernelFun<4, 0, 0, 2>},
        {{0, 4, 2, 0}, eri4KernelFun<0, 4, 2, 0>},
        {{4, 1, 0, 0}, eri4KernelFun<4, 1, 0, 0>},
        {{1, 4, 0, 0}, eri4KernelFun<1, 4, 0, 0>},
        {{4, 1, 1, 0}, eri4KernelFun<4, 1, 1, 0>},
        {{1, 4, 0, 1}, eri4KernelFun<1, 4, 0, 1>},
        {{4, 1, 0, 1}, eri4KernelFun<4, 1, 0, 1>},
        {{1, 4, 1, 0}, eri4KernelFun<1, 4, 1, 0>},
        {{4, 2, 0, 0}, eri4KernelFun<4, 2, 0, 0>},
        {{2, 4, 0, 0}, eri4KernelFun<2, 4, 0, 0>},
        {{5, 0, 0, 0}, eri4KernelFun<5, 0, 0, 0>},
        {{0, 5, 0, 0}, eri4KernelFun<0, 5, 0, 0>},
        {{5, 0, 1, 0}, eri4KernelFun<5, 0, 1, 0>},
        {{0, 5, 0, 1}, eri4KernelFun<0, 5, 0, 1>},
        {{5, 0, 0, 1}, eri4KernelFun<5, 0, 0, 1>},
        {{0, 5, 1, 0}, eri4KernelFun<0, 5, 1, 0>},
        {{5, 1, 0, 0}, eri4KernelFun<5, 1, 0, 0>},
        {{1, 5, 0, 0}, eri4KernelFun<1, 5, 0, 0>},
        {{6, 0, 0, 0}, eri4KernelFun<6, 0, 0, 0>},
        {{0, 6, 0, 0}, eri4KernelFun<0, 6, 0, 0>}
    };

    const std::map<std::tuple<int, int, int>, eri3_kernelfun_t> eri3_kernelfuns{
        {{0, 0, 0}, eri3KernelFun<0, 0, 0>},
        {{0, 0, 1}, eri3KernelFun<0, 0, 1>},
        {{0, 0, 2}, eri3KernelFun<0, 0, 2>},
        {{0, 0, 3}, eri3KernelFun<0, 0, 3>},
        {{0, 0, 4}, eri3KernelFun<0, 0, 4>},
        {{0, 0, 5}, eri3KernelFun<0, 0, 5>},
        {{0, 0, 6}, eri3KernelFun<0, 0, 6>},
        {{1, 0, 0}, eri3KernelFun<1, 0, 0>},
        {{0, 1, 0}, eri3KernelFun<0, 1, 0>},
        {{1, 0, 1}, eri3KernelFun<1, 0, 1>},
        {{0, 1, 1}, eri3KernelFun<0, 1, 1>},
        {{1, 0, 2}, eri3KernelFun<1, 0, 2>},
        {{0, 1, 2}, eri3KernelFun<0, 1, 2>},
        {{1, 0, 3}, eri3KernelFun<1, 0, 3>},
        {{0, 1, 3}, eri3KernelFun<0, 1, 3>},
        {{1, 0, 4}, eri3KernelFun<1, 0, 4>},
        {{0, 1, 4}, eri3KernelFun<0, 1, 4>},
        {{1, 0, 5}, eri3KernelFun<1, 0, 5>},
        {{0, 1, 5}, eri3KernelFun<0, 1, 5>},
        {{1, 1, 0}, eri3KernelFun<1, 1, 0>},
        {{1, 1, 1}, eri3KernelFun<1, 1, 1>},
        {{1, 1, 2}, eri3KernelFun<1, 1, 2>},
        {{1, 1, 3}, eri3KernelFun<1, 1, 3>},
        {{1, 1, 4}, eri3KernelFun<1, 1, 4>},
        {{2, 0, 0}, eri3KernelFun<2, 0, 0>},
        {{0, 2, 0}, eri3KernelFun<0, 2, 0>},
        {{2, 0, 1}, eri3KernelFun<2, 0, 1>},
        {{0, 2, 1}, eri3KernelFun<0, 2, 1>},
        {{2, 0, 2}, eri3KernelFun<2, 0, 2>},
        {{0, 2, 2}, eri3KernelFun<0, 2, 2>},
        {{2, 0, 3}, eri3KernelFun<2, 0, 3>},
        {{0, 2, 3}, eri3KernelFun<0, 2, 3>},
        {{2, 0, 4}, eri3KernelFun<2, 0, 4>},
        {{0, 2, 4}, eri3KernelFun<0, 2, 4>},
        {{2, 1, 0}, eri3KernelFun<2, 1, 0>},
        {{1, 2, 0}, eri3KernelFun<1, 2, 0>},
        {{2, 1, 1}, eri3KernelFun<2, 1, 1>},
        {{1, 2, 1}, eri3KernelFun<1, 2, 1>},
        {{2, 1, 2}, eri3KernelFun<2, 1, 2>},
        {{1, 2, 2}, eri3KernelFun<1, 2, 2>},
        {{2, 1, 3}, eri3KernelFun<2, 1, 3>},
        {{1, 2, 3}, eri3KernelFun<1, 2, 3>},
        {{2, 2, 0}, eri3KernelFun<2, 2, 0>},
        {{2, 2, 1}, eri3KernelFun<2, 2, 1>},
        {{2, 2, 2}, eri3KernelFun<2, 2, 2>},
        {{3, 0, 0}, eri3KernelFun<3, 0, 0>},
        {{0, 3, 0}, eri3KernelFun<0, 3, 0>},
        {{3, 0, 1}, eri3KernelFun<3, 0, 1>},
        {{0, 3, 1}, eri3KernelFun<0, 3, 1>},
        {{3, 0, 2}, eri3KernelFun<3, 0, 2>},
        {{0, 3, 2}, eri3KernelFun<0, 3, 2>},
        {{3, 0, 3}, eri3KernelFun<3, 0, 3>},
        {{0, 3, 3}, eri3KernelFun<0, 3, 3>},
        {{3, 1, 0}, eri3KernelFun<3, 1, 0>},
        {{1, 3, 0}, eri3KernelFun<1, 3, 0>},
        {{3, 1, 1}, eri3KernelFun<3, 1, 1>},
        {{1, 3, 1}, eri3KernelFun<1, 3, 1>},
        {{3, 1, 2}, eri3KernelFun<3, 1, 2>},
        {{1, 3, 2}, eri3KernelFun<1, 3, 2>},
        {{3, 2, 0}, eri3KernelFun<3, 2, 0>},
        {{2, 3, 0}, eri3KernelFun<2, 3, 0>},
        {{3, 2, 1}, eri3KernelFun<3, 2, 1>},
        {{2, 3, 1}, eri3KernelFun<2, 3, 1>},
        {{3, 3, 0}, eri3KernelFun<3, 3, 0>},
        {{4, 0, 0}, eri3KernelFun<4, 0, 0>},
        {{0, 4, 0}, eri3KernelFun<0, 4, 0>},
        {{4, 0, 1}, eri3KernelFun<4, 0, 1>},
        {{0, 4, 1}, eri3KernelFun<0, 4, 1>},
        {{4, 0, 2}, eri3KernelFun<4, 0, 2>},
        {{0, 4, 2}, eri3KernelFun<0, 4, 2>},
        {{4, 1, 0}, eri3KernelFun<4, 1, 0>},
        {{1, 4, 0}, eri3KernelFun<1, 4, 0>},
        {{4, 1, 1}, eri3KernelFun<4, 1, 1>},
        {{1, 4, 1}, eri3KernelFun<1, 4, 1>},
        {{4, 2, 0}, eri3KernelFun<4, 2, 0>},
        {{2, 4, 0}, eri3KernelFun<2, 4, 0>},
        {{5, 0, 0}, eri3KernelFun<5, 0, 0>},
        {{0, 5, 0}, eri3KernelFun<0, 5, 0>},
        {{5, 0, 1}, eri3KernelFun<5, 0, 1>},
        {{0, 5, 1}, eri3KernelFun<0, 5, 1>},
        {{5, 1, 0}, eri3KernelFun<5, 1, 0>},
        {{1, 5, 0}, eri3KernelFun<1, 5, 0>},
        {{6, 0, 0}, eri3KernelFun<6, 0, 0>},
        {{0, 6, 0}, eri3KernelFun<0, 6, 0>}
    };

    const std::map<std::tuple<int, int>, eri2_kernelfun_t> eri2_kernelfuns{
        {{0, 0}, eri2KernelFun<0, 0>},
        {{0, 1}, eri2KernelFun<0, 1>},
        {{0, 2}, eri2KernelFun<0, 2>},
        {{0, 3}, eri2KernelFun<0, 3>},
        {{0, 4}, eri2KernelFun<0, 4>},
        {{0, 5}, eri2KernelFun<0, 5>},
        {{0, 6}, eri2KernelFun<0, 6>},
        {{1, 0}, eri2KernelFun<1, 0>},
        {{1, 1}, eri2KernelFun<1, 1>},
        {{1, 2}, eri2KernelFun<1, 2>},
        {{1, 3}, eri2KernelFun<1, 3>},
        {{1, 4}, eri2KernelFun<1, 4>},
        {{1, 5}, eri2KernelFun<1, 5>},
        {{2, 0}, eri2KernelFun<2, 0>},
        {{2, 1}, eri2KernelFun<2, 1>},
        {{2, 2}, eri2KernelFun<2, 2>},
        {{2, 3}, eri2KernelFun<2, 3>},
        {{2, 4}, eri2KernelFun<2, 4>},
        {{3, 0}, eri2KernelFun<3, 0>},
        {{3, 1}, eri2KernelFun<3, 1>},
        {{3, 2}, eri2KernelFun<3, 2>},
        {{3, 3}, eri2KernelFun<3, 3>},
        {{4, 0}, eri2KernelFun<4, 0>},
        {{4, 1}, eri2KernelFun<4, 1>},
        {{4, 2}, eri2KernelFun<4, 2>},
        {{5, 0}, eri2KernelFun<5, 0>},
        {{5, 1}, eri2KernelFun<5, 1>},
        {{6, 0}, eri2KernelFun<6, 0>}
    };
}

namespace lible::ints
{
    // Derivative ERI kernels

    // 2-center
    template <int la, int lb>
    std::array<vec2d, 6> eri2d1KernelFun(size_t ishell_a, size_t ishell_b,
                                         const ShellData &sh_data_a, const ShellData &sh_data_b,
                                         const ERI2D1Kernel *eri2d1_kernel);

    std::array<vec2d, 6> eri2d1KernelFun(size_t ishell_a, size_t ishell_b,
                                         const ShellData &sh_data_a, const ShellData &sh_data_b,
                                         const ERI2D1Kernel *eri2d1_kernel);

    template <int la, int lb>
    arr2d<vec2d, 6, 6> eri2d2KernelFun(size_t ishell_a, size_t ishell_b,
                                       const ShellData &sh_data_a,
                                       const ShellData &sh_data_b,
                                       const ERI2D2Kernel *eri2d2_kernel);

    arr2d<vec2d, 6, 6> eri2d2KernelFun(size_t ishell_a, size_t ishell_b,
                                       const ShellData &sh_data_a,
                                       const ShellData &sh_data_b,
                                       const ERI2D2Kernel *eri2d2_kernel);

    // 3-center
    template <int la, int lb, int lc>
    std::array<vec3d, 9> eri3d1KernelFun(size_t ipair_ab, size_t ishell_c,
                                         const ShellPairData &sp_data_ab,
                                         const ShellData &sh_data_c,
                                         const ERI3D1Kernel *eri3d1_kernel);

    std::array<vec3d, 9> eri3d1KernelFun(size_t ipair_ab, size_t ishell_c,
                                         const ShellPairData &sp_data_ab,
                                         const ShellData &sh_data_c,
                                         const ERI3D1Kernel *eri3d1_kernel);

    template <int la, int lb, int lc>
    std::array<vec3d, 3> eri3socKernelFun(size_t ipair_ab, size_t ishell_c,
                                          const ShellPairData &sp_data_ab,
                                          const ShellData &sh_data_c,
                                          const ERI3SOCKernel *eri3soc_kernel);

    std::array<vec3d, 3> eri3socKernelFun(size_t ipair_ab, size_t ishell_c,
                                          const ShellPairData &sp_data_ab,
                                          const ShellData &sh_data_c,
                                          const ERI3SOCKernel *eri3soc_kernel);

    // 4-center
    template <int la, int lb, int lc, int ld>
    std::array<vec4d, 12> eri4d1KernelFun(size_t ipair_ab, size_t ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4D1Kernel *eri4d1_kernel);

    std::array<vec4d, 12> eri4d1KernelFun(size_t ipair_ab, size_t ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4D1Kernel *eri4d1_kernel);

    template <int la, int lb, int lc, int ld>
    std::array<vec4d, 3> eri4socKernelFun(size_t ipair_ab, size_t ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4SOCKernel *eri4soc_kernel);

    std::array<vec4d, 3> eri4socKernelFun(size_t ipair_ab, size_t ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4SOCKernel *eri4soc_kernel);

    const std::map<std::tuple<int, int>, eri2d1_kernelfun_t> eri2d1_kernelfuns{
        {{0, 0}, eri2d1KernelFun<0, 0>},
        {{0, 1}, eri2d1KernelFun<0, 1>},
        {{0, 2}, eri2d1KernelFun<0, 2>},
        {{0, 3}, eri2d1KernelFun<0, 3>},
        {{0, 4}, eri2d1KernelFun<0, 4>},
        {{0, 5}, eri2d1KernelFun<0, 5>},
        {{0, 6}, eri2d1KernelFun<0, 6>},
        {{1, 0}, eri2d1KernelFun<1, 0>},
        {{1, 1}, eri2d1KernelFun<1, 1>},
        {{1, 2}, eri2d1KernelFun<1, 2>},
        {{1, 3}, eri2d1KernelFun<1, 3>},
        {{1, 4}, eri2d1KernelFun<1, 4>},
        {{1, 5}, eri2d1KernelFun<1, 5>},
        {{2, 0}, eri2d1KernelFun<2, 0>},
        {{2, 1}, eri2d1KernelFun<2, 1>},
        {{2, 2}, eri2d1KernelFun<2, 2>},
        {{2, 3}, eri2d1KernelFun<2, 3>},
        {{2, 4}, eri2d1KernelFun<2, 4>},
        {{3, 0}, eri2d1KernelFun<3, 0>},
        {{3, 1}, eri2d1KernelFun<3, 1>},
        {{3, 2}, eri2d1KernelFun<3, 2>},
        {{3, 3}, eri2d1KernelFun<3, 3>},
        {{4, 0}, eri2d1KernelFun<4, 0>},
        {{4, 1}, eri2d1KernelFun<4, 1>},
        {{4, 2}, eri2d1KernelFun<4, 2>},
        {{5, 0}, eri2d1KernelFun<5, 0>},
        {{5, 1}, eri2d1KernelFun<5, 1>},
        {{6, 0}, eri2d1KernelFun<6, 0>}
    };

    const std::map<std::tuple<int, int>, eri2d2_kernelfun_t> eri2d2_kernelfuns{
        {{0, 0}, eri2d2KernelFun<0, 0>},
        {{0, 1}, eri2d2KernelFun<0, 1>},
        {{0, 2}, eri2d2KernelFun<0, 2>},
        {{0, 3}, eri2d2KernelFun<0, 3>},
        {{0, 4}, eri2d2KernelFun<0, 4>},
        {{0, 5}, eri2d2KernelFun<0, 5>},
        {{0, 6}, eri2d2KernelFun<0, 6>},
        {{1, 0}, eri2d2KernelFun<1, 0>},
        {{1, 1}, eri2d2KernelFun<1, 1>},
        {{1, 2}, eri2d2KernelFun<1, 2>},
        {{1, 3}, eri2d2KernelFun<1, 3>},
        {{1, 4}, eri2d2KernelFun<1, 4>},
        {{1, 5}, eri2d2KernelFun<1, 5>},
        {{2, 0}, eri2d2KernelFun<2, 0>},
        {{2, 1}, eri2d2KernelFun<2, 1>},
        {{2, 2}, eri2d2KernelFun<2, 2>},
        {{2, 3}, eri2d2KernelFun<2, 3>},
        {{2, 4}, eri2d2KernelFun<2, 4>},
        {{3, 0}, eri2d2KernelFun<3, 0>},
        {{3, 1}, eri2d2KernelFun<3, 1>},
        {{3, 2}, eri2d2KernelFun<3, 2>},
        {{3, 3}, eri2d2KernelFun<3, 3>},
        {{4, 0}, eri2d2KernelFun<4, 0>},
        {{4, 1}, eri2d2KernelFun<4, 1>},
        {{4, 2}, eri2d2KernelFun<4, 2>},
        {{5, 0}, eri2d2KernelFun<5, 0>},
        {{5, 1}, eri2d2KernelFun<5, 1>},
        {{6, 0}, eri2d2KernelFun<6, 0>}
    };

    const std::map<std::tuple<int, int, int>, eri3d1_kernelfun_t> eri3d1_kernelfuns{
        {{0, 0, 0}, eri3d1KernelFun<0, 0, 0>},
        {{0, 0, 1}, eri3d1KernelFun<0, 0, 1>},
        {{0, 0, 2}, eri3d1KernelFun<0, 0, 2>},
        {{0, 0, 3}, eri3d1KernelFun<0, 0, 3>},
        {{0, 0, 4}, eri3d1KernelFun<0, 0, 4>},
        {{0, 0, 5}, eri3d1KernelFun<0, 0, 5>},
        {{0, 0, 6}, eri3d1KernelFun<0, 0, 6>},
        {{1, 0, 0}, eri3d1KernelFun<1, 0, 0>},
        {{0, 1, 0}, eri3d1KernelFun<0, 1, 0>},
        {{1, 0, 1}, eri3d1KernelFun<1, 0, 1>},
        {{0, 1, 1}, eri3d1KernelFun<0, 1, 1>},
        {{1, 0, 2}, eri3d1KernelFun<1, 0, 2>},
        {{0, 1, 2}, eri3d1KernelFun<0, 1, 2>},
        {{1, 0, 3}, eri3d1KernelFun<1, 0, 3>},
        {{0, 1, 3}, eri3d1KernelFun<0, 1, 3>},
        {{1, 0, 4}, eri3d1KernelFun<1, 0, 4>},
        {{0, 1, 4}, eri3d1KernelFun<0, 1, 4>},
        {{1, 0, 5}, eri3d1KernelFun<1, 0, 5>},
        {{0, 1, 5}, eri3d1KernelFun<0, 1, 5>},
        {{1, 1, 0}, eri3d1KernelFun<1, 1, 0>},
        {{1, 1, 1}, eri3d1KernelFun<1, 1, 1>},
        {{1, 1, 2}, eri3d1KernelFun<1, 1, 2>},
        {{1, 1, 3}, eri3d1KernelFun<1, 1, 3>},
        {{1, 1, 4}, eri3d1KernelFun<1, 1, 4>},
        {{2, 0, 0}, eri3d1KernelFun<2, 0, 0>},
        {{0, 2, 0}, eri3d1KernelFun<0, 2, 0>},
        {{2, 0, 1}, eri3d1KernelFun<2, 0, 1>},
        {{0, 2, 1}, eri3d1KernelFun<0, 2, 1>},
        {{2, 0, 2}, eri3d1KernelFun<2, 0, 2>},
        {{0, 2, 2}, eri3d1KernelFun<0, 2, 2>},
        {{2, 0, 3}, eri3d1KernelFun<2, 0, 3>},
        {{0, 2, 3}, eri3d1KernelFun<0, 2, 3>},
        {{2, 0, 4}, eri3d1KernelFun<2, 0, 4>},
        {{0, 2, 4}, eri3d1KernelFun<0, 2, 4>},
        {{2, 1, 0}, eri3d1KernelFun<2, 1, 0>},
        {{1, 2, 0}, eri3d1KernelFun<1, 2, 0>},
        {{2, 1, 1}, eri3d1KernelFun<2, 1, 1>},
        {{1, 2, 1}, eri3d1KernelFun<1, 2, 1>},
        {{2, 1, 2}, eri3d1KernelFun<2, 1, 2>},
        {{1, 2, 2}, eri3d1KernelFun<1, 2, 2>},
        {{2, 1, 3}, eri3d1KernelFun<2, 1, 3>},
        {{1, 2, 3}, eri3d1KernelFun<1, 2, 3>},
        {{2, 2, 0}, eri3d1KernelFun<2, 2, 0>},
        {{2, 2, 1}, eri3d1KernelFun<2, 2, 1>},
        {{2, 2, 2}, eri3d1KernelFun<2, 2, 2>},
        {{3, 0, 0}, eri3d1KernelFun<3, 0, 0>},
        {{0, 3, 0}, eri3d1KernelFun<0, 3, 0>},
        {{3, 0, 1}, eri3d1KernelFun<3, 0, 1>},
        {{0, 3, 1}, eri3d1KernelFun<0, 3, 1>},
        {{3, 0, 2}, eri3d1KernelFun<3, 0, 2>},
        {{0, 3, 2}, eri3d1KernelFun<0, 3, 2>},
        {{3, 0, 3}, eri3d1KernelFun<3, 0, 3>},
        {{0, 3, 3}, eri3d1KernelFun<0, 3, 3>},
        {{3, 1, 0}, eri3d1KernelFun<3, 1, 0>},
        {{1, 3, 0}, eri3d1KernelFun<1, 3, 0>},
        {{3, 1, 1}, eri3d1KernelFun<3, 1, 1>},
        {{1, 3, 1}, eri3d1KernelFun<1, 3, 1>},
        {{3, 1, 2}, eri3d1KernelFun<3, 1, 2>},
        {{1, 3, 2}, eri3d1KernelFun<1, 3, 2>},
        {{3, 2, 0}, eri3d1KernelFun<3, 2, 0>},
        {{2, 3, 0}, eri3d1KernelFun<2, 3, 0>},
        {{3, 2, 1}, eri3d1KernelFun<3, 2, 1>},
        {{2, 3, 1}, eri3d1KernelFun<2, 3, 1>},
        {{3, 3, 0}, eri3d1KernelFun<3, 3, 0>},
        {{4, 0, 0}, eri3d1KernelFun<4, 0, 0>},
        {{0, 4, 0}, eri3d1KernelFun<0, 4, 0>},
        {{4, 0, 1}, eri3d1KernelFun<4, 0, 1>},
        {{0, 4, 1}, eri3d1KernelFun<0, 4, 1>},
        {{4, 0, 2}, eri3d1KernelFun<4, 0, 2>},
        {{0, 4, 2}, eri3d1KernelFun<0, 4, 2>},
        {{4, 1, 0}, eri3d1KernelFun<4, 1, 0>},
        {{1, 4, 0}, eri3d1KernelFun<1, 4, 0>},
        {{4, 1, 1}, eri3d1KernelFun<4, 1, 1>},
        {{1, 4, 1}, eri3d1KernelFun<1, 4, 1>},
        {{4, 2, 0}, eri3d1KernelFun<4, 2, 0>},
        {{2, 4, 0}, eri3d1KernelFun<2, 4, 0>},
        {{5, 0, 0}, eri3d1KernelFun<5, 0, 0>},
        {{0, 5, 0}, eri3d1KernelFun<0, 5, 0>},
        {{5, 0, 1}, eri3d1KernelFun<5, 0, 1>},
        {{0, 5, 1}, eri3d1KernelFun<0, 5, 1>},
        {{5, 1, 0}, eri3d1KernelFun<5, 1, 0>},
        {{1, 5, 0}, eri3d1KernelFun<1, 5, 0>},
        {{6, 0, 0}, eri3d1KernelFun<6, 0, 0>},
        {{0, 6, 0}, eri3d1KernelFun<0, 6, 0>}
    };

    const std::map<std::tuple<int, int, int, int>, eri4d1_kernelfun_t> eri4d1_kernelfuns{
        {{0, 0, 0, 0}, eri4d1KernelFun<0, 0, 0, 0>},
        {{0, 0, 1, 0}, eri4d1KernelFun<0, 0, 1, 0>},
        {{0, 0, 0, 1}, eri4d1KernelFun<0, 0, 0, 1>},
        {{0, 0, 1, 1}, eri4d1KernelFun<0, 0, 1, 1>},
        {{0, 0, 2, 0}, eri4d1KernelFun<0, 0, 2, 0>},
        {{0, 0, 0, 2}, eri4d1KernelFun<0, 0, 0, 2>},
        {{0, 0, 2, 1}, eri4d1KernelFun<0, 0, 2, 1>},
        {{0, 0, 1, 2}, eri4d1KernelFun<0, 0, 1, 2>},
        {{0, 0, 2, 2}, eri4d1KernelFun<0, 0, 2, 2>},
        {{0, 0, 3, 0}, eri4d1KernelFun<0, 0, 3, 0>},
        {{0, 0, 0, 3}, eri4d1KernelFun<0, 0, 0, 3>},
        {{0, 0, 3, 1}, eri4d1KernelFun<0, 0, 3, 1>},
        {{0, 0, 1, 3}, eri4d1KernelFun<0, 0, 1, 3>},
        {{0, 0, 3, 2}, eri4d1KernelFun<0, 0, 3, 2>},
        {{0, 0, 2, 3}, eri4d1KernelFun<0, 0, 2, 3>},
        {{0, 0, 3, 3}, eri4d1KernelFun<0, 0, 3, 3>},
        {{0, 0, 4, 0}, eri4d1KernelFun<0, 0, 4, 0>},
        {{0, 0, 0, 4}, eri4d1KernelFun<0, 0, 0, 4>},
        {{0, 0, 4, 1}, eri4d1KernelFun<0, 0, 4, 1>},
        {{0, 0, 1, 4}, eri4d1KernelFun<0, 0, 1, 4>},
        {{0, 0, 4, 2}, eri4d1KernelFun<0, 0, 4, 2>},
        {{0, 0, 2, 4}, eri4d1KernelFun<0, 0, 2, 4>},
        {{0, 0, 5, 0}, eri4d1KernelFun<0, 0, 5, 0>},
        {{0, 0, 0, 5}, eri4d1KernelFun<0, 0, 0, 5>},
        {{0, 0, 5, 1}, eri4d1KernelFun<0, 0, 5, 1>},
        {{0, 0, 1, 5}, eri4d1KernelFun<0, 0, 1, 5>},
        {{0, 0, 6, 0}, eri4d1KernelFun<0, 0, 6, 0>},
        {{0, 0, 0, 6}, eri4d1KernelFun<0, 0, 0, 6>},
        {{1, 0, 0, 0}, eri4d1KernelFun<1, 0, 0, 0>},
        {{0, 1, 0, 0}, eri4d1KernelFun<0, 1, 0, 0>},
        {{1, 0, 1, 0}, eri4d1KernelFun<1, 0, 1, 0>},
        {{0, 1, 0, 1}, eri4d1KernelFun<0, 1, 0, 1>},
        {{1, 0, 0, 1}, eri4d1KernelFun<1, 0, 0, 1>},
        {{0, 1, 1, 0}, eri4d1KernelFun<0, 1, 1, 0>},
        {{1, 0, 1, 1}, eri4d1KernelFun<1, 0, 1, 1>},
        {{0, 1, 1, 1}, eri4d1KernelFun<0, 1, 1, 1>},
        {{1, 0, 2, 0}, eri4d1KernelFun<1, 0, 2, 0>},
        {{0, 1, 0, 2}, eri4d1KernelFun<0, 1, 0, 2>},
        {{1, 0, 0, 2}, eri4d1KernelFun<1, 0, 0, 2>},
        {{0, 1, 2, 0}, eri4d1KernelFun<0, 1, 2, 0>},
        {{1, 0, 2, 1}, eri4d1KernelFun<1, 0, 2, 1>},
        {{0, 1, 1, 2}, eri4d1KernelFun<0, 1, 1, 2>},
        {{1, 0, 1, 2}, eri4d1KernelFun<1, 0, 1, 2>},
        {{0, 1, 2, 1}, eri4d1KernelFun<0, 1, 2, 1>},
        {{1, 0, 2, 2}, eri4d1KernelFun<1, 0, 2, 2>},
        {{0, 1, 2, 2}, eri4d1KernelFun<0, 1, 2, 2>},
        {{1, 0, 3, 0}, eri4d1KernelFun<1, 0, 3, 0>},
        {{0, 1, 0, 3}, eri4d1KernelFun<0, 1, 0, 3>},
        {{1, 0, 0, 3}, eri4d1KernelFun<1, 0, 0, 3>},
        {{0, 1, 3, 0}, eri4d1KernelFun<0, 1, 3, 0>},
        {{1, 0, 3, 1}, eri4d1KernelFun<1, 0, 3, 1>},
        {{0, 1, 1, 3}, eri4d1KernelFun<0, 1, 1, 3>},
        {{1, 0, 1, 3}, eri4d1KernelFun<1, 0, 1, 3>},
        {{0, 1, 3, 1}, eri4d1KernelFun<0, 1, 3, 1>},
        {{1, 0, 3, 2}, eri4d1KernelFun<1, 0, 3, 2>},
        {{0, 1, 2, 3}, eri4d1KernelFun<0, 1, 2, 3>},
        {{1, 0, 2, 3}, eri4d1KernelFun<1, 0, 2, 3>},
        {{0, 1, 3, 2}, eri4d1KernelFun<0, 1, 3, 2>},
        {{1, 0, 4, 0}, eri4d1KernelFun<1, 0, 4, 0>},
        {{0, 1, 0, 4}, eri4d1KernelFun<0, 1, 0, 4>},
        {{1, 0, 0, 4}, eri4d1KernelFun<1, 0, 0, 4>},
        {{0, 1, 4, 0}, eri4d1KernelFun<0, 1, 4, 0>},
        {{1, 0, 4, 1}, eri4d1KernelFun<1, 0, 4, 1>},
        {{0, 1, 1, 4}, eri4d1KernelFun<0, 1, 1, 4>},
        {{1, 0, 1, 4}, eri4d1KernelFun<1, 0, 1, 4>},
        {{0, 1, 4, 1}, eri4d1KernelFun<0, 1, 4, 1>},
        {{1, 0, 5, 0}, eri4d1KernelFun<1, 0, 5, 0>},
        {{0, 1, 0, 5}, eri4d1KernelFun<0, 1, 0, 5>},
        {{1, 0, 0, 5}, eri4d1KernelFun<1, 0, 0, 5>},
        {{0, 1, 5, 0}, eri4d1KernelFun<0, 1, 5, 0>},
        {{1, 1, 0, 0}, eri4d1KernelFun<1, 1, 0, 0>},
        {{1, 1, 1, 0}, eri4d1KernelFun<1, 1, 1, 0>},
        {{1, 1, 0, 1}, eri4d1KernelFun<1, 1, 0, 1>},
        {{1, 1, 1, 1}, eri4d1KernelFun<1, 1, 1, 1>},
        {{1, 1, 2, 0}, eri4d1KernelFun<1, 1, 2, 0>},
        {{1, 1, 0, 2}, eri4d1KernelFun<1, 1, 0, 2>},
        {{1, 1, 2, 1}, eri4d1KernelFun<1, 1, 2, 1>},
        {{1, 1, 1, 2}, eri4d1KernelFun<1, 1, 1, 2>},
        {{1, 1, 2, 2}, eri4d1KernelFun<1, 1, 2, 2>},
        {{1, 1, 3, 0}, eri4d1KernelFun<1, 1, 3, 0>},
        {{1, 1, 0, 3}, eri4d1KernelFun<1, 1, 0, 3>},
        {{1, 1, 3, 1}, eri4d1KernelFun<1, 1, 3, 1>},
        {{1, 1, 1, 3}, eri4d1KernelFun<1, 1, 1, 3>},
        {{1, 1, 4, 0}, eri4d1KernelFun<1, 1, 4, 0>},
        {{1, 1, 0, 4}, eri4d1KernelFun<1, 1, 0, 4>},
        {{2, 0, 0, 0}, eri4d1KernelFun<2, 0, 0, 0>},
        {{0, 2, 0, 0}, eri4d1KernelFun<0, 2, 0, 0>},
        {{2, 0, 1, 0}, eri4d1KernelFun<2, 0, 1, 0>},
        {{0, 2, 0, 1}, eri4d1KernelFun<0, 2, 0, 1>},
        {{2, 0, 0, 1}, eri4d1KernelFun<2, 0, 0, 1>},
        {{0, 2, 1, 0}, eri4d1KernelFun<0, 2, 1, 0>},
        {{2, 0, 1, 1}, eri4d1KernelFun<2, 0, 1, 1>},
        {{0, 2, 1, 1}, eri4d1KernelFun<0, 2, 1, 1>},
        {{2, 0, 2, 0}, eri4d1KernelFun<2, 0, 2, 0>},
        {{0, 2, 0, 2}, eri4d1KernelFun<0, 2, 0, 2>},
        {{2, 0, 0, 2}, eri4d1KernelFun<2, 0, 0, 2>},
        {{0, 2, 2, 0}, eri4d1KernelFun<0, 2, 2, 0>},
        {{2, 0, 2, 1}, eri4d1KernelFun<2, 0, 2, 1>},
        {{0, 2, 1, 2}, eri4d1KernelFun<0, 2, 1, 2>},
        {{2, 0, 1, 2}, eri4d1KernelFun<2, 0, 1, 2>},
        {{0, 2, 2, 1}, eri4d1KernelFun<0, 2, 2, 1>},
        {{2, 0, 2, 2}, eri4d1KernelFun<2, 0, 2, 2>},
        {{0, 2, 2, 2}, eri4d1KernelFun<0, 2, 2, 2>},
        {{2, 0, 3, 0}, eri4d1KernelFun<2, 0, 3, 0>},
        {{0, 2, 0, 3}, eri4d1KernelFun<0, 2, 0, 3>},
        {{2, 0, 0, 3}, eri4d1KernelFun<2, 0, 0, 3>},
        {{0, 2, 3, 0}, eri4d1KernelFun<0, 2, 3, 0>},
        {{2, 0, 3, 1}, eri4d1KernelFun<2, 0, 3, 1>},
        {{0, 2, 1, 3}, eri4d1KernelFun<0, 2, 1, 3>},
        {{2, 0, 1, 3}, eri4d1KernelFun<2, 0, 1, 3>},
        {{0, 2, 3, 1}, eri4d1KernelFun<0, 2, 3, 1>},
        {{2, 0, 4, 0}, eri4d1KernelFun<2, 0, 4, 0>},
        {{0, 2, 0, 4}, eri4d1KernelFun<0, 2, 0, 4>},
        {{2, 0, 0, 4}, eri4d1KernelFun<2, 0, 0, 4>},
        {{0, 2, 4, 0}, eri4d1KernelFun<0, 2, 4, 0>},
        {{2, 1, 0, 0}, eri4d1KernelFun<2, 1, 0, 0>},
        {{1, 2, 0, 0}, eri4d1KernelFun<1, 2, 0, 0>},
        {{2, 1, 1, 0}, eri4d1KernelFun<2, 1, 1, 0>},
        {{1, 2, 0, 1}, eri4d1KernelFun<1, 2, 0, 1>},
        {{2, 1, 0, 1}, eri4d1KernelFun<2, 1, 0, 1>},
        {{1, 2, 1, 0}, eri4d1KernelFun<1, 2, 1, 0>},
        {{2, 1, 1, 1}, eri4d1KernelFun<2, 1, 1, 1>},
        {{1, 2, 1, 1}, eri4d1KernelFun<1, 2, 1, 1>},
        {{2, 1, 2, 0}, eri4d1KernelFun<2, 1, 2, 0>},
        {{1, 2, 0, 2}, eri4d1KernelFun<1, 2, 0, 2>},
        {{2, 1, 0, 2}, eri4d1KernelFun<2, 1, 0, 2>},
        {{1, 2, 2, 0}, eri4d1KernelFun<1, 2, 2, 0>},
        {{2, 1, 2, 1}, eri4d1KernelFun<2, 1, 2, 1>},
        {{1, 2, 1, 2}, eri4d1KernelFun<1, 2, 1, 2>},
        {{2, 1, 1, 2}, eri4d1KernelFun<2, 1, 1, 2>},
        {{1, 2, 2, 1}, eri4d1KernelFun<1, 2, 2, 1>},
        {{2, 1, 3, 0}, eri4d1KernelFun<2, 1, 3, 0>},
        {{1, 2, 0, 3}, eri4d1KernelFun<1, 2, 0, 3>},
        {{2, 1, 0, 3}, eri4d1KernelFun<2, 1, 0, 3>},
        {{1, 2, 3, 0}, eri4d1KernelFun<1, 2, 3, 0>},
        {{2, 2, 0, 0}, eri4d1KernelFun<2, 2, 0, 0>},
        {{2, 2, 1, 0}, eri4d1KernelFun<2, 2, 1, 0>},
        {{2, 2, 0, 1}, eri4d1KernelFun<2, 2, 0, 1>},
        {{2, 2, 1, 1}, eri4d1KernelFun<2, 2, 1, 1>},
        {{2, 2, 2, 0}, eri4d1KernelFun<2, 2, 2, 0>},
        {{2, 2, 0, 2}, eri4d1KernelFun<2, 2, 0, 2>},
        {{3, 0, 0, 0}, eri4d1KernelFun<3, 0, 0, 0>},
        {{0, 3, 0, 0}, eri4d1KernelFun<0, 3, 0, 0>},
        {{3, 0, 1, 0}, eri4d1KernelFun<3, 0, 1, 0>},
        {{0, 3, 0, 1}, eri4d1KernelFun<0, 3, 0, 1>},
        {{3, 0, 0, 1}, eri4d1KernelFun<3, 0, 0, 1>},
        {{0, 3, 1, 0}, eri4d1KernelFun<0, 3, 1, 0>},
        {{3, 0, 1, 1}, eri4d1KernelFun<3, 0, 1, 1>},
        {{0, 3, 1, 1}, eri4d1KernelFun<0, 3, 1, 1>},
        {{3, 0, 2, 0}, eri4d1KernelFun<3, 0, 2, 0>},
        {{0, 3, 0, 2}, eri4d1KernelFun<0, 3, 0, 2>},
        {{3, 0, 0, 2}, eri4d1KernelFun<3, 0, 0, 2>},
        {{0, 3, 2, 0}, eri4d1KernelFun<0, 3, 2, 0>},
        {{3, 0, 2, 1}, eri4d1KernelFun<3, 0, 2, 1>},
        {{0, 3, 1, 2}, eri4d1KernelFun<0, 3, 1, 2>},
        {{3, 0, 1, 2}, eri4d1KernelFun<3, 0, 1, 2>},
        {{0, 3, 2, 1}, eri4d1KernelFun<0, 3, 2, 1>},
        {{3, 0, 3, 0}, eri4d1KernelFun<3, 0, 3, 0>},
        {{0, 3, 0, 3}, eri4d1KernelFun<0, 3, 0, 3>},
        {{3, 0, 0, 3}, eri4d1KernelFun<3, 0, 0, 3>},
        {{0, 3, 3, 0}, eri4d1KernelFun<0, 3, 3, 0>},
        {{3, 1, 0, 0}, eri4d1KernelFun<3, 1, 0, 0>},
        {{1, 3, 0, 0}, eri4d1KernelFun<1, 3, 0, 0>},
        {{3, 1, 1, 0}, eri4d1KernelFun<3, 1, 1, 0>},
        {{1, 3, 0, 1}, eri4d1KernelFun<1, 3, 0, 1>},
        {{3, 1, 0, 1}, eri4d1KernelFun<3, 1, 0, 1>},
        {{1, 3, 1, 0}, eri4d1KernelFun<1, 3, 1, 0>},
        {{3, 1, 1, 1}, eri4d1KernelFun<3, 1, 1, 1>},
        {{1, 3, 1, 1}, eri4d1KernelFun<1, 3, 1, 1>},
        {{3, 1, 2, 0}, eri4d1KernelFun<3, 1, 2, 0>},
        {{1, 3, 0, 2}, eri4d1KernelFun<1, 3, 0, 2>},
        {{3, 1, 0, 2}, eri4d1KernelFun<3, 1, 0, 2>},
        {{1, 3, 2, 0}, eri4d1KernelFun<1, 3, 2, 0>},
        {{3, 2, 0, 0}, eri4d1KernelFun<3, 2, 0, 0>},
        {{2, 3, 0, 0}, eri4d1KernelFun<2, 3, 0, 0>},
        {{3, 2, 1, 0}, eri4d1KernelFun<3, 2, 1, 0>},
        {{2, 3, 0, 1}, eri4d1KernelFun<2, 3, 0, 1>},
        {{3, 2, 0, 1}, eri4d1KernelFun<3, 2, 0, 1>},
        {{2, 3, 1, 0}, eri4d1KernelFun<2, 3, 1, 0>},
        {{3, 3, 0, 0}, eri4d1KernelFun<3, 3, 0, 0>},
        {{4, 0, 0, 0}, eri4d1KernelFun<4, 0, 0, 0>},
        {{0, 4, 0, 0}, eri4d1KernelFun<0, 4, 0, 0>},
        {{4, 0, 1, 0}, eri4d1KernelFun<4, 0, 1, 0>},
        {{0, 4, 0, 1}, eri4d1KernelFun<0, 4, 0, 1>},
        {{4, 0, 0, 1}, eri4d1KernelFun<4, 0, 0, 1>},
        {{0, 4, 1, 0}, eri4d1KernelFun<0, 4, 1, 0>},
        {{4, 0, 1, 1}, eri4d1KernelFun<4, 0, 1, 1>},
        {{0, 4, 1, 1}, eri4d1KernelFun<0, 4, 1, 1>},
        {{4, 0, 2, 0}, eri4d1KernelFun<4, 0, 2, 0>},
        {{0, 4, 0, 2}, eri4d1KernelFun<0, 4, 0, 2>},
        {{4, 0, 0, 2}, eri4d1KernelFun<4, 0, 0, 2>},
        {{0, 4, 2, 0}, eri4d1KernelFun<0, 4, 2, 0>},
        {{4, 1, 0, 0}, eri4d1KernelFun<4, 1, 0, 0>},
        {{1, 4, 0, 0}, eri4d1KernelFun<1, 4, 0, 0>},
        {{4, 1, 1, 0}, eri4d1KernelFun<4, 1, 1, 0>},
        {{1, 4, 0, 1}, eri4d1KernelFun<1, 4, 0, 1>},
        {{4, 1, 0, 1}, eri4d1KernelFun<4, 1, 0, 1>},
        {{1, 4, 1, 0}, eri4d1KernelFun<1, 4, 1, 0>},
        {{4, 2, 0, 0}, eri4d1KernelFun<4, 2, 0, 0>},
        {{2, 4, 0, 0}, eri4d1KernelFun<2, 4, 0, 0>},
        {{5, 0, 0, 0}, eri4d1KernelFun<5, 0, 0, 0>},
        {{0, 5, 0, 0}, eri4d1KernelFun<0, 5, 0, 0>},
        {{5, 0, 1, 0}, eri4d1KernelFun<5, 0, 1, 0>},
        {{0, 5, 0, 1}, eri4d1KernelFun<0, 5, 0, 1>},
        {{5, 0, 0, 1}, eri4d1KernelFun<5, 0, 0, 1>},
        {{0, 5, 1, 0}, eri4d1KernelFun<0, 5, 1, 0>},
        {{5, 1, 0, 0}, eri4d1KernelFun<5, 1, 0, 0>},
        {{1, 5, 0, 0}, eri4d1KernelFun<1, 5, 0, 0>},
        {{6, 0, 0, 0}, eri4d1KernelFun<6, 0, 0, 0>},
        {{0, 6, 0, 0}, eri4d1KernelFun<0, 6, 0, 0>}
    };

    const std::map<std::tuple<int, int, int, int>, eri4soc_kernelfun_t> eri4soc_kernelfuns{
        {{0, 0, 0, 0}, eri4socKernelFun<0, 0, 0, 0>},
        {{0, 0, 1, 0}, eri4socKernelFun<0, 0, 1, 0>},
        {{0, 0, 0, 1}, eri4socKernelFun<0, 0, 0, 1>},
        {{0, 0, 1, 1}, eri4socKernelFun<0, 0, 1, 1>},
        {{0, 0, 2, 0}, eri4socKernelFun<0, 0, 2, 0>},
        {{0, 0, 0, 2}, eri4socKernelFun<0, 0, 0, 2>},
        {{0, 0, 2, 1}, eri4socKernelFun<0, 0, 2, 1>},
        {{0, 0, 1, 2}, eri4socKernelFun<0, 0, 1, 2>},
        {{0, 0, 2, 2}, eri4socKernelFun<0, 0, 2, 2>},
        {{0, 0, 3, 0}, eri4socKernelFun<0, 0, 3, 0>},
        {{0, 0, 0, 3}, eri4socKernelFun<0, 0, 0, 3>},
        {{0, 0, 3, 1}, eri4socKernelFun<0, 0, 3, 1>},
        {{0, 0, 1, 3}, eri4socKernelFun<0, 0, 1, 3>},
        {{0, 0, 3, 2}, eri4socKernelFun<0, 0, 3, 2>},
        {{0, 0, 2, 3}, eri4socKernelFun<0, 0, 2, 3>},
        {{0, 0, 3, 3}, eri4socKernelFun<0, 0, 3, 3>},
        {{0, 0, 4, 0}, eri4socKernelFun<0, 0, 4, 0>},
        {{0, 0, 0, 4}, eri4socKernelFun<0, 0, 0, 4>},
        {{0, 0, 4, 1}, eri4socKernelFun<0, 0, 4, 1>},
        {{0, 0, 1, 4}, eri4socKernelFun<0, 0, 1, 4>},
        {{0, 0, 4, 2}, eri4socKernelFun<0, 0, 4, 2>},
        {{0, 0, 2, 4}, eri4socKernelFun<0, 0, 2, 4>},
        {{0, 0, 5, 0}, eri4socKernelFun<0, 0, 5, 0>},
        {{0, 0, 0, 5}, eri4socKernelFun<0, 0, 0, 5>},
        {{0, 0, 5, 1}, eri4socKernelFun<0, 0, 5, 1>},
        {{0, 0, 1, 5}, eri4socKernelFun<0, 0, 1, 5>},
        {{0, 0, 6, 0}, eri4socKernelFun<0, 0, 6, 0>},
        {{0, 0, 0, 6}, eri4socKernelFun<0, 0, 0, 6>},
        {{1, 0, 0, 0}, eri4socKernelFun<1, 0, 0, 0>},
        {{0, 1, 0, 0}, eri4socKernelFun<0, 1, 0, 0>},
        {{1, 0, 1, 0}, eri4socKernelFun<1, 0, 1, 0>},
        {{0, 1, 0, 1}, eri4socKernelFun<0, 1, 0, 1>},
        {{1, 0, 0, 1}, eri4socKernelFun<1, 0, 0, 1>},
        {{0, 1, 1, 0}, eri4socKernelFun<0, 1, 1, 0>},
        {{1, 0, 1, 1}, eri4socKernelFun<1, 0, 1, 1>},
        {{0, 1, 1, 1}, eri4socKernelFun<0, 1, 1, 1>},
        {{1, 0, 2, 0}, eri4socKernelFun<1, 0, 2, 0>},
        {{0, 1, 0, 2}, eri4socKernelFun<0, 1, 0, 2>},
        {{1, 0, 0, 2}, eri4socKernelFun<1, 0, 0, 2>},
        {{0, 1, 2, 0}, eri4socKernelFun<0, 1, 2, 0>},
        {{1, 0, 2, 1}, eri4socKernelFun<1, 0, 2, 1>},
        {{0, 1, 1, 2}, eri4socKernelFun<0, 1, 1, 2>},
        {{1, 0, 1, 2}, eri4socKernelFun<1, 0, 1, 2>},
        {{0, 1, 2, 1}, eri4socKernelFun<0, 1, 2, 1>},
        {{1, 0, 2, 2}, eri4socKernelFun<1, 0, 2, 2>},
        {{0, 1, 2, 2}, eri4socKernelFun<0, 1, 2, 2>},
        {{1, 0, 3, 0}, eri4socKernelFun<1, 0, 3, 0>},
        {{0, 1, 0, 3}, eri4socKernelFun<0, 1, 0, 3>},
        {{1, 0, 0, 3}, eri4socKernelFun<1, 0, 0, 3>},
        {{0, 1, 3, 0}, eri4socKernelFun<0, 1, 3, 0>},
        {{1, 0, 3, 1}, eri4socKernelFun<1, 0, 3, 1>},
        {{0, 1, 1, 3}, eri4socKernelFun<0, 1, 1, 3>},
        {{1, 0, 1, 3}, eri4socKernelFun<1, 0, 1, 3>},
        {{0, 1, 3, 1}, eri4socKernelFun<0, 1, 3, 1>},
        {{1, 0, 3, 2}, eri4socKernelFun<1, 0, 3, 2>},
        {{0, 1, 2, 3}, eri4socKernelFun<0, 1, 2, 3>},
        {{1, 0, 2, 3}, eri4socKernelFun<1, 0, 2, 3>},
        {{0, 1, 3, 2}, eri4socKernelFun<0, 1, 3, 2>},
        {{1, 0, 4, 0}, eri4socKernelFun<1, 0, 4, 0>},
        {{0, 1, 0, 4}, eri4socKernelFun<0, 1, 0, 4>},
        {{1, 0, 0, 4}, eri4socKernelFun<1, 0, 0, 4>},
        {{0, 1, 4, 0}, eri4socKernelFun<0, 1, 4, 0>},
        {{1, 0, 4, 1}, eri4socKernelFun<1, 0, 4, 1>},
        {{0, 1, 1, 4}, eri4socKernelFun<0, 1, 1, 4>},
        {{1, 0, 1, 4}, eri4socKernelFun<1, 0, 1, 4>},
        {{0, 1, 4, 1}, eri4socKernelFun<0, 1, 4, 1>},
        {{1, 0, 5, 0}, eri4socKernelFun<1, 0, 5, 0>},
        {{0, 1, 0, 5}, eri4socKernelFun<0, 1, 0, 5>},
        {{1, 0, 0, 5}, eri4socKernelFun<1, 0, 0, 5>},
        {{0, 1, 5, 0}, eri4socKernelFun<0, 1, 5, 0>},
        {{1, 1, 0, 0}, eri4socKernelFun<1, 1, 0, 0>},
        {{1, 1, 1, 0}, eri4socKernelFun<1, 1, 1, 0>},
        {{1, 1, 0, 1}, eri4socKernelFun<1, 1, 0, 1>},
        {{1, 1, 1, 1}, eri4socKernelFun<1, 1, 1, 1>},
        {{1, 1, 2, 0}, eri4socKernelFun<1, 1, 2, 0>},
        {{1, 1, 0, 2}, eri4socKernelFun<1, 1, 0, 2>},
        {{1, 1, 2, 1}, eri4socKernelFun<1, 1, 2, 1>},
        {{1, 1, 1, 2}, eri4socKernelFun<1, 1, 1, 2>},
        {{1, 1, 2, 2}, eri4socKernelFun<1, 1, 2, 2>},
        {{1, 1, 3, 0}, eri4socKernelFun<1, 1, 3, 0>},
        {{1, 1, 0, 3}, eri4socKernelFun<1, 1, 0, 3>},
        {{1, 1, 3, 1}, eri4socKernelFun<1, 1, 3, 1>},
        {{1, 1, 1, 3}, eri4socKernelFun<1, 1, 1, 3>},
        {{1, 1, 4, 0}, eri4socKernelFun<1, 1, 4, 0>},
        {{1, 1, 0, 4}, eri4socKernelFun<1, 1, 0, 4>},
        {{2, 0, 0, 0}, eri4socKernelFun<2, 0, 0, 0>},
        {{0, 2, 0, 0}, eri4socKernelFun<0, 2, 0, 0>},
        {{2, 0, 1, 0}, eri4socKernelFun<2, 0, 1, 0>},
        {{0, 2, 0, 1}, eri4socKernelFun<0, 2, 0, 1>},
        {{2, 0, 0, 1}, eri4socKernelFun<2, 0, 0, 1>},
        {{0, 2, 1, 0}, eri4socKernelFun<0, 2, 1, 0>},
        {{2, 0, 1, 1}, eri4socKernelFun<2, 0, 1, 1>},
        {{0, 2, 1, 1}, eri4socKernelFun<0, 2, 1, 1>},
        {{2, 0, 2, 0}, eri4socKernelFun<2, 0, 2, 0>},
        {{0, 2, 0, 2}, eri4socKernelFun<0, 2, 0, 2>},
        {{2, 0, 0, 2}, eri4socKernelFun<2, 0, 0, 2>},
        {{0, 2, 2, 0}, eri4socKernelFun<0, 2, 2, 0>},
        {{2, 0, 2, 1}, eri4socKernelFun<2, 0, 2, 1>},
        {{0, 2, 1, 2}, eri4socKernelFun<0, 2, 1, 2>},
        {{2, 0, 1, 2}, eri4socKernelFun<2, 0, 1, 2>},
        {{0, 2, 2, 1}, eri4socKernelFun<0, 2, 2, 1>},
        {{2, 0, 2, 2}, eri4socKernelFun<2, 0, 2, 2>},
        {{0, 2, 2, 2}, eri4socKernelFun<0, 2, 2, 2>},
        {{2, 0, 3, 0}, eri4socKernelFun<2, 0, 3, 0>},
        {{0, 2, 0, 3}, eri4socKernelFun<0, 2, 0, 3>},
        {{2, 0, 0, 3}, eri4socKernelFun<2, 0, 0, 3>},
        {{0, 2, 3, 0}, eri4socKernelFun<0, 2, 3, 0>},
        {{2, 0, 3, 1}, eri4socKernelFun<2, 0, 3, 1>},
        {{0, 2, 1, 3}, eri4socKernelFun<0, 2, 1, 3>},
        {{2, 0, 1, 3}, eri4socKernelFun<2, 0, 1, 3>},
        {{0, 2, 3, 1}, eri4socKernelFun<0, 2, 3, 1>},
        {{2, 0, 4, 0}, eri4socKernelFun<2, 0, 4, 0>},
        {{0, 2, 0, 4}, eri4socKernelFun<0, 2, 0, 4>},
        {{2, 0, 0, 4}, eri4socKernelFun<2, 0, 0, 4>},
        {{0, 2, 4, 0}, eri4socKernelFun<0, 2, 4, 0>},
        {{2, 1, 0, 0}, eri4socKernelFun<2, 1, 0, 0>},
        {{1, 2, 0, 0}, eri4socKernelFun<1, 2, 0, 0>},
        {{2, 1, 1, 0}, eri4socKernelFun<2, 1, 1, 0>},
        {{1, 2, 0, 1}, eri4socKernelFun<1, 2, 0, 1>},
        {{2, 1, 0, 1}, eri4socKernelFun<2, 1, 0, 1>},
        {{1, 2, 1, 0}, eri4socKernelFun<1, 2, 1, 0>},
        {{2, 1, 1, 1}, eri4socKernelFun<2, 1, 1, 1>},
        {{1, 2, 1, 1}, eri4socKernelFun<1, 2, 1, 1>},
        {{2, 1, 2, 0}, eri4socKernelFun<2, 1, 2, 0>},
        {{1, 2, 0, 2}, eri4socKernelFun<1, 2, 0, 2>},
        {{2, 1, 0, 2}, eri4socKernelFun<2, 1, 0, 2>},
        {{1, 2, 2, 0}, eri4socKernelFun<1, 2, 2, 0>},
        {{2, 1, 2, 1}, eri4socKernelFun<2, 1, 2, 1>},
        {{1, 2, 1, 2}, eri4socKernelFun<1, 2, 1, 2>},
        {{2, 1, 1, 2}, eri4socKernelFun<2, 1, 1, 2>},
        {{1, 2, 2, 1}, eri4socKernelFun<1, 2, 2, 1>},
        {{2, 1, 3, 0}, eri4socKernelFun<2, 1, 3, 0>},
        {{1, 2, 0, 3}, eri4socKernelFun<1, 2, 0, 3>},
        {{2, 1, 0, 3}, eri4socKernelFun<2, 1, 0, 3>},
        {{1, 2, 3, 0}, eri4socKernelFun<1, 2, 3, 0>},
        {{2, 2, 0, 0}, eri4socKernelFun<2, 2, 0, 0>},
        {{2, 2, 1, 0}, eri4socKernelFun<2, 2, 1, 0>},
        {{2, 2, 0, 1}, eri4socKernelFun<2, 2, 0, 1>},
        {{2, 2, 1, 1}, eri4socKernelFun<2, 2, 1, 1>},
        {{2, 2, 2, 0}, eri4socKernelFun<2, 2, 2, 0>},
        {{2, 2, 0, 2}, eri4socKernelFun<2, 2, 0, 2>},
        {{3, 0, 0, 0}, eri4socKernelFun<3, 0, 0, 0>},
        {{0, 3, 0, 0}, eri4socKernelFun<0, 3, 0, 0>},
        {{3, 0, 1, 0}, eri4socKernelFun<3, 0, 1, 0>},
        {{0, 3, 0, 1}, eri4socKernelFun<0, 3, 0, 1>},
        {{3, 0, 0, 1}, eri4socKernelFun<3, 0, 0, 1>},
        {{0, 3, 1, 0}, eri4socKernelFun<0, 3, 1, 0>},
        {{3, 0, 1, 1}, eri4socKernelFun<3, 0, 1, 1>},
        {{0, 3, 1, 1}, eri4socKernelFun<0, 3, 1, 1>},
        {{3, 0, 2, 0}, eri4socKernelFun<3, 0, 2, 0>},
        {{0, 3, 0, 2}, eri4socKernelFun<0, 3, 0, 2>},
        {{3, 0, 0, 2}, eri4socKernelFun<3, 0, 0, 2>},
        {{0, 3, 2, 0}, eri4socKernelFun<0, 3, 2, 0>},
        {{3, 0, 2, 1}, eri4socKernelFun<3, 0, 2, 1>},
        {{0, 3, 1, 2}, eri4socKernelFun<0, 3, 1, 2>},
        {{3, 0, 1, 2}, eri4socKernelFun<3, 0, 1, 2>},
        {{0, 3, 2, 1}, eri4socKernelFun<0, 3, 2, 1>},
        {{3, 0, 3, 0}, eri4socKernelFun<3, 0, 3, 0>},
        {{0, 3, 0, 3}, eri4socKernelFun<0, 3, 0, 3>},
        {{3, 0, 0, 3}, eri4socKernelFun<3, 0, 0, 3>},
        {{0, 3, 3, 0}, eri4socKernelFun<0, 3, 3, 0>},
        {{3, 1, 0, 0}, eri4socKernelFun<3, 1, 0, 0>},
        {{1, 3, 0, 0}, eri4socKernelFun<1, 3, 0, 0>},
        {{3, 1, 1, 0}, eri4socKernelFun<3, 1, 1, 0>},
        {{1, 3, 0, 1}, eri4socKernelFun<1, 3, 0, 1>},
        {{3, 1, 0, 1}, eri4socKernelFun<3, 1, 0, 1>},
        {{1, 3, 1, 0}, eri4socKernelFun<1, 3, 1, 0>},
        {{3, 1, 1, 1}, eri4socKernelFun<3, 1, 1, 1>},
        {{1, 3, 1, 1}, eri4socKernelFun<1, 3, 1, 1>},
        {{3, 1, 2, 0}, eri4socKernelFun<3, 1, 2, 0>},
        {{1, 3, 0, 2}, eri4socKernelFun<1, 3, 0, 2>},
        {{3, 1, 0, 2}, eri4socKernelFun<3, 1, 0, 2>},
        {{1, 3, 2, 0}, eri4socKernelFun<1, 3, 2, 0>},
        {{3, 2, 0, 0}, eri4socKernelFun<3, 2, 0, 0>},
        {{2, 3, 0, 0}, eri4socKernelFun<2, 3, 0, 0>},
        {{3, 2, 1, 0}, eri4socKernelFun<3, 2, 1, 0>},
        {{2, 3, 0, 1}, eri4socKernelFun<2, 3, 0, 1>},
        {{3, 2, 0, 1}, eri4socKernelFun<3, 2, 0, 1>},
        {{2, 3, 1, 0}, eri4socKernelFun<2, 3, 1, 0>},
        {{3, 3, 0, 0}, eri4socKernelFun<3, 3, 0, 0>},
        {{4, 0, 0, 0}, eri4socKernelFun<4, 0, 0, 0>},
        {{0, 4, 0, 0}, eri4socKernelFun<0, 4, 0, 0>},
        {{4, 0, 1, 0}, eri4socKernelFun<4, 0, 1, 0>},
        {{0, 4, 0, 1}, eri4socKernelFun<0, 4, 0, 1>},
        {{4, 0, 0, 1}, eri4socKernelFun<4, 0, 0, 1>},
        {{0, 4, 1, 0}, eri4socKernelFun<0, 4, 1, 0>},
        {{4, 0, 1, 1}, eri4socKernelFun<4, 0, 1, 1>},
        {{0, 4, 1, 1}, eri4socKernelFun<0, 4, 1, 1>},
        {{4, 0, 2, 0}, eri4socKernelFun<4, 0, 2, 0>},
        {{0, 4, 0, 2}, eri4socKernelFun<0, 4, 0, 2>},
        {{4, 0, 0, 2}, eri4socKernelFun<4, 0, 0, 2>},
        {{0, 4, 2, 0}, eri4socKernelFun<0, 4, 2, 0>},
        {{4, 1, 0, 0}, eri4socKernelFun<4, 1, 0, 0>},
        {{1, 4, 0, 0}, eri4socKernelFun<1, 4, 0, 0>},
        {{4, 1, 1, 0}, eri4socKernelFun<4, 1, 1, 0>},
        {{1, 4, 0, 1}, eri4socKernelFun<1, 4, 0, 1>},
        {{4, 1, 0, 1}, eri4socKernelFun<4, 1, 0, 1>},
        {{1, 4, 1, 0}, eri4socKernelFun<1, 4, 1, 0>},
        {{4, 2, 0, 0}, eri4socKernelFun<4, 2, 0, 0>},
        {{2, 4, 0, 0}, eri4socKernelFun<2, 4, 0, 0>},
        {{5, 0, 0, 0}, eri4socKernelFun<5, 0, 0, 0>},
        {{0, 5, 0, 0}, eri4socKernelFun<0, 5, 0, 0>},
        {{5, 0, 1, 0}, eri4socKernelFun<5, 0, 1, 0>},
        {{0, 5, 0, 1}, eri4socKernelFun<0, 5, 0, 1>},
        {{5, 0, 0, 1}, eri4socKernelFun<5, 0, 0, 1>},
        {{0, 5, 1, 0}, eri4socKernelFun<0, 5, 1, 0>},
        {{5, 1, 0, 0}, eri4socKernelFun<5, 1, 0, 0>},
        {{1, 5, 0, 0}, eri4socKernelFun<1, 5, 0, 0>},
        {{6, 0, 0, 0}, eri4socKernelFun<6, 0, 0, 0>},
        {{0, 6, 0, 0}, eri4socKernelFun<0, 6, 0, 0>}
    };

    const std::map<std::tuple<int, int, int>, eri3soc_kernelfun_t> eri3soc_kernelfuns{
        {{0, 0, 0}, eri3socKernelFun<0, 0, 0>},
        {{0, 0, 1}, eri3socKernelFun<0, 0, 1>},
        {{0, 0, 2}, eri3socKernelFun<0, 0, 2>},
        {{0, 0, 3}, eri3socKernelFun<0, 0, 3>},
        {{0, 0, 4}, eri3socKernelFun<0, 0, 4>},
        {{0, 0, 5}, eri3socKernelFun<0, 0, 5>},
        {{0, 0, 6}, eri3socKernelFun<0, 0, 6>},
        {{1, 0, 0}, eri3socKernelFun<1, 0, 0>},
        {{0, 1, 0}, eri3socKernelFun<0, 1, 0>},
        {{1, 0, 1}, eri3socKernelFun<1, 0, 1>},
        {{0, 1, 1}, eri3socKernelFun<0, 1, 1>},
        {{1, 0, 2}, eri3socKernelFun<1, 0, 2>},
        {{0, 1, 2}, eri3socKernelFun<0, 1, 2>},
        {{1, 0, 3}, eri3socKernelFun<1, 0, 3>},
        {{0, 1, 3}, eri3socKernelFun<0, 1, 3>},
        {{1, 0, 4}, eri3socKernelFun<1, 0, 4>},
        {{0, 1, 4}, eri3socKernelFun<0, 1, 4>},
        {{1, 0, 5}, eri3socKernelFun<1, 0, 5>},
        {{0, 1, 5}, eri3socKernelFun<0, 1, 5>},
        {{1, 1, 0}, eri3socKernelFun<1, 1, 0>},
        {{1, 1, 1}, eri3socKernelFun<1, 1, 1>},
        {{1, 1, 2}, eri3socKernelFun<1, 1, 2>},
        {{1, 1, 3}, eri3socKernelFun<1, 1, 3>},
        {{1, 1, 4}, eri3socKernelFun<1, 1, 4>},
        {{2, 0, 0}, eri3socKernelFun<2, 0, 0>},
        {{0, 2, 0}, eri3socKernelFun<0, 2, 0>},
        {{2, 0, 1}, eri3socKernelFun<2, 0, 1>},
        {{0, 2, 1}, eri3socKernelFun<0, 2, 1>},
        {{2, 0, 2}, eri3socKernelFun<2, 0, 2>},
        {{0, 2, 2}, eri3socKernelFun<0, 2, 2>},
        {{2, 0, 3}, eri3socKernelFun<2, 0, 3>},
        {{0, 2, 3}, eri3socKernelFun<0, 2, 3>},
        {{2, 0, 4}, eri3socKernelFun<2, 0, 4>},
        {{0, 2, 4}, eri3socKernelFun<0, 2, 4>},
        {{2, 1, 0}, eri3socKernelFun<2, 1, 0>},
        {{1, 2, 0}, eri3socKernelFun<1, 2, 0>},
        {{2, 1, 1}, eri3socKernelFun<2, 1, 1>},
        {{1, 2, 1}, eri3socKernelFun<1, 2, 1>},
        {{2, 1, 2}, eri3socKernelFun<2, 1, 2>},
        {{1, 2, 2}, eri3socKernelFun<1, 2, 2>},
        {{2, 1, 3}, eri3socKernelFun<2, 1, 3>},
        {{1, 2, 3}, eri3socKernelFun<1, 2, 3>},
        {{2, 2, 0}, eri3socKernelFun<2, 2, 0>},
        {{2, 2, 1}, eri3socKernelFun<2, 2, 1>},
        {{2, 2, 2}, eri3socKernelFun<2, 2, 2>},
        {{3, 0, 0}, eri3socKernelFun<3, 0, 0>},
        {{0, 3, 0}, eri3socKernelFun<0, 3, 0>},
        {{3, 0, 1}, eri3socKernelFun<3, 0, 1>},
        {{0, 3, 1}, eri3socKernelFun<0, 3, 1>},
        {{3, 0, 2}, eri3socKernelFun<3, 0, 2>},
        {{0, 3, 2}, eri3socKernelFun<0, 3, 2>},
        {{3, 0, 3}, eri3socKernelFun<3, 0, 3>},
        {{0, 3, 3}, eri3socKernelFun<0, 3, 3>},
        {{3, 1, 0}, eri3socKernelFun<3, 1, 0>},
        {{1, 3, 0}, eri3socKernelFun<1, 3, 0>},
        {{3, 1, 1}, eri3socKernelFun<3, 1, 1>},
        {{1, 3, 1}, eri3socKernelFun<1, 3, 1>},
        {{3, 1, 2}, eri3socKernelFun<3, 1, 2>},
        {{1, 3, 2}, eri3socKernelFun<1, 3, 2>},
        {{3, 2, 0}, eri3socKernelFun<3, 2, 0>},
        {{2, 3, 0}, eri3socKernelFun<2, 3, 0>},
        {{3, 2, 1}, eri3socKernelFun<3, 2, 1>},
        {{2, 3, 1}, eri3socKernelFun<2, 3, 1>},
        {{3, 3, 0}, eri3socKernelFun<3, 3, 0>},
        {{4, 0, 0}, eri3socKernelFun<4, 0, 0>},
        {{0, 4, 0}, eri3socKernelFun<0, 4, 0>},
        {{4, 0, 1}, eri3socKernelFun<4, 0, 1>},
        {{0, 4, 1}, eri3socKernelFun<0, 4, 1>},
        {{4, 0, 2}, eri3socKernelFun<4, 0, 2>},
        {{0, 4, 2}, eri3socKernelFun<0, 4, 2>},
        {{4, 1, 0}, eri3socKernelFun<4, 1, 0>},
        {{1, 4, 0}, eri3socKernelFun<1, 4, 0>},
        {{4, 1, 1}, eri3socKernelFun<4, 1, 1>},
        {{1, 4, 1}, eri3socKernelFun<1, 4, 1>},
        {{4, 2, 0}, eri3socKernelFun<4, 2, 0>},
        {{2, 4, 0}, eri3socKernelFun<2, 4, 0>},
        {{5, 0, 0}, eri3socKernelFun<5, 0, 0>},
        {{0, 5, 0}, eri3socKernelFun<0, 5, 0>},
        {{5, 0, 1}, eri3socKernelFun<5, 0, 1>},
        {{0, 5, 1}, eri3socKernelFun<0, 5, 1>},
        {{5, 1, 0}, eri3socKernelFun<5, 1, 0>},
        {{1, 5, 0}, eri3socKernelFun<1, 5, 0>},
        {{6, 0, 0}, eri3socKernelFun<6, 0, 0>},
        {{0, 6, 0}, eri3socKernelFun<0, 6, 0>}
    };
}

lints::ERI4Kernel::ERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd)
{
    ecoeffs_bra_ = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs_ket_ = ecoeffsSHARK(sp_data_cd, true);

    auto [la, lb] = sp_data_ab.getLPair();
    auto [lc, ld] = sp_data_cd.getLPair();
    int labcd = la + lb + lc + ld;
    boys_grid_ = BoysGrid(labcd);

    if (labcd <= _max_l_rollout_)
        eri4_kernelfun_ = eri4_kernelfuns.at({la, lb, lc, ld});
    else
        eri4_kernelfun_ = [](const size_t ipair_ab, const size_t ipair_cd,
                             const ShellPairData &spd_ab, const ShellPairData &spd_cd,
                             const ERI4Kernel *eri4_kernel) -> vec4d
        {
            return eri4KernelFun(ipair_ab, ipair_cd, spd_ab, spd_cd, eri4_kernel);
        };
}

lints::ERI3Kernel::ERI3Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c)
{
    ecoeffs_bra_ = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs_ket_ = ecoeffsSHARK(sh_data_c, true);

    auto [la, lb] = sp_data_ab.getLPair();
    int lc = sh_data_c.l_;
    int labc = la + lb + lc;
    boys_grid_ = BoysGrid(labc);

    if (labc <= _max_l_rollout_)
        eri3_kernelfun_ = eri3_kernelfuns.at({la, lb, lc});
    else
        eri3_kernelfun_ = [](const size_t ipair_ab, const size_t ish_c,
                             const ShellPairData &spd_ab, const ShellData &shd_c,
                             const ERI3Kernel *eri3_kernel) -> vec3d
                          {
                              return eri3KernelFun(ipair_ab, ish_c, spd_ab, shd_c, eri3_kernel);
                          };
}

lints::ERI2Kernel::ERI2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b)
{
    ecoeffs_bra_ = ecoeffsSHARK(sh_data_a, false);
    ecoeffs_ket_ = ecoeffsSHARK(sh_data_b, true);

    int la = sh_data_a.l_;
    int lb = sh_data_b.l_;
    int lab = la + lb;
    boys_grid_ = BoysGrid(lab);

    if (lab <= _max_l_rollout_)
        eri2_kernelfun_ = eri2_kernelfuns.at({la, lb});
    else
        eri2_kernelfun_ = [](const size_t ish_a, const size_t ish_b,
                             const ShellData &shd_a, const ShellData &shd_b,
                             const ERI2Kernel *eri2_kernel) -> vec2d
        {
            return eri2KernelFun(ish_a, ish_b, shd_a, shd_b, eri2_kernel);
        };
}

lints::ERI4D1Kernel::ERI4D1Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd)
{
    ecoeffs0_bra_ = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs1_bra_ = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket_ = ecoeffsSHARK(sp_data_cd, true);
    ecoeffs1_ket_ = ecoeffsD1SHARK(sp_data_cd, true);

    auto [la, lb] = sp_data_ab.getLPair();
    auto [lc, ld] = sp_data_cd.getLPair();
    int labcd = la + lb + lc + ld;
    boys_grid_ = BoysGrid(labcd + 1);

    if (labcd <= _max_l_rollout_)
        eri4d1_kernelfun_ = eri4d1_kernelfuns.at({la, lb, lc, ld});
    else
        eri4d1_kernelfun_ = [](const size_t ipair_ab, const size_t ipair_cd,
                               const ShellPairData &spd_ab, const ShellPairData &spd_cd,
                               const ERI4D1Kernel *eri4d1_kernel) -> std::array<vec4d, 12>
        {
            return eri4d1KernelFun(ipair_ab, ipair_cd, spd_ab, spd_cd, eri4d1_kernel);
        };
}

lints::ERI3D1Kernel::ERI3D1Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c)
{
    ecoeffs0_bra_ = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs1_bra_ = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket_ = ecoeffsSHARK(sh_data_c, true);

    auto [la, lb] = sp_data_ab.getLPair();
    int lc = sh_data_c.l_;
    int labc = la + lb + lc;
    boys_grid_ = BoysGrid(labc + 1);

    if (labc <= _max_l_rollout_)
        eri3d1_kernelfun_ = eri3d1_kernelfuns.at({la, lb, lc});
    else
        eri3d1_kernelfun_ = [](const size_t ipair_ab, const size_t ish_c,
                               const ShellPairData &spd_ab, const ShellData &shd_c,
                               const ERI3D1Kernel *eri3d1_kernel) -> std::array<vec3d, 9>
        {
            return eri3d1KernelFun(ipair_ab, ish_c, spd_ab, shd_c, eri3d1_kernel);
        };
}

lints::ERI2D1Kernel::ERI2D1Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b)
{
    ecoeffs_bra_ = ecoeffsSHARK(sh_data_a, false);
    ecoeffs_ket_ = ecoeffsSHARK(sh_data_b, true);

    int la = sh_data_a.l_;
    int lb = sh_data_b.l_;
    int lab = la + lb;
    boys_grid_ = BoysGrid(lab + 1);

    if (lab <= _max_l_rollout_)
        eri2d1_kernelfun_ = eri2d1_kernelfuns.at({la, lb});
    else
        eri2d1_kernelfun_ = [](const size_t ish_a, const size_t ish_b,
                               const ShellData &shd_a, const ShellData &shd_b,
                               const ERI2D1Kernel *eri2d1_kernel) -> std::array<vec2d, 6>
        {
            return eri2d1KernelFun(ish_a, ish_b, shd_a, shd_b, eri2d1_kernel);
        };
}

lints::ERI2D2Kernel::ERI2D2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b)
{
    ecoeffs_bra_ = ecoeffsSHARK(sh_data_a, false);
    ecoeffs_ket_ = ecoeffsSHARK(sh_data_b, true);

    int la = sh_data_a.l_;
    int lb = sh_data_b.l_;
    int lab = la + lb;
    boys_grid_ = BoysGrid(lab + 2);

    if (lab <= _max_l_rollout_)
        eri2d2_kernelfun_ = eri2d2_kernelfuns.at({la, lb});
    else
        eri2d2_kernelfun_ = [](const size_t ish_a, const size_t ish_b,
                               const ShellData &shd_a, const ShellData &shd_b,
                               const ERI2D2Kernel *eri2d2_kernel) -> arr2d<vec2d, 6, 6>
        {
            return eri2d2KernelFun(ish_a, ish_b, shd_a, shd_b, eri2d2_kernel);
        };
}

lints::ERI4SOCKernel::ERI4SOCKernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd)
{
    ecoeffs1_bra_ = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket_ = ecoeffsSHARK(sp_data_cd, true);

    auto [la, lb] = sp_data_ab.getLPair();
    auto [lc, ld] = sp_data_cd.getLPair();
    int labcd = la + lb + lc + ld;
    boys_grid_ = BoysGrid(labcd + 1);

    if (labcd <= _max_l_rollout_)
        eri4soc_kernelfun_ = eri4soc_kernelfuns.at({la, lb, lc, ld});
    else
        eri4soc_kernelfun_ = [](const size_t ipair_ab, const size_t ipair_cd,
                                const ShellPairData &spd_ab, const ShellPairData &spd_cd,
                                const ERI4SOCKernel *eri4soc_kernel) -> std::array<vec4d, 3>
        {
            return eri4socKernelFun(ipair_ab, ipair_cd, spd_ab, spd_cd, eri4soc_kernel);
        };
}

lints::ERI3SOCKernel::ERI3SOCKernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c)
{
    ecoeffs1_bra_ = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket_ = ecoeffsSHARK(sh_data_c, true);

    auto [la, lb] = sp_data_ab.getLPair();
    int lc = sh_data_c.l_;
    int labc = la + lb + lc;
    boys_grid_ = BoysGrid(labc + 1);

    if (labc <= _max_l_rollout_)
        eri3soc_kernelfun_ = eri3soc_kernelfuns.at({la, lb, lc});
    else
        eri3soc_kernelfun_ = [](const size_t ipair_ab, const size_t ish_c,
                                const ShellPairData &spd_ab, const ShellData &shd_c,
                                const ERI3SOCKernel *eri3soc_kernel) -> std::array<vec3d, 3>
        {
            return eri3socKernelFun(ipair_ab, ish_c, spd_ab, shd_c, eri3soc_kernel);
        };
}