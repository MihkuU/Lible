#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/defs.hpp>
#include <lible/ints/utils.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <format>
#include <tuple>

namespace LI = lible::ints;

using std::array, std::string, std::vector;

namespace lible::ints
{
    // regular ERI kernels

    template <int la, int lb, int lc, int ld>
    vec4d eri4KernelFun(const int ipair_ab, const int ipair_cd,
                        const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                        const ERI4Kernel *eri4_kernel);

    vec4d eri4KernelFun(const int ipair_ab, const int ipair_cd,
                        const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                        const ERI4Kernel *eri4_kernel);

    template <int la, int lb, int lc>
    vec3d eri3KernelFun(const int ipair_ab, const int ishell_c,
                        const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                        const ERI3Kernel *eri3_kernel);

    vec3d eri3KernelFun(const int ipair_ab, const int ishell_c,
                        const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                        const ERI3Kernel *eri3_kernel);

    template <int la, int lb>
    vec2d eri2KernelFun(const int ishell_a, const int ishell_b,
                        const ShellData &sh_data_a, const ShellData &sh_data_b,
                        const ERI2Kernel *eri2_kernel);

    vec2d eri2KernelFun(const int ishell_a, const int ishell_b,
                        const ShellData &sh_data_a, const ShellData &sh_data_b,
                        const ERI2Kernel *eri2_kernel);

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
        {{0, 6, 0, 0}, eri4KernelFun<0, 6, 0, 0>}};

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
        {{0, 6, 0}, eri3KernelFun<0, 6, 0>}};

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
        {{6, 0}, eri2KernelFun<6, 0>}};
}

namespace lible::ints
{
    // Derivative ERI kernels

    // 2-center
    template <int la, int lb>
    std::array<vec2d, 6> eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                         const ShellData &sh_data_a, const ShellData &sh_data_b,
                                         const ERI2D1Kernel *eri2d1_kernel);

    std::array<vec2d, 6> eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                         const ShellData &sh_data_a, const ShellData &sh_data_b,
                                         const ERI2D1Kernel *eri2d1_kernel);

    template <int la, int lb>
    arr2d<vec2d, 6, 6> eri2d2KernelFun(const int ishell_a, const int ishell_b,
                                       const ShellData &sh_data_a,
                                       const ShellData &sh_data_b,
                                       const ERI2D2Kernel *eri2d2_kernel);

    arr2d<vec2d, 6, 6> eri2d2KernelFun(const int ishell_a, const int ishell_b,
                                       const ShellData &sh_data_a,
                                       const ShellData &sh_data_b,
                                       const ERI2D2Kernel *eri2d2_kernel);

    // 3-center
    template <int la, int lb, int lc>
    std::array<vec3d, 9> eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                         const ShellPairData &sp_data_ab,
                                         const ShellData &sh_data_c,
                                         const ERI3D1Kernel *eri3d1_kernel);

    std::array<vec3d, 9> eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                         const ShellPairData &sp_data_ab,
                                         const ShellData &sh_data_c,
                                         const ERI3D1Kernel *eri3d1_kernel);

    // template <int la, int lb, int lc>
    // arr2d<vec3d, 9, 9> eri3d2KernelFun(const int ipair_ab, const int ishell_c,
    //                                    const ShellPairData &sp_data_ab,
    //                                    const ShellData &sh_data_c,
    //                                    const ERI3D2Kernel *eri3d2_kernel);

    // arr2d<vec3d, 9, 9> eri3d2KernelFun(const int ipair_ab, const int ishell_c,
    //                                    const ShellPairData &sp_data_ab,
    //                                    const ShellData &sh_data_c,
    //                                    const ERI3D2Kernel *eri3d2_kernel);

    template <int la, int lb, int lc>
    std::array<vec3d, 3> eri3socKernelFun(const int ipair_ab, const int ishell_c,
                                          const ShellPairData &sp_data_ab,
                                          const ShellData &sh_data_c,
                                          const ERI3SOCKernel *eri3soc_kernel);

    std::array<vec3d, 3> eri3socKernelFun(const int ipair_ab, const int ishell_c,
                                          const ShellPairData &sp_data_ab,
                                          const ShellData &sh_data_c,
                                          const ERI3SOCKernel *eri3soc_kernel);

    // 4-center
    template <int la, int lb, int lc, int ld>
    std::array<vec4d, 12> eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4D1Kernel *eri4d1_kernel);

    std::array<vec4d, 12> eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4D1Kernel *eri4d1_kernel);

    // template <int la, int lb, int lc, int ld>
    // arr2d<vec4d, 12, 12> eri4d2KernelFun(const int ipair_ab, const int ipair_cd,
    //                                      const ShellPairData &sp_data_ab,
    //                                      const ShellPairData &sp_data_cd,
    //                                      const ERI4D2Kernel *eri4d2_kernel);

    // arr2d<vec4d, 12, 12> eri4d2KernelFun(const int ipair_ab, const int ipair_cd,
    //                                      const ShellPairData &sp_data_ab,
    //                                      const ShellPairData &sp_data_cd,
    //                                      const ERI4D2Kernel *eri4d2_kernel);

    template <int la, int lb, int lc, int ld>
    std::array<vec4d, 3> eri4socKernelFun(const int ipair_ab, const int ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4SOCKernel *eri4soc_kernel);

    std::array<vec4d, 3> eri4socKernelFun(const int ipair_ab, const int ipair_cd,
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
        {{6, 0}, eri2d1KernelFun<6, 0>}};

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
        {{6, 0}, eri2d2KernelFun<6, 0>}};

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
        {{0, 6, 0}, eri3d1KernelFun<0, 6, 0>}};

    // const std::map<std::tuple<int, int, int>, eri3d2_kernelfun_t> eri3d2_kernelfuns{
    //     {{0, 0, 0}, eri3d2KernelFun<0, 0, 0>},
    //     {{0, 0, 1}, eri3d2KernelFun<0, 0, 1>},
    //     {{0, 0, 2}, eri3d2KernelFun<0, 0, 2>},
    //     {{0, 0, 3}, eri3d2KernelFun<0, 0, 3>},
    //     {{0, 0, 4}, eri3d2KernelFun<0, 0, 4>},
    //     {{0, 0, 5}, eri3d2KernelFun<0, 0, 5>},
    //     {{0, 0, 6}, eri3d2KernelFun<0, 0, 6>},
    //     {{1, 0, 0}, eri3d2KernelFun<1, 0, 0>},
    //     {{0, 1, 0}, eri3d2KernelFun<0, 1, 0>},
    //     {{1, 0, 1}, eri3d2KernelFun<1, 0, 1>},
    //     {{0, 1, 1}, eri3d2KernelFun<0, 1, 1>},
    //     {{1, 0, 2}, eri3d2KernelFun<1, 0, 2>},
    //     {{0, 1, 2}, eri3d2KernelFun<0, 1, 2>},
    //     {{1, 0, 3}, eri3d2KernelFun<1, 0, 3>},
    //     {{0, 1, 3}, eri3d2KernelFun<0, 1, 3>},
    //     {{1, 0, 4}, eri3d2KernelFun<1, 0, 4>},
    //     {{0, 1, 4}, eri3d2KernelFun<0, 1, 4>},
    //     {{1, 0, 5}, eri3d2KernelFun<1, 0, 5>},
    //     {{0, 1, 5}, eri3d2KernelFun<0, 1, 5>},
    //     {{1, 1, 0}, eri3d2KernelFun<1, 1, 0>},
    //     {{1, 1, 1}, eri3d2KernelFun<1, 1, 1>},
    //     {{1, 1, 2}, eri3d2KernelFun<1, 1, 2>},
    //     {{1, 1, 3}, eri3d2KernelFun<1, 1, 3>},
    //     {{1, 1, 4}, eri3d2KernelFun<1, 1, 4>},
    //     {{2, 0, 0}, eri3d2KernelFun<2, 0, 0>},
    //     {{0, 2, 0}, eri3d2KernelFun<0, 2, 0>},
    //     {{2, 0, 1}, eri3d2KernelFun<2, 0, 1>},
    //     {{0, 2, 1}, eri3d2KernelFun<0, 2, 1>},
    //     {{2, 0, 2}, eri3d2KernelFun<2, 0, 2>},
    //     {{0, 2, 2}, eri3d2KernelFun<0, 2, 2>},
    //     {{2, 0, 3}, eri3d2KernelFun<2, 0, 3>},
    //     {{0, 2, 3}, eri3d2KernelFun<0, 2, 3>},
    //     {{2, 0, 4}, eri3d2KernelFun<2, 0, 4>},
    //     {{0, 2, 4}, eri3d2KernelFun<0, 2, 4>},
    //     {{2, 1, 0}, eri3d2KernelFun<2, 1, 0>},
    //     {{1, 2, 0}, eri3d2KernelFun<1, 2, 0>},
    //     {{2, 1, 1}, eri3d2KernelFun<2, 1, 1>},
    //     {{1, 2, 1}, eri3d2KernelFun<1, 2, 1>},
    //     {{2, 1, 2}, eri3d2KernelFun<2, 1, 2>},
    //     {{1, 2, 2}, eri3d2KernelFun<1, 2, 2>},
    //     {{2, 1, 3}, eri3d2KernelFun<2, 1, 3>},
    //     {{1, 2, 3}, eri3d2KernelFun<1, 2, 3>},
    //     {{2, 2, 0}, eri3d2KernelFun<2, 2, 0>},
    //     {{2, 2, 1}, eri3d2KernelFun<2, 2, 1>},
    //     {{2, 2, 2}, eri3d2KernelFun<2, 2, 2>},
    //     {{3, 0, 0}, eri3d2KernelFun<3, 0, 0>},
    //     {{0, 3, 0}, eri3d2KernelFun<0, 3, 0>},
    //     {{3, 0, 1}, eri3d2KernelFun<3, 0, 1>},
    //     {{0, 3, 1}, eri3d2KernelFun<0, 3, 1>},
    //     {{3, 0, 2}, eri3d2KernelFun<3, 0, 2>},
    //     {{0, 3, 2}, eri3d2KernelFun<0, 3, 2>},
    //     {{3, 0, 3}, eri3d2KernelFun<3, 0, 3>},
    //     {{0, 3, 3}, eri3d2KernelFun<0, 3, 3>},
    //     {{3, 1, 0}, eri3d2KernelFun<3, 1, 0>},
    //     {{1, 3, 0}, eri3d2KernelFun<1, 3, 0>},
    //     {{3, 1, 1}, eri3d2KernelFun<3, 1, 1>},
    //     {{1, 3, 1}, eri3d2KernelFun<1, 3, 1>},
    //     {{3, 1, 2}, eri3d2KernelFun<3, 1, 2>},
    //     {{1, 3, 2}, eri3d2KernelFun<1, 3, 2>},
    //     {{3, 2, 0}, eri3d2KernelFun<3, 2, 0>},
    //     {{2, 3, 0}, eri3d2KernelFun<2, 3, 0>},
    //     {{3, 2, 1}, eri3d2KernelFun<3, 2, 1>},
    //     {{2, 3, 1}, eri3d2KernelFun<2, 3, 1>},
    //     {{3, 3, 0}, eri3d2KernelFun<3, 3, 0>},
    //     {{4, 0, 0}, eri3d2KernelFun<4, 0, 0>},
    //     {{0, 4, 0}, eri3d2KernelFun<0, 4, 0>},
    //     {{4, 0, 1}, eri3d2KernelFun<4, 0, 1>},
    //     {{0, 4, 1}, eri3d2KernelFun<0, 4, 1>},
    //     {{4, 0, 2}, eri3d2KernelFun<4, 0, 2>},
    //     {{0, 4, 2}, eri3d2KernelFun<0, 4, 2>},
    //     {{4, 1, 0}, eri3d2KernelFun<4, 1, 0>},
    //     {{1, 4, 0}, eri3d2KernelFun<1, 4, 0>},
    //     {{4, 1, 1}, eri3d2KernelFun<4, 1, 1>},
    //     {{1, 4, 1}, eri3d2KernelFun<1, 4, 1>},
    //     {{4, 2, 0}, eri3d2KernelFun<4, 2, 0>},
    //     {{2, 4, 0}, eri3d2KernelFun<2, 4, 0>},
    //     {{5, 0, 0}, eri3d2KernelFun<5, 0, 0>},
    //     {{0, 5, 0}, eri3d2KernelFun<0, 5, 0>},
    //     {{5, 0, 1}, eri3d2KernelFun<5, 0, 1>},
    //     {{0, 5, 1}, eri3d2KernelFun<0, 5, 1>},
    //     {{5, 1, 0}, eri3d2KernelFun<5, 1, 0>},
    //     {{1, 5, 0}, eri3d2KernelFun<1, 5, 0>},
    //     {{6, 0, 0}, eri3d2KernelFun<6, 0, 0>},
    //     {{0, 6, 0}, eri3d2KernelFun<0, 6, 0>}};

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
        {{0, 6, 0, 0}, eri4d1KernelFun<0, 6, 0, 0>}};

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
        {{0, 6, 0, 0}, eri4socKernelFun<0, 6, 0, 0>}};

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
        {{0, 6, 0}, eri3socKernelFun<0, 6, 0>}};
}

LI::ERI4Kernel::ERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                           const eri4_kernelfun_t &eri4_kernelfun)
    : eri4_kernelfun(eri4_kernelfun)
{

    ecoeffs_bra = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs_ket = ecoeffsSHARK(sp_data_cd, true);

    int labcd = sp_data_ab.la + sp_data_ab.lb + sp_data_cd.la + sp_data_cd.lb;
    boys_grid = BoysGrid(labcd);
}

LI::ERI3Kernel::ERI3Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                           const eri3_kernelfun_t &eri3_kernelfun)
    : eri3_kernelfun(eri3_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs_ket = ecoeffsSHARK(sh_data_c, true);

    int labc = sp_data_ab.la + sp_data_ab.lb + sh_data_c.l;
    boys_grid = BoysGrid(labc);
}

LI::ERI2Kernel::ERI2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b,
                           const eri2_kernelfun_t &eri2_kernelfun)
    : eri2_kernelfun(eri2_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARK(sh_data_a, false);
    ecoeffs_ket = ecoeffsSHARK(sh_data_b, true);

    int lab = sh_data_a.l + sh_data_b.l;
    boys_grid = BoysGrid(lab);
}

LI::ERI4D1Kernel::ERI4D1Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                               const eri4d1_kernelfun_t &eri4d1_kernelfun)
    : eri4d1_kernelfun(eri4d1_kernelfun)
{
    ecoeffs0_bra = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs1_bra = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket = ecoeffsSHARK(sp_data_cd, true);
    ecoeffs1_ket = ecoeffsD1SHARK(sp_data_cd, true);

    int labcd = sp_data_ab.la + sp_data_ab.lb + sp_data_cd.la + sp_data_cd.lb;
    boys_grid = BoysGrid(labcd + 1);
}

LI::ERI3D1Kernel::ERI3D1Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                               const eri3d1_kernelfun_t &eri3d1_kernelfun)
    : eri3d1_kernelfun(eri3d1_kernelfun)
{
    ecoeffs0_bra = ecoeffsSHARK(sp_data_ab, false);
    ecoeffs1_bra = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket = ecoeffsSHARK(sh_data_c, true);

    int labc = sp_data_ab.la + sp_data_ab.lb + sh_data_c.l;
    boys_grid = BoysGrid(labc + 1);
}

LI::ERI2D1Kernel::ERI2D1Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b,
                               const eri2d1_kernelfun_t &eri2d1_kernelfun)
    : eri2d1_kernelfun(eri2d1_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARK(sh_data_a, false);
    ecoeffs_ket = ecoeffsSHARK(sh_data_b, true);

    int lab = sh_data_a.l + sh_data_b.l;
    boys_grid = BoysGrid(lab + 1);
}

LI::ERI2D2Kernel::ERI2D2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b,
                               const eri2d2_kernelfun_t &eri2d2_kernelfun)
    : eri2d2_kernelfun(eri2d2_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARK(sh_data_a, false);
    ecoeffs_ket = ecoeffsSHARK(sh_data_b, true);

    int lab = sh_data_a.l + sh_data_b.l;
    boys_grid = BoysGrid(lab + 2);
}

LI::ERI4SOCKernel::ERI4SOCKernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                                 const eri4soc_kernelfun_t &eri4soc_kernelfun)
    : eri4soc_kernelfun(eri4soc_kernelfun)
{
    ecoeffs1_bra = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket = ecoeffsSHARK(sp_data_cd, true);

    int labcd = sp_data_ab.la + sp_data_ab.lb + sp_data_cd.la + sp_data_cd.lb;
    boys_grid = BoysGrid(labcd + 1);
}

LI::ERI3SOCKernel::ERI3SOCKernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                                 const eri3soc_kernelfun_t &eri3soc_kernelfun)
    : eri3soc_kernelfun(eri3soc_kernelfun)
{
    ecoeffs1_bra = ecoeffsD1SHARK(sp_data_ab, false);
    ecoeffs0_ket = ecoeffsSHARK(sh_data_c, true);

    int labc = sp_data_ab.la + sp_data_ab.lb + sh_data_c.l;
    boys_grid = BoysGrid(labc + 1);
}

LI::ERI4Kernel LI::deployERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sp_data_cd.la;
    int ld = sp_data_cd.lb;

    int labcd = la + lb + lc + ld;
    if (labcd <= _max_l_rollout_)
        return ERI4Kernel(sp_data_ab, sp_data_cd, eri4_kernelfuns.at({la, lb, lc, ld}));
    else
        return ERI4Kernel(sp_data_ab, sp_data_cd,
                          [](const int ipair_ab, const int ipair_cd,
                             const ShellPairData &sp_data_ab,
                             const ShellPairData &sp_data_cd,
                             const ERI4Kernel *eri4_kernel) -> vec4d
                          {
                              return eri4KernelFun(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                                   eri4_kernel);
                          });
}

LI::ERI3Kernel LI::deployERI3Kernel(const ShellPairData &sp_data_ab,
                                    const ShellData &sh_data_c)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sh_data_c.l;

    int labc = la + lb + lc;
    if (labc <= _max_l_rollout_)
        return ERI3Kernel(sp_data_ab, sh_data_c, eri3_kernelfuns.at({la, lb, lc}));
    else
        return ERI3Kernel(sp_data_ab, sh_data_c,
                          [](const int ipair_ab, const int ish_c,
                             const ShellPairData &sp_data_ab,
                             const ShellData &sh_data_c,
                             const ERI3Kernel *eri3_kernel) -> vec3d
                          {
                              return eri3KernelFun(ipair_ab, ish_c, sp_data_ab, sh_data_c,
                                                   eri3_kernel);
                          });
}

LI::ERI2Kernel LI::deployERI2Kernel(const ShellData &sh_data_a,
                                    const ShellData &sh_data_b)
{
    int la = sh_data_a.l;
    int lb = sh_data_b.l;

    int lab = la + lb;
    if (lab <= _max_l_rollout_)
        return ERI2Kernel(sh_data_a, sh_data_b, eri2_kernelfuns.at({la, lb}));
    else
        return ERI2Kernel(sh_data_a, sh_data_b,
                          [](const int ish_a, const int ish_b,
                             const ShellData &sh_data_a,
                             const ShellData &sh_data_b,
                             const ERI2Kernel *eri2_kernel) -> vec2d
                          {
                              return eri2KernelFun(ish_a, ish_b, sh_data_a, sh_data_b,
                                                   eri2_kernel);
                          });
}

LI::ERI4D1Kernel LI::deployERI4D1Kernel(const ShellPairData &sp_data_ab,
                                        const ShellPairData &sp_data_cd)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sp_data_cd.la;
    int ld = sp_data_cd.lb;

    int labcd = la + lb + lc + ld;
    if (labcd <= _max_l_rollout_)
        return ERI4D1Kernel(sp_data_ab, sp_data_cd, eri4d1_kernelfuns.at({la, lb, lc, ld}));
    else
        return ERI4D1Kernel(sp_data_ab, sp_data_cd,
                            [](const int ipair_ab, const int ipair_cd,
                               const ShellPairData &sp_data_ab,
                               const ShellPairData &sp_data_cd,
                               const ERI4D1Kernel *eri4d1_kernel) -> std::array<vec4d, 12>
                            {
                                return eri4d1KernelFun(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                                       eri4d1_kernel);
                            });
}

LI::ERI3D1Kernel LI::deployERI3D1Kernel(const ShellPairData &sp_data_ab,
                                        const ShellData &sh_data_c)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sh_data_c.l;

    int labc = la + lb + lc;
    if (labc <= _max_l_rollout_)
        return ERI3D1Kernel(sp_data_ab, sh_data_c, eri3d1_kernelfuns.at({la, lb, lc}));
    else
        return ERI3D1Kernel(sp_data_ab, sh_data_c,
                            [](const int ipair_ab, const int ish_c,
                               const ShellPairData &sp_data_ab,
                               const ShellData &sh_data_c,
                               const ERI3D1Kernel *eri3d1_kernel) -> std::array<vec3d, 9>
                            {
                                return eri3d1KernelFun(ipair_ab, ish_c, sp_data_ab, sh_data_c,
                                                       eri3d1_kernel);
                            });
}

LI::ERI2D1Kernel LI::deployERI2D1Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b)
{
    int la = sh_data_a.l;
    int lb = sh_data_b.l;

    int lab = la + lb;
    if (lab <= _max_l_rollout_)
        return ERI2D1Kernel(sh_data_a, sh_data_b, eri2d1_kernelfuns.at({la, lb}));
    else
        return ERI2D1Kernel(sh_data_a, sh_data_b,
                            [](const int ish_a, const int ish_b,
                               const ShellData &sh_data_a,
                               const ShellData &sh_data_b,
                               const ERI2D1Kernel *eri2d1_kernel) -> std::array<vec2d, 6>
                            {
                                return eri2d1KernelFun(ish_a, ish_b, sh_data_a, sh_data_b,
                                                       eri2d1_kernel);
                            });
}

LI::ERI2D2Kernel LI::deployERI2D2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b)
{
    int la = sh_data_a.l;
    int lb = sh_data_b.l;

    int lab = la + lb;
    if (lab <= _max_l_rollout_)
        return ERI2D2Kernel(sh_data_a, sh_data_b, eri2d2_kernelfuns.at({la, lb}));
    else
        return ERI2D2Kernel(sh_data_a, sh_data_b,
                            [](const int ish_a, const int ish_b,
                               const ShellData &sh_data_a,
                               const ShellData &sh_data_b,
                               const ERI2D2Kernel *eri2d2_kernel) -> arr2d<vec2d, 6, 6>
                            {
                                return eri2d2KernelFun(ish_a, ish_b, sh_data_a, sh_data_b,
                                                       eri2d2_kernel);
                            });
}

LI::ERI4SOCKernel LI::deployERI4SOCKernel(const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sp_data_cd.la;
    int ld = sp_data_cd.lb;

    int labcd = la + lb + lc + ld;
    if (labcd <= _max_l_rollout_)
        return ERI4SOCKernel(sp_data_ab, sp_data_cd, eri4soc_kernelfuns.at({la, lb, lc, ld}));
    else
        return ERI4SOCKernel(sp_data_ab, sp_data_cd,
                             [](const int ipair_ab, const int ipair_cd,
                                const ShellPairData &sp_data_ab,
                                const ShellPairData &sp_data_cd,
                                const ERI4SOCKernel *eri4soc_kernel) -> std::array<vec4d, 3>
                             {
                                 return eri4socKernelFun(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd,
                                                         eri4soc_kernel);
                             });
}

LI::ERI3SOCKernel LI::deployERI3SOCKernel(const ShellPairData &sp_data_ab,
                                          const ShellData &sh_data_c)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sh_data_c.l;

    int labc = la + lb + lc;
    if (labc <= _max_l_rollout_)
        return ERI3SOCKernel(sp_data_ab, sh_data_c, eri3soc_kernelfuns.at({la, lb, lc}));
    else
        return ERI3SOCKernel(sp_data_ab, sh_data_c,
                             [](const int ipair_ab, const int ish_c,
                                const ShellPairData &sp_data_ab,
                                const ShellData &sh_data_c,
                                const ERI3SOCKernel *eri3soc_kernel) -> std::array<vec3d, 3>
                             {
                                 return eri3socKernelFun(ipair_ab, ish_c, sp_data_ab, sh_data_c,
                                                         eri3soc_kernel);
                             });
}