#include <lible/ints/defs.hpp>
#include <lible/ints/ecoeffs.hpp>
#include <lible/ints/twoel/eri_kernels.hpp>

#include <format>
#include <tuple>

namespace LI = lible::ints;

using std::string;

namespace lible::ints
{
    template <int la, int lb>
    std::array<vec2d, 6> eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                         const ShellData &sh_data_a, const ShellData &sh_data_b,
                                         const ERI2D1Kernel *eri2d1_kernel);

    std::array<vec2d, 6> eri2d1KernelFun(const int ishell_a, const int ishell_b,
                                         const ShellData &sh_data_a, const ShellData &sh_data_b,
                                         const ERI2D1Kernel *eri2d1_kernel);

    template <int la, int lb, int lc>
    std::array<vec3d, 9> eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                         const ShellPairData &sp_data_ab,
                                         const ShellData &sh_data_c,
                                         const ERI3D1Kernel *eri3d1_kernel);

    std::array<vec3d, 9> eri3d1KernelFun(const int ipair_ab, const int ishell_c,
                                         const ShellPairData &sp_data_ab,
                                         const ShellData &sh_data_c,
                                         const ERI3D1Kernel *eri3d1_kernel);

    template <int la, int lb, int lc, int ld>
    std::array<vec4d, 12> eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4D1Kernel *eri4d1_kernel);

    std::array<vec4d, 12> eri4d1KernelFun(const int ipair_ab, const int ipair_cd,
                                          const ShellPairData &sp_data_ab,
                                          const ShellPairData &sp_data_cd,
                                          const ERI4D1Kernel *eri4d1_kernel);

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

    const std::map<std::tuple<int, int, int>, eri3d1_kernelfun_t> eri3d1_kernelfuns{
        {{0, 0, 0}, eri3d1KernelFun<0, 0, 0>},
        {{0, 0, 0}, eri3d1KernelFun<0, 0, 0>},
        {{0, 0, 1}, eri3d1KernelFun<0, 0, 1>},
        {{0, 0, 1}, eri3d1KernelFun<0, 0, 1>},
        {{0, 0, 2}, eri3d1KernelFun<0, 0, 2>},
        {{0, 0, 2}, eri3d1KernelFun<0, 0, 2>},
        {{0, 0, 3}, eri3d1KernelFun<0, 0, 3>},
        {{0, 0, 3}, eri3d1KernelFun<0, 0, 3>},
        {{0, 0, 4}, eri3d1KernelFun<0, 0, 4>},
        {{0, 0, 4}, eri3d1KernelFun<0, 0, 4>},
        {{0, 0, 5}, eri3d1KernelFun<0, 0, 5>},
        {{0, 0, 5}, eri3d1KernelFun<0, 0, 5>},
        {{0, 0, 6}, eri3d1KernelFun<0, 0, 6>},
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
        {{1, 1, 0}, eri3d1KernelFun<1, 1, 0>},
        {{1, 1, 1}, eri3d1KernelFun<1, 1, 1>},
        {{1, 1, 1}, eri3d1KernelFun<1, 1, 1>},
        {{1, 1, 2}, eri3d1KernelFun<1, 1, 2>},
        {{1, 1, 2}, eri3d1KernelFun<1, 1, 2>},
        {{1, 1, 3}, eri3d1KernelFun<1, 1, 3>},
        {{1, 1, 3}, eri3d1KernelFun<1, 1, 3>},
        {{1, 1, 4}, eri3d1KernelFun<1, 1, 4>},
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
        {{2, 2, 0}, eri3d1KernelFun<2, 2, 0>},
        {{2, 2, 1}, eri3d1KernelFun<2, 2, 1>},
        {{2, 2, 1}, eri3d1KernelFun<2, 2, 1>},
        {{2, 2, 2}, eri3d1KernelFun<2, 2, 2>},
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
}

LI::ERI4D1Kernel::ERI4D1Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                               const eri4d1_kernelfun_t &eri4d1_kernelfun)
    : eri4d1_kernelfun(eri4d1_kernelfun)
{
    ecoeffs0_bra = ecoeffsSHARKBig(sp_data_ab);
    ecoeffs1_bra = ecoeffsD1SHARK(sp_data_ab);
    ecoeffs0_ket = ecoeffsSHARKBig(sp_data_cd);
    ecoeffs1_ket = ecoeffsD1SHARK(sp_data_cd);

    int labcd = sp_data_ab.la + sp_data_ab.lb + sp_data_cd.la + sp_data_cd.lb;
    boys_grid = BoysGrid(labcd + 1);
}

LI::ERI3D1Kernel::ERI3D1Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                               const eri3d1_kernelfun_t &eri3d1_kernelfun)
    : eri3d1_kernelfun(eri3d1_kernelfun)
{
    ecoeffs0_bra = ecoeffsSHARKBig(sp_data_ab);
    ecoeffs1_bra = ecoeffsD1SHARK(sp_data_ab);
    ecoeffs0_ket = ecoeffsSHARKBig(sh_data_c);

    int labc = sp_data_ab.la + sp_data_ab.lb + sh_data_c.l;
    boys_grid = BoysGrid(labc + 1);
}

LI::ERI2D1Kernel::ERI2D1Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b,
                               const eri2d1_kernelfun_t &eri2d1_kernelfun)
    : eri2d1_kernelfun(eri2d1_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARKBig(sh_data_a);
    ecoeffs_ket = ecoeffsSHARKBig(sh_data_b);

    int lab = sh_data_a.l + sh_data_b.l;
    boys_grid = BoysGrid(lab + 1);
}

LI::ERI4D1Kernel LI::deployERI4D1Kernel(const ShellPairData &sp_data_ab,
                                        const ShellPairData &sp_data_cd)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sp_data_cd.la;
    int ld = sp_data_cd.lb;

    int labcd = la + lb + lc + ld;
    int l_max = _eri_kernel_max_l_;
    if (labcd > _eri_kernel_max_l_)
    {
        string msg = std::format("deployERI4D1Kernel(): lab + lcd = {} is larger than the allowed max: {}!\n",
                                 labcd, l_max);
        throw std::runtime_error(msg);
    }

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
    int l_max = _eri_kernel_max_l_;
    if (labc > _eri_kernel_max_l_)
    {
        string msg = std::format("deployERI3D1Kernel(): lab + lc = {} is larger than the allowed max: {}!\n",
                                 labc, l_max);
        throw std::runtime_error(msg);
    }

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

LI::ERI2D1Kernel LI::deployERI2D1Kernel(const ShellData &sh_data_a,
                                        const ShellData &sh_data_b)
{
    int la = sh_data_a.l;
    int lb = sh_data_b.l;

    int lab = la + lb;
    int l_max = _eri_kernel_max_l_;
    if (lab > _eri_kernel_max_l_)
    {
        string msg = std::format("deployERI2D1Kernel(): la + lb = {} is larger than the allowed max: {}!\n",
                                 lab, l_max);
        throw std::runtime_error(msg);
    }

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
