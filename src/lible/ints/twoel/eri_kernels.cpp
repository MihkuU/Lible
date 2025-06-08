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
        {{0, 0, 1, 1}, eri4KernelFun<0, 0, 1, 1>},
        {{0, 0, 2, 0}, eri4KernelFun<0, 0, 2, 0>},
        {{0, 0, 2, 1}, eri4KernelFun<0, 0, 2, 1>},
        {{0, 0, 2, 2}, eri4KernelFun<0, 0, 2, 2>},
        {{0, 0, 3, 0}, eri4KernelFun<0, 0, 3, 0>},
        {{0, 0, 3, 1}, eri4KernelFun<0, 0, 3, 1>},
        {{0, 0, 3, 2}, eri4KernelFun<0, 0, 3, 2>},
        {{0, 0, 3, 3}, eri4KernelFun<0, 0, 3, 3>},
        {{0, 0, 4, 0}, eri4KernelFun<0, 0, 4, 0>},
        {{0, 0, 4, 1}, eri4KernelFun<0, 0, 4, 1>},
        {{0, 0, 4, 2}, eri4KernelFun<0, 0, 4, 2>},
        {{0, 0, 5, 0}, eri4KernelFun<0, 0, 5, 0>},
        {{0, 0, 5, 1}, eri4KernelFun<0, 0, 5, 1>},
        {{0, 0, 6, 0}, eri4KernelFun<0, 0, 6, 0>},
        {{1, 0, 0, 0}, eri4KernelFun<1, 0, 0, 0>},
        {{1, 0, 1, 0}, eri4KernelFun<1, 0, 1, 0>},
        {{1, 0, 1, 1}, eri4KernelFun<1, 0, 1, 1>},
        {{1, 0, 2, 0}, eri4KernelFun<1, 0, 2, 0>},
        {{1, 0, 2, 1}, eri4KernelFun<1, 0, 2, 1>},
        {{1, 0, 2, 2}, eri4KernelFun<1, 0, 2, 2>},
        {{1, 0, 3, 0}, eri4KernelFun<1, 0, 3, 0>},
        {{1, 0, 3, 1}, eri4KernelFun<1, 0, 3, 1>},
        {{1, 0, 3, 2}, eri4KernelFun<1, 0, 3, 2>},
        {{1, 0, 4, 0}, eri4KernelFun<1, 0, 4, 0>},
        {{1, 0, 4, 1}, eri4KernelFun<1, 0, 4, 1>},
        {{1, 0, 5, 0}, eri4KernelFun<1, 0, 5, 0>},
        {{1, 1, 0, 0}, eri4KernelFun<1, 1, 0, 0>},
        {{1, 1, 1, 0}, eri4KernelFun<1, 1, 1, 0>},
        {{1, 1, 1, 1}, eri4KernelFun<1, 1, 1, 1>},
        {{1, 1, 2, 0}, eri4KernelFun<1, 1, 2, 0>},
        {{1, 1, 2, 1}, eri4KernelFun<1, 1, 2, 1>},
        {{1, 1, 2, 2}, eri4KernelFun<1, 1, 2, 2>},
        {{1, 1, 3, 0}, eri4KernelFun<1, 1, 3, 0>},
        {{1, 1, 3, 1}, eri4KernelFun<1, 1, 3, 1>},
        {{1, 1, 4, 0}, eri4KernelFun<1, 1, 4, 0>},
        {{2, 0, 0, 0}, eri4KernelFun<2, 0, 0, 0>},
        {{2, 0, 1, 0}, eri4KernelFun<2, 0, 1, 0>},
        {{2, 0, 1, 1}, eri4KernelFun<2, 0, 1, 1>},
        {{2, 0, 2, 0}, eri4KernelFun<2, 0, 2, 0>},
        {{2, 0, 2, 1}, eri4KernelFun<2, 0, 2, 1>},
        {{2, 0, 2, 2}, eri4KernelFun<2, 0, 2, 2>},
        {{2, 0, 3, 0}, eri4KernelFun<2, 0, 3, 0>},
        {{2, 0, 3, 1}, eri4KernelFun<2, 0, 3, 1>},
        {{2, 0, 4, 0}, eri4KernelFun<2, 0, 4, 0>},
        {{2, 1, 0, 0}, eri4KernelFun<2, 1, 0, 0>},
        {{2, 1, 1, 0}, eri4KernelFun<2, 1, 1, 0>},
        {{2, 1, 1, 1}, eri4KernelFun<2, 1, 1, 1>},
        {{2, 1, 2, 0}, eri4KernelFun<2, 1, 2, 0>},
        {{2, 1, 2, 1}, eri4KernelFun<2, 1, 2, 1>},
        {{2, 1, 3, 0}, eri4KernelFun<2, 1, 3, 0>},
        {{2, 2, 0, 0}, eri4KernelFun<2, 2, 0, 0>},
        {{2, 2, 1, 0}, eri4KernelFun<2, 2, 1, 0>},
        {{2, 2, 1, 1}, eri4KernelFun<2, 2, 1, 1>},
        {{2, 2, 2, 0}, eri4KernelFun<2, 2, 2, 0>},
        {{3, 0, 0, 0}, eri4KernelFun<3, 0, 0, 0>},
        {{3, 0, 1, 0}, eri4KernelFun<3, 0, 1, 0>},
        {{3, 0, 1, 1}, eri4KernelFun<3, 0, 1, 1>},
        {{3, 0, 2, 0}, eri4KernelFun<3, 0, 2, 0>},
        {{3, 0, 2, 1}, eri4KernelFun<3, 0, 2, 1>},
        {{3, 0, 3, 0}, eri4KernelFun<3, 0, 3, 0>},
        {{3, 1, 0, 0}, eri4KernelFun<3, 1, 0, 0>},
        {{3, 1, 1, 0}, eri4KernelFun<3, 1, 1, 0>},
        {{3, 1, 1, 1}, eri4KernelFun<3, 1, 1, 1>},
        {{3, 1, 2, 0}, eri4KernelFun<3, 1, 2, 0>},
        {{3, 2, 0, 0}, eri4KernelFun<3, 2, 0, 0>},
        {{3, 2, 1, 0}, eri4KernelFun<3, 2, 1, 0>},
        {{3, 3, 0, 0}, eri4KernelFun<3, 3, 0, 0>},
        {{4, 0, 0, 0}, eri4KernelFun<4, 0, 0, 0>},
        {{4, 0, 1, 0}, eri4KernelFun<4, 0, 1, 0>},
        {{4, 0, 1, 1}, eri4KernelFun<4, 0, 1, 1>},
        {{4, 0, 2, 0}, eri4KernelFun<4, 0, 2, 0>},
        {{4, 1, 0, 0}, eri4KernelFun<4, 1, 0, 0>},
        {{4, 1, 1, 0}, eri4KernelFun<4, 1, 1, 0>},
        {{4, 2, 0, 0}, eri4KernelFun<4, 2, 0, 0>},
        {{5, 0, 0, 0}, eri4KernelFun<5, 0, 0, 0>},
        {{5, 0, 1, 0}, eri4KernelFun<5, 0, 1, 0>},
        {{5, 1, 0, 0}, eri4KernelFun<5, 1, 0, 0>},
        {{6, 0, 0, 0}, eri4KernelFun<6, 0, 0, 0>}};

    const std::map<std::tuple<int, int, int>, eri3_kernelfun_t> eri3_kernelfuns{
        {{0, 0, 0}, eri3KernelFun<0, 0, 0>},
        {{0, 0, 1}, eri3KernelFun<0, 0, 1>},
        {{0, 0, 2}, eri3KernelFun<0, 0, 2>},
        {{0, 0, 3}, eri3KernelFun<0, 0, 3>},
        {{0, 0, 4}, eri3KernelFun<0, 0, 4>},
        {{0, 0, 5}, eri3KernelFun<0, 0, 5>},
        {{0, 0, 6}, eri3KernelFun<0, 0, 6>},
        {{1, 0, 0}, eri3KernelFun<1, 0, 0>},
        {{1, 0, 1}, eri3KernelFun<1, 0, 1>},
        {{1, 0, 2}, eri3KernelFun<1, 0, 2>},
        {{1, 0, 3}, eri3KernelFun<1, 0, 3>},
        {{1, 0, 4}, eri3KernelFun<1, 0, 4>},
        {{1, 0, 5}, eri3KernelFun<1, 0, 5>},
        {{1, 1, 0}, eri3KernelFun<1, 1, 0>},
        {{1, 1, 1}, eri3KernelFun<1, 1, 1>},
        {{1, 1, 2}, eri3KernelFun<1, 1, 2>},
        {{1, 1, 3}, eri3KernelFun<1, 1, 3>},
        {{1, 1, 4}, eri3KernelFun<1, 1, 4>},
        {{2, 0, 0}, eri3KernelFun<2, 0, 0>},
        {{2, 0, 1}, eri3KernelFun<2, 0, 1>},
        {{2, 0, 2}, eri3KernelFun<2, 0, 2>},
        {{2, 0, 3}, eri3KernelFun<2, 0, 3>},
        {{2, 0, 4}, eri3KernelFun<2, 0, 4>},
        {{2, 1, 0}, eri3KernelFun<2, 1, 0>},
        {{2, 1, 1}, eri3KernelFun<2, 1, 1>},
        {{2, 1, 2}, eri3KernelFun<2, 1, 2>},
        {{2, 1, 3}, eri3KernelFun<2, 1, 3>},
        {{2, 2, 0}, eri3KernelFun<2, 2, 0>},
        {{2, 2, 1}, eri3KernelFun<2, 2, 1>},
        {{2, 2, 2}, eri3KernelFun<2, 2, 2>},
        {{3, 0, 0}, eri3KernelFun<3, 0, 0>},
        {{3, 0, 1}, eri3KernelFun<3, 0, 1>},
        {{3, 0, 2}, eri3KernelFun<3, 0, 2>},
        {{3, 0, 3}, eri3KernelFun<3, 0, 3>},
        {{3, 1, 0}, eri3KernelFun<3, 1, 0>},
        {{3, 1, 1}, eri3KernelFun<3, 1, 1>},
        {{3, 1, 2}, eri3KernelFun<3, 1, 2>},
        {{3, 2, 0}, eri3KernelFun<3, 2, 0>},
        {{3, 2, 1}, eri3KernelFun<3, 2, 1>},
        {{3, 3, 0}, eri3KernelFun<3, 3, 0>},
        {{4, 0, 0}, eri3KernelFun<4, 0, 0>},
        {{4, 0, 1}, eri3KernelFun<4, 0, 1>},
        {{4, 0, 2}, eri3KernelFun<4, 0, 2>},
        {{4, 1, 0}, eri3KernelFun<4, 1, 0>},
        {{4, 1, 1}, eri3KernelFun<4, 1, 1>},
        {{4, 2, 0}, eri3KernelFun<4, 2, 0>},
        {{5, 0, 0}, eri3KernelFun<5, 0, 0>},
        {{5, 0, 1}, eri3KernelFun<5, 0, 1>},
        {{5, 1, 0}, eri3KernelFun<5, 1, 0>},
        {{6, 0, 0}, eri3KernelFun<6, 0, 0>}};

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

LI::ERI4Kernel::ERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                           const eri4_kernelfun_t &eri4_kernelfun)
    : eri4_kernelfun(eri4_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARK(sp_data_ab);
    ecoeffs_ket = ecoeffsSHARK(sp_data_cd);

    int labcd = sp_data_ab.la + sp_data_ab.lb + sp_data_cd.la + sp_data_cd.lb;
    boys_grid = BoysGrid(labcd);
}

LI::ERI3Kernel::ERI3Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                           const eri3_kernelfun_t &eri3_kernelfun)
    : eri3_kernelfun(eri3_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARK(sp_data_ab);
    ecoeffs_ket = ecoeffsSHARK(sh_data_c);

    int labc = sp_data_ab.la + sp_data_ab.lb + sh_data_c.l;
    boys_grid = BoysGrid(labc);
}

LI::ERI2Kernel::ERI2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b,
                           const eri2_kernelfun_t &eri2_kernelfun)
    : eri2_kernelfun(eri2_kernelfun)
{
    ecoeffs_bra = ecoeffsSHARK(sh_data_a);
    ecoeffs_ket = ecoeffsSHARK(sh_data_b);

    int lab = sh_data_a.l + sh_data_b.l;
    boys_grid = BoysGrid(lab);    
}

LI::ERI4Kernel LI::deployERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd)
{
    int la = sp_data_ab.la;
    int lb = sp_data_ab.lb;
    int lc = sp_data_cd.la;
    int ld = sp_data_cd.lb;

    int labcd = la + lb + lc + ld;
    int l_max = _eri_kernel_max_l_;
    if (labcd > _eri_kernel_max_l_)
    {
        string msg = std::format("deployERI4Kernel(): lab + lcd = {} is larger than the allowed max: {}!\n",
                                 labcd, l_max);
        throw std::runtime_error(msg);
    }

    if (labcd <= 6)
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
    int l_max = _eri_kernel_max_l_;
    if (labc > _eri_kernel_max_l_)
    {
        string msg = std::format("deployERI3Kernel(): lab + lc = {} is larger than the allowed max: {}!\n",
                                 labc, l_max);
        throw std::runtime_error(msg);
    }

    if (labc <= 6)
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
    int l_max = _eri_kernel_max_l_;
    if (lab > _eri_kernel_max_l_)
    {
        string msg = std::format("deployERI2Kernel(): la + lb = {} is larger than the allowed max: {}!\n",
                                 lab, l_max);
        throw std::runtime_error(msg);
    }

    if (lab <= 6)
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
