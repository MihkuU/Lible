#pragma once

#include <lible/ints/shell_pair_data.hpp>

#include <functional>

namespace lible
{
    namespace ints
    {
        struct ERI4Kernel;
        struct ERI3Kernel;
        struct ERI2Kernel;
        struct ERI4KernelD1;
        struct ERI3KernelD1;
        struct ERI2KernelD1;
        struct ERI4KernelD2;
        struct ERI3KernelD2;
        struct ERI2KernelD2;        

        using eri4_kernelfun_t = std::function<vec4d(const int ipair_ab, const int ipair_cd,
                                                     const ShellPairData &sp_data_ab,
                                                     const ShellPairData &sp_data_cd,
                                                     const ERI4Kernel *eri4_kernel)>;

        using eri3_kernelfun_t = std::function<vec3d(const int ipair_ab, const int ishell_c,
                                                     const ShellPairData &sp_data_ab,
                                                     const ShellData &sh_data_c,
                                                     const ERI3Kernel *eri3_kernel)>;

        using eri2_kernelfun_t = std::function<vec2d(const int ishell_a, const int ishell_b,
                                                     const ShellData &sp_data_a, 
                                                     const ShellData &sh_data_b,
                                                     const ERI2Kernel *eri2_kernel)>;

        struct ERI4Kernel
        {
            ERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                       const eri4_kernelfun_t &eri4_kernelfun);

            vec4d operator()(const int ipair_ab, const int ipair_cd,
                             const ShellPairData &sp_data_ab,
                             const ShellPairData &sp_data_cd) const
            {
                return eri4_kernelfun(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd, this);
            }

            std::vector<double> ecoeffs_bra;
            std::vector<double> ecoeffs_ket;
            eri4_kernelfun_t eri4_kernelfun;
        };

        struct ERI3Kernel
        {
            ERI3Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                       const eri3_kernelfun_t &eri3_kernelfun);

            vec3d operator()(const int ipair_ab, const int ishell_c,
                             const ShellPairData &sp_data_ab,
                             const ShellData &sh_data_c) const
            {
                return eri3_kernelfun(ipair_ab, ishell_c, sp_data_ab, sh_data_c, this);
            }

            std::vector<double> ecoeffs_bra;
            std::vector<double> ecoeffs_ket;
            eri3_kernelfun_t eri3_kernelfun;
        };

        struct ERI2Kernel
        {
            ERI2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b,
                       const eri2_kernelfun_t &eri2_kernelfun);

            vec2d operator()(const int ishell_a, const int ishell_b, 
                             const ShellData &sh_data_a,
                             const ShellData &sh_data_b) const
            {
                return eri2_kernelfun(ishell_a, ishell_b, sh_data_a, sh_data_b, this);
            }

            std::vector<double> ecoeffs_bra;
            std::vector<double> ecoeffs_ket;
            eri2_kernelfun_t eri2_kernelfun;
        };

        ERI4Kernel deployERI4Kernel(const ShellPairData &sp_data_ab,
                                    const ShellPairData &sp_data_cd);

        ERI3Kernel deployERI3Kernel(const ShellPairData &sp_data_ab,
                                    const ShellData &sh_data_c);

        ERI2Kernel deployERI2Kernel(const ShellData &sh_data_a,
                                    const ShellData &sh_data_b);
    }
}