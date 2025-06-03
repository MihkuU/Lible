#pragma once

#include <lible/ints/boys_function.hpp>
#include <lible/ints/rints.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <functional>

namespace lible
{
    namespace ints
    {
        struct ERI4Kernel;
        struct ERI3Kernel;
        struct ERI2Kernel;
        struct ERI4D1Kernel;
        struct ERI3D1Kernel;
        struct ERI2D1Kernel;
        struct ERI4D2Kernel;
        struct ERI3D2Kernel;
        struct ERI2D2Kernel;

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

        using eri4d1_kernelfun_t = std::function<std::array<vec4d, 12>(
            const int ipair_ab, const int ipair_cd, const ShellPairData &sp_data_ab,
            const ShellPairData &sp_data_cd, const ERI4D1Kernel *eri4d1_kernel)>;

        using eri3d1_kernelfun_t = std::function<std::array<vec3d, 9>(
            const int ipair_ab, const int ishell_c, const ShellPairData &sp_data_ab,
            const ShellData &sh_data_c, const ERI3D1Kernel *eri3d1_kernel)>;

        using eri2d1_kernelfun_t = std::function<std::array<vec2d, 6>(
            const int ishell_a, const int ishell_b, const ShellData &sh_data_a,
            const ShellData &sh_data_b, const ERI2D1Kernel *eri2d1_kernel)>;

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

            BoysGrid boys_grid;
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

            BoysGrid boys_grid;
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

            BoysGrid boys_grid;
        };

        struct ERI4D1Kernel
        {
            ERI4D1Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd,
                         const eri4d1_kernelfun_t &eri4d1_kernelfun);

            std::array<vec4d, 12> operator()(const int ipair_ab, const int ipair_cd,
                                             const ShellPairData &sp_data_ab,
                                             const ShellPairData &sp_data_cd) const
            {
                return eri4d1_kernelfun(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd, this);
            }

            std::vector<double> ecoeffs0_bra;
            std::vector<double> ecoeffs1_bra;
            std::vector<double> ecoeffs0_ket;
            std::vector<double> ecoeffs1_ket;
            eri4d1_kernelfun_t eri4d1_kernelfun;
        };

        struct ERI3D1Kernel
        {
            ERI3D1Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c,
                         const eri3d1_kernelfun_t &eri3d1_kernelfun);

            std::array<vec3d, 9> operator()(const int ipair_ab, const int ishell_c,
                                            const ShellPairData &sp_data_ab,
                                            const ShellData &sh_data_c) const
            {
                return eri3d1_kernelfun(ipair_ab, ishell_c, sp_data_ab, sh_data_c, this);
            }

            std::vector<double> ecoeffs0_bra;
            std::vector<double> ecoeffs1_bra;
            std::vector<double> ecoeffs0_ket;
            eri3d1_kernelfun_t eri3d1_kernelfun;

            BoysGrid boys_grid;
        };

        struct ERI2D1Kernel
        {
            ERI2D1Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b,
                         const eri2d1_kernelfun_t &eri2d1_kernelfun);

            std::array<vec2d, 6> operator()(const int ishell_a, const int ishell_b,
                                            const ShellData &sh_data_a,
                                            const ShellData &sh_data_b) const
            {
                return eri2d1_kernelfun(ishell_a, ishell_b, sh_data_a, sh_data_b, this);
            }

            std::vector<double> ecoeffs_bra;
            std::vector<double> ecoeffs_ket;
            eri2d1_kernelfun_t eri2d1_kernelfun;

            BoysGrid boys_grid;
        };

        ERI4Kernel deployERI4Kernel(const ShellPairData &sp_data_ab,
                                    const ShellPairData &sp_data_cd);

        ERI3Kernel deployERI3Kernel(const ShellPairData &sp_data_ab,
                                    const ShellData &sh_data_c);

        ERI2Kernel deployERI2Kernel(const ShellData &sh_data_a,
                                    const ShellData &sh_data_b);

        ERI4D1Kernel deployERI4D1Kernel(const ShellPairData &sp_data_ab,
                                        const ShellPairData &sp_data_cd);

        ERI3D1Kernel deployERI3D1Kernel(const ShellPairData &sp_data_ab,
                                        const ShellData &sh_data_c);

        ERI2D1Kernel deployERI2D1Kernel(const ShellData &sh_data_a,
                                        const ShellData &sh_data_b);

    }
}