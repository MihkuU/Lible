#pragma once

#include <lible/ints/boys_function.hpp>
#include <lible/ints/shell_pair_data.hpp>

#include <functional>

namespace lible::ints
{
    struct ERI4Kernel;
    struct ERI3Kernel;
    struct ERI2Kernel;

    struct ERI4D1Kernel;
    struct ERI3D1Kernel;
    struct ERI2D1Kernel;

    struct ERI2D2Kernel;

    struct ERI4SOCKernel;
    struct ERI3SOCKernel;

    using eri4_kernelfun_t = std::function<vec4d(
        size_t ipair_ab, size_t ipair_cd, const ShellPairData &sp_data_ab,
        const ShellPairData &sp_data_cd, const ERI4Kernel *eri4_kernel)>;

    using eri3_kernelfun_t = std::function<vec3d(
        size_t ipair_ab, size_t ishell_c, const ShellPairData &sp_data_ab,
        const ShellData &sh_data_c, const ERI3Kernel *eri3_kernel)>;

    using eri2_kernelfun_t = std::function<vec2d(
        size_t ishell_a, size_t ishell_b, const ShellData &sp_data_a,
        const ShellData &sh_data_b, const ERI2Kernel *eri2_kernel)>;

    using eri4d1_kernelfun_t = std::function<std::array<vec4d, 12>(
        size_t ipair_ab, size_t ipair_cd, const ShellPairData &sp_data_ab,
        const ShellPairData &sp_data_cd, const ERI4D1Kernel *eri4d1_kernel)>;

    using eri3d1_kernelfun_t = std::function<std::array<vec3d, 9>(
        size_t ipair_ab, size_t ishell_c, const ShellPairData &sp_data_ab,
        const ShellData &sh_data_c, const ERI3D1Kernel *eri3d1_kernel)>;

    using eri2d1_kernelfun_t = std::function<std::array<vec2d, 6>(
        size_t ishell_a, size_t ishell_b, const ShellData &sh_data_a,
        const ShellData &sh_data_b, const ERI2D1Kernel *eri2d1_kernel)>;

    using eri2d2_kernelfun_t = std::function<arr2d<vec2d, 6, 6>(
        size_t ishell_a, size_t ishell_b, const ShellData &sh_data_a,
        const ShellData &sh_data_b, const ERI2D2Kernel *eri2d2_kernel)>;

    using eri4soc_kernelfun_t = std::function<std::array<vec4d, 3>(
        size_t ipair_ab, size_t ipair_cd, const ShellPairData &sp_data_ab,
        const ShellPairData &sp_data_cd, const ERI4SOCKernel *eri4soc_kernel)>;

    using eri3soc_kernelfun_t = std::function<std::array<vec3d, 3>(
        size_t ipair_ab, size_t ishell_c, const ShellPairData &sp_data_ab,
        const ShellData &sh_data_c, const ERI3SOCKernel *eri3soc_kernel)>;

    struct ERI4Kernel
    {
        ERI4Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd);

        vec4d operator()(const size_t ipair_ab, const size_t ipair_cd,
                         const ShellPairData &sp_data_ab,
                         const ShellPairData &sp_data_cd) const
        {
            return eri4_kernelfun_(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd, this);
        }

        std::vector<double> ecoeffs_bra_;
        std::vector<double> ecoeffs_ket_;
        eri4_kernelfun_t eri4_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI3Kernel
    {
        ERI3Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c);

        vec3d operator()(const size_t ipair_ab, const size_t ishell_c,
                         const ShellPairData &sp_data_ab,
                         const ShellData &sh_data_c) const
        {
            return eri3_kernelfun_(ipair_ab, ishell_c, sp_data_ab, sh_data_c, this);
        }

        std::vector<double> ecoeffs_bra_;
        std::vector<double> ecoeffs_ket_;
        eri3_kernelfun_t eri3_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI2Kernel
    {
        ERI2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b);

        vec2d operator()(const size_t ishell_a, const size_t ishell_b,
                         const ShellData &sh_data_a,
                         const ShellData &sh_data_b) const
        {
            return eri2_kernelfun_(ishell_a, ishell_b, sh_data_a, sh_data_b, this);
        }

        std::vector<double> ecoeffs_bra_;
        std::vector<double> ecoeffs_ket_;
        eri2_kernelfun_t eri2_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI4D1Kernel
    {
        ERI4D1Kernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd);

        std::array<vec4d, 12> operator()(const size_t ipair_ab, const size_t ipair_cd,
                                         const ShellPairData &sp_data_ab,
                                         const ShellPairData &sp_data_cd) const
        {
            return eri4d1_kernelfun_(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd, this);
        }

        std::vector<double> ecoeffs0_bra_;
        std::vector<double> ecoeffs1_bra_;
        std::vector<double> ecoeffs0_ket_;
        std::vector<double> ecoeffs1_ket_;
        eri4d1_kernelfun_t eri4d1_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI3D1Kernel
    {
        ERI3D1Kernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c);

        std::array<vec3d, 9> operator()(const size_t ipair_ab, const size_t ishell_c,
                                        const ShellPairData &sp_data_ab,
                                        const ShellData &sh_data_c) const
        {
            return eri3d1_kernelfun_(ipair_ab, ishell_c, sp_data_ab, sh_data_c, this);
        }

        std::vector<double> ecoeffs0_bra_;
        std::vector<double> ecoeffs1_bra_;
        std::vector<double> ecoeffs0_ket_;
        eri3d1_kernelfun_t eri3d1_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI2D1Kernel
    {
        ERI2D1Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b);

        std::array<vec2d, 6> operator()(const size_t ishell_a, const size_t ishell_b,
                                        const ShellData &sh_data_a,
                                        const ShellData &sh_data_b) const
        {
            return eri2d1_kernelfun_(ishell_a, ishell_b, sh_data_a, sh_data_b, this);
        }

        std::vector<double> ecoeffs_bra_;
        std::vector<double> ecoeffs_ket_;
        eri2d1_kernelfun_t eri2d1_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI2D2Kernel
    {
        ERI2D2Kernel(const ShellData &sh_data_a, const ShellData &sh_data_b);

        arr2d<vec2d, 6, 6> operator()(const size_t ishell_a, const size_t ishell_b,
                                      const ShellData &sh_data_a,
                                      const ShellData &sh_data_b) const
        {
            return eri2d2_kernelfun_(ishell_a, ishell_b, sh_data_a, sh_data_b, this);
        }

        std::vector<double> ecoeffs_bra_;
        std::vector<double> ecoeffs_ket_;
        eri2d2_kernelfun_t eri2d2_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI4SOCKernel
    {
        ERI4SOCKernel(const ShellPairData &sp_data_ab, const ShellPairData &sp_data_cd);

        std::array<vec4d, 3> operator()(const size_t ipair_ab, const size_t ipair_cd,
                                        const ShellPairData &sp_data_ab,
                                        const ShellPairData &sp_data_cd) const
        {
            return eri4soc_kernelfun_(ipair_ab, ipair_cd, sp_data_ab, sp_data_cd, this);
        }

        std::vector<double> ecoeffs1_bra_;
        std::vector<double> ecoeffs0_ket_;
        eri4soc_kernelfun_t eri4soc_kernelfun_;

        BoysGrid boys_grid_;
    };

    struct ERI3SOCKernel
    {
        ERI3SOCKernel(const ShellPairData &sp_data_ab, const ShellData &sh_data_c);

        std::array<vec3d, 3> operator()(const size_t ipair_ab, const size_t ishell_c,
                                        const ShellPairData &sp_data_ab,
                                        const ShellData &sh_data_c) const
        {
            return eri3soc_kernelfun_(ipair_ab, ishell_c, sp_data_ab, sh_data_c, this);
        }

        std::vector<double> ecoeffs1_bra_;
        std::vector<double> ecoeffs0_ket_;
        eri3soc_kernelfun_t eri3soc_kernelfun_;

        BoysGrid boys_grid_;
    };
}
