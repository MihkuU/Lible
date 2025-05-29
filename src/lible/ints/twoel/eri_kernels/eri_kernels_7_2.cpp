#include <lible/ints/twoel/eri_kernel_funs.hpp>

template lible::vec4d lible::ints::eri4KernelFun<4, 3, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 3, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<4, 3, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 3, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<5, 2, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<5, 2, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<5, 2, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<5, 2, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 1, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 1, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 1, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 1, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 0, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 0, 1, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 0, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 0, 2, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec3d lible::ints::eri3KernelFun<4, 3, 2>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<4, 3, 2>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<5, 2, 2>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<5, 2, 2>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<6, 1, 2>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<6, 1, 2>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<7, 0, 2>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<7, 0, 2>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec2d lible::ints::eri2KernelFun<7, 2>(const int ishell_a, const int ishell_b,
                                                       const ShellData &sh_data_a,
                                                       const ShellData &sh_data_b,
                                                       const ERI2Kernel *eri2_kernel);

template lible::vec2d lible::ints::two::eri2Kernel<7, 2>(const int ishell_a, const int ishell_b,
                                                         const std::vector<double> &ecoeffs_a,
                                                         const std::vector<double> &ecoeffs_b_tsp,
                                                         const ShellData &sh_data_a,
                                                         const ShellData &sh_data_b);

template std::array<lible::vec2d, 6> lible::ints::two::eri2d1Kernel<7, 2>(const int ishell_a, const int ishell_b,
                                                                        const std::vector<double> &ecoeffs_a,
                                                                        const std::vector<double> &ecoeffs_b_tsp,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b);

