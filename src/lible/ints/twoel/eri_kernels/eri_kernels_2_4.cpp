#include <lible/ints/twoel/eri_kernel_funs.hpp>

template lible::vec4d lible::ints::eri4KernelFun<1, 1, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 1, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<1, 1, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<1, 1, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 1, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 1, 1, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<1, 1, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<1, 1, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 1, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 1, 0, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<1, 1, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 0, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 0, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 2, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 0, 2, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 0, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 0, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 2, 1, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 2, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 0, 1, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 0, 3, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 0, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 0, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 2, 0, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 2, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 0, 0, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 0, 4, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec3d lible::ints::eri3KernelFun<1, 1, 4>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<1, 1, 4>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<1, 1, 4>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<2, 0, 4>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<2, 0, 4>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<0, 2, 4>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<2, 0, 4>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec2d lible::ints::eri2KernelFun<2, 4>(const int ishell_a, const int ishell_b,
                                                       const ShellData &sh_data_a,
                                                       const ShellData &sh_data_b,
                                                       const ERI2Kernel *eri2_kernel);

template std::array<lible::vec2d, 6> lible::ints::eri2d1KernelFun<2, 4>(const int ishell_a, const int ishell_b,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b,
                                                                        const ERI2D1Kernel *eri2d1_kernel);

template lible::vec2d lible::ints::two::eri2Kernel<2, 4>(const int ishell_a, const int ishell_b,
                                                         const std::vector<double> &ecoeffs_a,
                                                         const std::vector<double> &ecoeffs_b_tsp,
                                                         const ShellData &sh_data_a,
                                                         const ShellData &sh_data_b);

template std::array<lible::vec2d, 6> lible::ints::two::eri2d1Kernel<2, 4>(const int ishell_a, const int ishell_b,
                                                                        const std::vector<double> &ecoeffs_a,
                                                                        const std::vector<double> &ecoeffs_b_tsp,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b);

