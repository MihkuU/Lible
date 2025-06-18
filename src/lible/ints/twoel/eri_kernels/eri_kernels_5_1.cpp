#include <lible/ints/twoel/eri_kernel_funs.hpp>

template lible::vec4d lible::ints::eri4KernelFun<3, 2, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 2, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 3, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 3, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 2, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::eri4KernelFun<2, 3, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<2, 3, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<3, 2, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<4, 1, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 1, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 4, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 4, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 1, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::eri4KernelFun<1, 4, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<1, 4, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<4, 1, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<5, 0, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 0, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 5, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 5, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 0, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::eri4KernelFun<0, 5, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<0, 5, 1, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec4d lible::ints::eri4KernelFun<5, 0, 0, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template lible::vec3d lible::ints::eri3KernelFun<3, 2, 1>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<3, 2, 1>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<2, 3, 1>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::eri3KernelFun<4, 1, 1>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<4, 1, 1>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<1, 4, 1>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::eri3KernelFun<5, 0, 1>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<5, 0, 1>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<0, 5, 1>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec2d lible::ints::eri2KernelFun<5, 1>(const int ishell_a, const int ishell_b,
                                                       const ShellData &sh_data_a,
                                                       const ShellData &sh_data_b,
                                                       const ERI2Kernel *eri2_kernel);

template std::array<lible::vec2d, 6> lible::ints::eri2d1KernelFun<5, 1>(const int ishell_a, const int ishell_b,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b,
                                                                        const ERI2D1Kernel *eri2d1_kernel);

