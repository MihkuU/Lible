#include <lible/ints/twoel/eri_kernel_funs.hpp>

template lible::vec4d lible::ints::eri4KernelFun<6, 6, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 6, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 6, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 6, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 6, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 6, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 6, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 6, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 6, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 6, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 6, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 6, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<6, 6, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<6, 6, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<6, 6, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 5, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 5, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 5, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 5, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 5, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 5, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 5, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 5, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 5, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 5, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 5, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 5, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<7, 5, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<5, 7, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<7, 5, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<7, 5, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<8, 4, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<8, 4, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<8, 4, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<8, 4, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<8, 4, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<8, 4, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<8, 4, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<8, 4, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<8, 4, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<8, 4, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<8, 4, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<8, 4, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<8, 4, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 8, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<8, 4, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<8, 4, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<9, 3, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<9, 3, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<9, 3, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<9, 3, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<9, 3, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<9, 3, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<9, 3, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<9, 3, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<9, 3, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<9, 3, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<9, 3, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<9, 3, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<9, 3, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 9, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<9, 3, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<9, 3, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<10, 2, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<10, 2, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<10, 2, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<10, 2, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<10, 2, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<10, 2, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<10, 2, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<10, 2, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<10, 2, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<10, 2, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<10, 2, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<10, 2, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<10, 2, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 10, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<10, 2, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<10, 2, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<11, 1, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<11, 1, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<11, 1, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<11, 1, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<11, 1, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<11, 1, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<11, 1, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<11, 1, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<11, 1, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<11, 1, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<11, 1, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<11, 1, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<11, 1, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 11, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<11, 1, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<11, 1, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<12, 0, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<12, 0, 6, 6>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<12, 0, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 5, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<12, 0, 7, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<12, 0, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 4, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<12, 0, 8, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<12, 0, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 3, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<12, 0, 9, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<12, 0, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 2, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<12, 0, 10, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<12, 0, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 1, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<12, 0, 11, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<12, 0, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 12, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<12, 0, 0, 12>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<12, 0, 12, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec3d lible::ints::eri3KernelFun<6, 6, 12>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<6, 6, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<6, 6, 12>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<7, 5, 12>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<7, 5, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<5, 7, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<7, 5, 12>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<8, 4, 12>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<8, 4, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<4, 8, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<8, 4, 12>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<9, 3, 12>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<9, 3, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<3, 9, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<9, 3, 12>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<10, 2, 12>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<10, 2, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<2, 10, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<10, 2, 12>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<11, 1, 12>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<11, 1, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<1, 11, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<11, 1, 12>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<12, 0, 12>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<12, 0, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<0, 12, 12>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<12, 0, 12>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec2d lible::ints::eri2KernelFun<12, 12>(const int ishell_a, const int ishell_b,
                                                       const ShellData &sh_data_a,
                                                       const ShellData &sh_data_b,
                                                       const ERI2Kernel *eri2_kernel);

template std::array<lible::vec2d, 6> lible::ints::eri2d1KernelFun<12, 12>(const int ishell_a, const int ishell_b,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b,
                                                                        const ERI2D1Kernel *eri2d1_kernel);

template lible::vec2d lible::ints::two::eri2Kernel<12, 12>(const int ishell_a, const int ishell_b,
                                                         const std::vector<double> &ecoeffs_a,
                                                         const std::vector<double> &ecoeffs_b_tsp,
                                                         const ShellData &sh_data_a,
                                                         const ShellData &sh_data_b);

template std::array<lible::vec2d, 6> lible::ints::two::eri2d1Kernel<12, 12>(const int ishell_a, const int ishell_b,
                                                                        const std::vector<double> &ecoeffs_a,
                                                                        const std::vector<double> &ecoeffs_b_tsp,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b);

