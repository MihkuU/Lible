#include <lible/ints/twoel/eri_kernel_funs.hpp>

template lible::vec4d lible::ints::eri4KernelFun<2, 2, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 5, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 2, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 2, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 4, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 2, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 2, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 3, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 2, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 2, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 2, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 2, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 2, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 1, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 2, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<2, 2, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<2, 2, 0, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<2, 2, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<3, 1, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 5, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 5, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<3, 1, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<3, 1, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 4, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 4, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<3, 1, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<3, 1, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 3, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 3, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<3, 1, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<3, 1, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 2, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 2, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<3, 1, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<3, 1, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 1, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 1, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<3, 1, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<3, 1, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 0, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<1, 3, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<3, 1, 0, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<3, 1, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<4, 0, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 5, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 5, 6>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 0, 6, 5>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<4, 0, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 4, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 4, 7>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 0, 7, 4>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<4, 0, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 3, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 3, 8>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 0, 8, 3>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<4, 0, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 2, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 2, 9>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 0, 9, 2>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<4, 0, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 1, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 1, 10>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 0, 10, 1>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec4d lible::ints::eri4KernelFun<4, 0, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellPairData &sp_data_cd,
                                                             const ERI4Kernel *eri4_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 0, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<0, 4, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template std::array<lible::vec4d, 12> lible::ints::eri4d1KernelFun<4, 0, 0, 11>(const int ipair_ab, const int ipair_cd,
                                                                               const ShellPairData &sh_data_ab,
                                                                               const ShellPairData &sp_data_cd,
                                                                               const ERI4D1Kernel *eri4d1_kernel);

template lible::vec4d lible::ints::two::eri4Kernel<4, 0, 11, 0>(const int ipair_ab, const int ipair_cd,
                                                               const std::vector<double> &ecoeffs_ab,
                                                               const std::vector<double> &ecoeffs_cd_tsp,
                                                               const ShellPairData &sp_data_ab,
                                                               const ShellPairData &sp_data_cd);

template lible::vec3d lible::ints::eri3KernelFun<2, 2, 11>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<2, 2, 11>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<2, 2, 11>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<3, 1, 11>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<3, 1, 11>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<1, 3, 11>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<3, 1, 11>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec3d lible::ints::eri3KernelFun<4, 0, 11>(const int ipair_ab, const int ishell_c,
                                                          const ShellPairData &sp_data_ab,
                                                          const ShellData &sh_data_c,
                                                          const ERI3Kernel *eri3_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<4, 0, 11>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template std::array<lible::vec3d, 9> lible::ints::eri3d1KernelFun<0, 4, 11>(const int ipair_ab, const int ishell_c,
                                                                           const ShellPairData &sh_data_ab,
                                                                           const ShellData &sh_data_c,
                                                                           const ERI3D1Kernel *eri3d1_kernel);

template lible::vec3d lible::ints::two::eri3Kernel<4, 0, 11>(const int ipair_ab, const int ishell_c,
                                                             const std::vector<double> &ecoeffs_ab,
                                                             const std::vector<double> &ecoeffs_c,
                                                             const ShellPairData &sp_data_ab,
                                                             const ShellData &sh_data_c);

template lible::vec2d lible::ints::eri2KernelFun<4, 11>(const int ishell_a, const int ishell_b,
                                                       const ShellData &sh_data_a,
                                                       const ShellData &sh_data_b,
                                                       const ERI2Kernel *eri2_kernel);

template std::array<lible::vec2d, 6> lible::ints::eri2d1KernelFun<4, 11>(const int ishell_a, const int ishell_b,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b,
                                                                        const ERI2D1Kernel *eri2d1_kernel);

template lible::vec2d lible::ints::two::eri2Kernel<4, 11>(const int ishell_a, const int ishell_b,
                                                         const std::vector<double> &ecoeffs_a,
                                                         const std::vector<double> &ecoeffs_b_tsp,
                                                         const ShellData &sh_data_a,
                                                         const ShellData &sh_data_b);

template std::array<lible::vec2d, 6> lible::ints::two::eri2d1Kernel<4, 11>(const int ishell_a, const int ishell_b,
                                                                        const std::vector<double> &ecoeffs_a,
                                                                        const std::vector<double> &ecoeffs_b_tsp,
                                                                        const ShellData &sh_data_a,
                                                                        const ShellData &sh_data_b);

