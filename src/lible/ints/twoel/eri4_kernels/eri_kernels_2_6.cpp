#include <lible/ints/twoel/eri_kernels.hpp>

template void lible::ints::two::eri4Kernel<1, 1, 3, 3>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<1, 1, 4, 2>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<1, 1, 5, 1>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<1, 1, 6, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<2, 0, 3, 3>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<2, 0, 4, 2>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<2, 0, 5, 1>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<2, 0, 6, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri3Kernel<1, 1, 6>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<2, 0, 6>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri2Kernel<2, 6>(const int, const int,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 double*);

