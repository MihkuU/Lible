#include <lible/ints/twoel/eri_kernels.hpp>

template void lible::ints::two::eri4Kernel<1, 1, 2, 1>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<1, 1, 3, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<2, 0, 2, 1>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<2, 0, 3, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri3Kernel<1, 1, 3>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<2, 0, 3>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri2Kernel<2, 3>(const int, const int,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 double*);

