#include <lible/ints/twoel/eri_kernels.hpp>

template void lible::ints::two::eri4Kernel<5, 5, 0, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<6, 4, 0, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<7, 3, 0, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<8, 2, 0, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<9, 1, 0, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<10, 0, 0, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri3Kernel<5, 5, 0>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<6, 4, 0>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<7, 3, 0>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<8, 2, 0>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<9, 1, 0>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<10, 0, 0>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri2Kernel<10, 0>(const int, const int,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 double*);
