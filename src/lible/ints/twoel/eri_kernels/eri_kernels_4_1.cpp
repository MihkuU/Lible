#include <lible/ints/twoel/eri_kernels.hpp>

template void lible::ints::two::eri4Kernel<2, 2, 1, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<3, 1, 1, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri4Kernel<4, 0, 1, 0>(const int, const int, const int, const int,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       const double*, const double*,
                                                       double*);

template void lible::ints::two::eri3Kernel<2, 2, 1>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<3, 1, 1>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri3Kernel<4, 0, 1>(const int, const int, const int,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, const double*,
                                                    const double*, const double*, double*);

template void lible::ints::two::eri2Kernel<4, 1>(const int, const int,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 const double*, const double*,
                                                 double*);

