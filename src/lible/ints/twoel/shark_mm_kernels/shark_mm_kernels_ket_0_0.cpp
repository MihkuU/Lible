#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_ket1<0, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
}

template<> void lible::ints::shark_mm_ket2<0, 0, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
}

