#include <lible/ints/twoel/shark_mm_kernels.hpp>

template<> void lible::ints::shark_mm_ket1<1, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[1] += R[1] * ET[0];
    R_x_ET[2] += R[2] * ET[0];
    R_x_ET[3] += R[3] * ET[0];
}

template<> void lible::ints::shark_mm_ket2<1, 0, 0>(const double *R, const double *ET, double *R_x_ET)
{
    R_x_ET[0] += R[0] * ET[0];
    R_x_ET[1] += R[1] * ET[0];
    R_x_ET[2] += R[2] * ET[0];
    R_x_ET[3] += R[3] * ET[0];
}

